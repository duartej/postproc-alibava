// ***********************************************
// 
// Pool of functions to post-process alibava data
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "AlibavaPostProcessor.h"

// The alibava run and event headers
#include "AuxiliaryStructures.h"
#include "ALIBAVA.h"

// ROOT 
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"

// System headers
#include <fstream>
#include <map>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace auxmem
{
    // deallocating  helper
    void deallocate_memory(std::map<int,std::vector<float>*> & m)
    {
        for(auto & i_v: m)
        {
            if(i_v.second != nullptr)
            {
                delete i_v.second;
                i_v.second = nullptr;
            }
        }
    }
    template <class T>
        void deallocate_memory(std::map<int,std::map<int,T*> > & histos)
        {
            for(auto & i_m2: histos)
            {
                for(auto & i_mI: i_m2.second)
                {
                    if(i_mI.second != nullptr)
                    {
                        delete i_mI.second;
                        i_mI.second = nullptr;
                    }
                }
            }
        }
}

AlibavaPostProcessor::AlibavaPostProcessor():
    _pedestal_subtracted(false),
    _postproc_treename("postproc_Events"),
    _channel_mask({{0,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)},
            {1,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)}}),
    _automasking(true),
    _mask_criterium(2.5)
{
}

AlibavaPostProcessor::AlibavaPostProcessor(const std::map<int,std::vector<int> > & channel_mask):
    _pedestal_subtracted(false),
    _postproc_treename("postproc_Events"),
    _channel_mask({{0,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)},
            {1,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)}}),
    _automasking(true),
    _mask_criterium(2.5)
{
    // All channels initialize to 0 (masked)    
    for(auto & chip_vect: channel_mask)
    {
        // Valid only chip 0 and 1
        if(chip_vect.first != 0 && chip_vect.first != 1)
        {
            std::cerr << "[AlibavaPostProcessor::AlibavaPostProcessor ERROR] Invalid chip"
                << " number [" << chip_vect.first << "] " << std::endl;
            // Exception??
            return;
        }
        // Be sure the user introduced the vector, then mask all channels
        // so activate them only if requested
        if(chip_vect.second.size() > 0)
        {
            std::fill(_channel_mask[chip_vect.first].begin(),_channel_mask[chip_vect.first].end(),0);
        }
        for(auto & ch: chip_vect.second)
        {
            if(ch >= static_cast<int>(_channel_mask[chip_vect.first].size()))
            {
                std::cerr << "[AlibavaPostProcessor::AlibavaPostProcessor WARNING] Invalid "
                    << " channel number [" << ch << "] " << std::endl;
            }
            _channel_mask[chip_vect.first][ch]=1;
        }
    }
}

// Alternative algorithm using the calibration_charge from the Events tree
CalibrateBeetleMap AlibavaPostProcessor::calibrate(IOManager & gauge)
{
    std::string data_b1_brname("data_beetle1");
    std::string data_b2_brname("data_beetle2");
    
    std::string charge_brname("calibration_charge");
    
    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { data_b1_brname, data_b2_brname, charge_brname };
    gauge.set_events_tree_access(data_names);
    
    // for the charge
    float thecharge = -99999.9;
    gauge.set_events_tree_branch_address(charge_brname,&thecharge);
    // Helper map 
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    gauge.set_events_tree_branch_address(data_b1_brname,&(thedata[0]));
    gauge.set_events_tree_branch_address(data_b2_brname,&(thedata[1]));
    // The histograms (XXX: Hardcoded bins and ranges)
    std::map<int,std::map<int,TProfile*> > histos;
    std::map<int,std::map<int,TF1*> >  cal_curve;
    for(int chip = 0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        std::map<int,TProfile*> hchip;
        std::map<int,TF1*> fchip;
        for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS ; ++ichan)
        {
            hchip[ichan] = new TProfile(std::string("calibration"+std::to_string(ichan)).c_str(),"",
                    2*gauge.get_calibration_parameters()->nPulses+1,
                   (-1)*gauge.get_calibration_parameters()->finalCharge,
                   gauge.get_calibration_parameters()->finalCharge,
                   -1000.0,1000.0);
            hchip[ichan]->SetDirectory(0);
            fchip[ichan] = new TF1(std::string("calibration_curve"+std::to_string(ichan)).c_str(),"pol1");
        }
        histos[chip] = hchip;
        cal_curve[chip] = fchip;
    }
    // loop over the tree to fill the histograms
    // filling the histograms injected charge vs. ADC count 
    // per events 
    // NOTE the charge injection pattern follows:
    //
    //     EVENT   |   STRIP    | CHARGE SIGN
    // ------------+------------+-------------
    //  evt%2 == 0 |  i%2 == 0  |  NEGATIVE
    //  evt%2 == 0 |  i%2 == 1  |  POSITIVE
    //  ======================================
    //  evt%2 == 1 |  i%2 == 0  |  POSITIVE
    //  evt%2 == 1 |  i%2 == 1  |  NEGATIVE
    //  
    //  That conditions can be summarize with
    //  (-1)**[ (EVENT%2) + ((STRIP+1)%2) ]
    const int nentries = gauge.get_events_number_entries();
    for(int k = 0; k < nentries; ++k)
    {
        thecharge = -99999.9;
        gauge.get_events_entry(k);
        for(int beetle = 0; beetle < ALIBAVA::NOOFCHIPS; ++beetle)
        {
            // The injected charge is the same for all the strips, 
            // just taking the first one
            const int injected_pulse = static_cast<int>( thecharge );

            for(int istrip = 0; istrip < ALIBAVA::NOOFCHANNELS; ++istrip)
            {
                // Note that the injected charge will depend on event parity
                // and channel parity
                //const int sign = std::pow(-1,((k+1)%2))*std::pow(-1,(istrip%2));
                // or static_cast<bool>(k%2) XOR static_cast<bool>(i%2) 
                // i.e. a XNOR gate 
                const int sign = std::pow(-1,((k%2)+(istrip+1)%2));
                histos[beetle][istrip]->Fill(sign*injected_pulse, (*(thedata[beetle]))[istrip] );
            }
        }
    }
    // The return data
    CalibrateBeetleMap output;
    output[0].reserve(ALIBAVA::NOOFCHANNELS);
    output[1].reserve(ALIBAVA::NOOFCHANNELS);
    // The fit: a canvas for the fits
    // Loop over the channels to fit the calibration curve (stright line):
    //   - the slope of the line: [ADC counts/#electrons], therefore
    //   the number of electrons corresponding to a given ADC counts: 1/slope
    TCanvas * cc = new TCanvas("def");
    // An auxiliary map to check that the fit worked well
    std::map<int,std::vector<int>> fitstatus_m = { {0,{0}},{0,{0}} };
    for(auto & v: fitstatus_m)
    {
        v.second.reserve(ALIBAVA::NOOFCHANNELS);
    }
    for(int chip = 0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
        {
            // the fit
            const int fit_status = histos[chip][ichan]->Fit(cal_curve[chip][ichan],"Q");
            fitstatus_m[chip].push_back(fit_status);
            // Slope: [ADC counts/number of electrons*1e-3] --> 1.0/slope
            //const int elec_per_adc = std::round(1.0/(tempfit->GetParameter(1))); 
            //const float err = tempfit->GetParError(1)*elec_per_adc; 
            output[chip].push_back( static_cast<int>(std::round(1.0/(cal_curve[chip][ichan])->GetParameter(1))) );
        }
    }
    // Check if anything went wrong, 
    for(const auto & chip_v: fitstatus_m)
    {
        for(int ichan = 0; ichan < static_cast<int>(chip_v.second.size()); ++ichan)
        {
            if(chip_v.second[ichan] != 0)
            {
                std::cout << "\033[1;33mWARNING\033[1;m Problems with the "
                    << "line fit @ CHIP: " << chip_v.first 
                    << " , STRIP: " << ichan << " (Fit status:" 
                    << chip_v.second[ichan] << ")" << std::endl;
            }
        }
    }
    // Persistify the calibration curves (note the chip number must start with 1)
    for(const auto & hdictperchip: histos)
    {
        const int chipnumber = hdictperchip.first+1;
        for(const auto & hperchannel: hdictperchip.second)
        {
            const int ichannel = hperchannel.first;
            gauge.book_monitor_plot("calibration_profile_"+std::to_string(ichannel),hperchannel.second,chipnumber);
        }
    }
    // And fill the calibration monitor 3d histogram based on 
    // Deallocating memory
    delete cc;
    cc = nullptr;
    // Deallocate memory from ROOT stuff
    auxmem::deallocate_memory<TProfile>(histos);
    auxmem::deallocate_memory<TF1>(cal_curve);
    // vector and maps
    auxmem::deallocate_memory(thedata);
    // RE-activate all branches (and reset their object addresses)
    gauge.reset_events_tree();

    return output;
}

void AlibavaPostProcessor::mask_channels(const PedestalNoiseBeetleMap & ped_noise)
{
    // Asumming configuration
    if(not _automasking)
    {
        std::cout << "\033[1;33mWARNING\033[1;m Automasking is de-activated " << std::endl;
        // The user decide to put this masking 
        // As this algorithm needs to be run (in order to include the channels masked 
        // by the user), just put a never-filled criterium
        _mask_criterium = -99999.9;
    }
    
    for(int chip=0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        // 1. -- Define a map to keep only non-noisy channels
        // ichannel: noise
        std::map<int,float> non_noisy_map;
        
        // Get the noise vector and check the noisy channel condition
        // which is given by:
        //    i-strip is noisy if |N_{i} - <Noise> > Xsigma
        for(int ch=0; ch < static_cast<int>(ped_noise.at(chip).second.size()); ++ch)
        {
            // Not using already masked channels
            if(this->is_channel_masked(chip,ch))
            {
                continue;
            }
            non_noisy_map.emplace(ch,ped_noise.at(chip).second[ch]);
        }

        // 1. Obtain the mean and the standard deviation for the noise
        float mean_noise  = get_mean(non_noisy_map);
        float stddev_noise= get_std_dev(non_noisy_map,mean_noise);
        
        // 2. Use the mean and std. dev to find noise channels
        //    while using criteria defined by the user, by
        //    removing iteratively noisy channels until converge, i.e
        //    the non_noisy_map will not loose any element anymore
        unsigned int last_vec_size = ped_noise.at(chip).second.size();
        std::cout << "\033[1;34mINFO [AUTO-MASKING]\033[1;m CHIP: " << chip << std::endl;
        do
        {
            // after the check from the 'while' statement, update
            // the last know vector size
            last_vec_size = non_noisy_map.size();
            for(auto it = non_noisy_map.cbegin(); it != non_noisy_map.cend(); /* no increment*/)
            {
                if(this->is_channel_masked(chip,it->first))
                {
                    ++it;
                    continue;
                }
                // Remove the noisy channels and tag it as masked channel
                if( std::abs(it->second-mean_noise)/stddev_noise > _mask_criterium )
                {
                    std::cout << "\033[1;34mINFO [AUTO-MASKING]\033[1;m noisy channel: " 
                        << std::setw(3) << it->first << " with noise: " << std::setprecision(3) 
                        << std::setw(5) << it->second << " (<noise>_{all}: " 
                        << std::setprecision(3) << std::setw(5) << mean_noise 
                        << ", st. dev.: " << std::setprecision(3) << std::setw(5) << stddev_noise << ")" 
                        << std::endl;
                    _channel_mask[chip][it->first] = 0;
                    non_noisy_map.erase(it++);
                }
                else
                {
                    // or keep it
                    ++it;
                }
            }
            // 3. Recalculate the mean with the excluded signal channels
            mean_noise = get_mean(non_noisy_map);
            stddev_noise= get_std_dev(non_noisy_map,mean_noise);
        } while( non_noisy_map.size() != last_vec_size );
    }
}

void AlibavaPostProcessor::print_channels_mask() const
{
    std::string automasking_str("Yes");
    if(!_automasking)
    {
        automasking_str = "No";
    }
    // define columns for the channels
    const int channelprintnum = 16;
    
    std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m Automasking activated : " 
        << std::setw(3) << automasking_str << std::endl;
    std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m Automasking noisy channels if " 
        << "|noise_{i} - <noise>| > " << std::setprecision(2) << std::setw(4) 
        << _mask_criterium << " x sigma_{noise} " << std::endl;
    std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m ******************** Channel Masking ******************* " 
        << std::endl;
    for(int chip=0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m CHIP: " << chip << std::endl;
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ) 
        {
            if(ichan % channelprintnum ==0)
            {
                if(ichan !=0)
                {
                    std::cout << std::endl;
                }
                std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m   Channels "<< std::setw(3) << ichan 
                    << " - "<< std::setw(3) << ichan+channelprintnum-1  << " : ";
            }
            std::cout << !this->is_channel_masked(chip,ichan) << " ";
            ++ichan;
        }
        std::cout << std::endl;
    }
    std::cout << "\033[1;34mINFO [CHANNEL MASK]\033[1;m ******************** Channel Masking ******************* " 
        << std::endl;
}


PedestalNoiseBeetleMap AlibavaPostProcessor::calculate_pedestal_noise(const IOManager & pedestal)
{
    std::string data_b1_brname("data_beetle1");
    std::string data_b2_brname("data_beetle2");
    
    if(this->_pedestal_subtracted)
    {
        data_b1_brname = "postproc_data_beetle1";
        data_b2_brname = "postproc_data_beetle2";
    }

    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { data_b1_brname, data_b2_brname};
    pedestal.set_events_tree_access(data_names);
    
    // Helper map 
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    pedestal.set_events_tree_branch_address(data_b1_brname,&(thedata[0]));
    pedestal.set_events_tree_branch_address(data_b2_brname,&(thedata[1]));
    // The histograms (XXX: Hardcoded bins and ranges)
    std::map<int,std::map<int,TH1F*> > histos;
    std::map<int,std::map<int,TF1*> >  gausfunc;
    for(int chip = 0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        std::map<int,TH1F*> hchip;
        std::map<int,TF1*> fchip;
        for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS ; ++ichan)
        {
            hchip[ichan] = new TH1F(std::string("pedestal_noise"+std::to_string(ichan)).c_str(),"",
                   1000,0.0,1000.0);
            hchip[ichan]->SetDirectory(0);
            fchip[ichan] = new TF1(std::string("pedestal_noise_gaus"+std::to_string(ichan)).c_str(),"gaus");
        }
        histos[chip] = hchip;
        gausfunc[chip] = fchip;
    }

    // XXX: HERE IT SHOULD BE FILLED the Event histogram/TGraph --> PROBLEM
    //      THE CLIENT IS THE PEDESTAL Manager!!

    // loop over the tree to fill the histograms
    const int nentries = pedestal.get_events_number_entries();
    for(int k = 0; k < nentries; ++k)
    {
        pedestal.get_events_entry(k);
        for(int beetle = 0; beetle < ALIBAVA::NOOFCHIPS; ++beetle)
        {
            for(int istrip = 0; istrip < ALIBAVA::NOOFCHANNELS; ++istrip)
            {
                histos[beetle][istrip]->Fill( (*(thedata[beetle]))[istrip] );
            }
        }        
    }
    // The return data
    PedestalNoiseBeetleMap output;
    output[0].first.reserve(ALIBAVA::NOOFCHANNELS);
    output[0].second.reserve(ALIBAVA::NOOFCHANNELS);
    output[1].first.reserve(ALIBAVA::NOOFCHANNELS);
    output[1].second.reserve(ALIBAVA::NOOFCHANNELS);

    // The fit: a canvas for the fits
    // Loop over the channels to fit the gaussian:
    //  - mean: pedestals
    //  - sigma: noise (if mean=0, then the noise is the common 
    //           noise)
    TCanvas * cc = new TCanvas("def");
    // An auxiliary map to check that the fit worked well
    std::map<int,std::vector<int>> fitstatus_m = { {0,{}},{0,{0}} };
    for(auto & v: fitstatus_m)
    {
        v.second.reserve(ALIBAVA::NOOFCHANNELS);
    }
    for(int chip = 0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
        {
            const int fit_status = histos[chip][ichan]->Fit(gausfunc[chip][ichan],"Q");
            fitstatus_m[chip].push_back(fit_status);
            // Note this will be active after the pedestal is subtracted
            if(this->is_channel_masked(chip,ichan))
            {
                output[chip].first.push_back(0.0);
                output[chip].second.push_back(0.0);
                // doesn't matter the fit
                fitstatus_m[chip][ichan] = 0;
            }
            else
            {
                output[chip].first.push_back( (gausfunc[chip][ichan])->GetParameter(1) );
                output[chip].second.push_back( (gausfunc[chip][ichan])->GetParameter(2) );
            }
        }
    }
    // Check if anything went wrong, 
    for(const auto & chip_v: fitstatus_m)
    {
        for(int ichan = 0; ichan < static_cast<int>(chip_v.second.size()); ++ichan)
        {
            if(chip_v.second[ichan] != 0)
            {
                std::cout << "\033[1;33mWARNING\033[1;m Problems with the "
                    << "gaussian fit @ CHIP: " << chip_v.first 
                    << " , STRIP: " << ichan << " (Fit status:" 
                    << chip_v.second[ichan] << ")" << std::endl;
            }
        }
    }
    // Deallocating memory
    delete cc;
    cc = nullptr;
    // ROOT stuff
    auxmem::deallocate_memory<TH1F>(histos);
    auxmem::deallocate_memory<TF1>(gausfunc);
    // vector and maps
    auxmem::deallocate_memory(thedata);
    
    // RE-activate all branches (and reset their object addresses)
    pedestal.reset_events_tree();

    // Mask noisy channels only if this is the first call of this
    // function, before do the subtraction
    if(!this->_pedestal_subtracted)
    {
        this->mask_channels(output);
    }

    return output;
}


void AlibavaPostProcessor::get_pedestal_noise_free(const IOManager & pedestal, const PedestalNoiseBeetleMap & mean_ped_map)
{
    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { "data_beetle1", "data_beetle2"};
    pedestal.set_events_tree_access(data_names);
    
    // Helper map 
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    pedestal.set_events_tree_branch_address("data_beetle1",&(thedata[0]));
    pedestal.set_events_tree_branch_address("data_beetle2",&(thedata[1]));

    // Helper maps 
    std::map<int,std::vector<float>*> postproc_thedata = { {0,nullptr}, {1,nullptr} };

    // The new tree 
    TTree * t = new TTree(this->_postproc_treename.c_str(),"Post-processed Events");
    // Create the new branches
    // PEdestal with the common mode subtracted
    t->Branch("postproc_data_beetle1",&(postproc_thedata[0]));
    t->Branch("postproc_data_beetle2",&(postproc_thedata[1]));
    
    // loop over the tree to fill the histograms
    const int nentries = pedestal.get_events_number_entries();
    for(int k = 0; k < nentries; ++k)
    {
        for(auto & chip_v: postproc_thedata)
        {
            if(chip_v.second != nullptr)
            {
                chip_v.second->clear();
            }
        }
        /*if( postproc_thedata[0] != nullptr )
        {
            postproc_thedata[0]->clear();
        }
        if( postproc_thedata[1] != nullptr )
        {
            postproc_thedata[1]->clear();
        }*/

        pedestal.get_events_entry(k);
        
        std::vector<float> nullsignal;
        nullsignal.reserve(ALIBAVA::NOOFCHANNELS);
        for(int chip = 0; chip < ALIBAVA::NOOFCHIPS; ++chip)
        {
            nullsignal.clear();
            for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
            {
                // Obtain the null-signal ADCs, i.e subtract to each pedestal 
                // the calculated <pedestal>_evts for this channel)
                if(this->is_channel_masked(chip,ichan))
                {
                    nullsignal.push_back(0.0);
                }
                else
                {
                    nullsignal.push_back( (*(thedata[chip]))[ichan]-mean_ped_map.at(chip).first[ichan] );
                }
            }
            // Now calculate the common noise over the null-signal ADCs
            std::pair<float,float> cmmd_cerr= this->calculate_common_noise(nullsignal,this->_channel_mask[chip]);
            // And Subtract the common noise to the pedestals: 
            // pedestals without the common noise
            for(int ichan=0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
            {
                if(this->is_channel_masked(chip,ichan))
                {
                    postproc_thedata[chip]->push_back(0.0);
                }
                else
                {
                    postproc_thedata[chip]->push_back((*(thedata[chip]))[ichan]-cmmd_cerr.first);
                }
            }
        }
        t->Fill();
    }
    // To store it with it (don't deallocate it yet, it will be done
    // afterwards)
    pedestal.get_events_tree()->AddFriend(this->_postproc_treename.c_str());
    
    // RE-activate all branches (and reset their object addresses)
    pedestal.reset_events_tree();

    // Let the instance know the pedestal was subtracted, therefore
    // the postprocEvents tree exits
    this->_pedestal_subtracted = true;

    // deallocate memory
    auxmem::deallocate_memory(thedata);
    auxmem::deallocate_memory(postproc_thedata);
}

// Overloaded with no masked channels
std::pair<float,float> AlibavaPostProcessor::calculate_common_noise(const std::vector<float> & nullsignal)
{
    std::vector<int> ch_m(ALIBAVA::NOOFCHANNELS,1);
    return AlibavaPostProcessor::calculate_common_noise(nullsignal,ch_m);
}

std::pair<float,float> AlibavaPostProcessor::calculate_common_noise(const std::vector<float> & nullsignal,const std::vector<int> & channel_mask)
{
    // The criteria to discard signal ADC entries
    const float _NoiseDeviation = 2.5;

    // Map to keep only non-signal ADCs 
    // ichannel: ADC
    std::map<int,float> non_signal_map;
    
    // initialize the map
    for(int i=0; i < static_cast<int>(nullsignal.size()); ++i)
    {
        if(channel_mask[i] == 0)
        {
            continue;
        }
        non_signal_map.emplace(i,nullsignal[i]);
    }
    
    // 1. Obtain the mean and the standard deviation for the signal ADCs 
    //    (pedestal substracted) (i.e. the common noise)
    float mean_signal  = get_mean(nullsignal);
    float stddev_signal= get_std_dev(nullsignal,mean_signal);
    
    unsigned int last_vec_size = non_signal_map.size();
    do 
    {
        // after the check form the while statament, update the last 
        // vector size 
        last_vec_size = non_signal_map.size();
        // 2. Use the mean and standard deviation to exclude real
        //    signal presence. Criteria: anything above 2.5 sigmas from the mean
        //    is considered signal
        for(auto it = non_signal_map.cbegin(); it != non_signal_map.cend(); /* no increment*/)
        {
            if(channel_mask[it->first] == 0)
            {
                ++it;
                continue;
            }
            // Remove the signal channels 
            if( std::abs(it->second-mean_signal)/stddev_signal > _NoiseDeviation )
            {
                non_signal_map.erase(it++);
            }
            else
            {
                // or keep it
                ++it;
            }
        }
        // 3. Recalculate the mean with the excluded signal channels
        mean_signal  = get_mean(non_signal_map);
        stddev_signal= get_std_dev(non_signal_map,mean_signal);
        // and the new vector size
    } while( non_signal_map.size() != last_vec_size );
    
    // fill the vector with all the same element
    return std::make_pair(mean_signal,stddev_signal);    
}

std::vector<float> AlibavaPostProcessor::convert_map_in_vector(const std::map<int,float> & m)
{
    std::vector<float> v;
    v.reserve(m.size());
    for(const auto & _f: m)
    {
        v.push_back(_f.second);
    }
    return v;
}

float AlibavaPostProcessor::get_mean(const std::map<int,float> & m)
{
    return get_mean(convert_map_in_vector(m));
}

float AlibavaPostProcessor::get_mean(const std::vector<float> & v)
{
    if(v.size() == 0)
    {   
        return 0.0;
    }
    return std::accumulate(v.begin(),v.end(),0.0)/float(v.size());
}

float AlibavaPostProcessor::get_std_dev(const std::map<int,float> & m,const float & mean)
{
    return get_std_dev(convert_map_in_vector(m),mean);
}

float AlibavaPostProcessor::get_std_dev(const std::vector<float> & v,const float & mean)
{
    if(v.size()==0)
    {
        return 0.0;
    }
    // standard deviation = sqrt( E[(x-E[x])^2] )
    std::vector<float> diff(v.size());
    // Get x-E[x] in `diff`
    std::transform(v.begin(),v.end(),diff.begin(),
            [mean](const float & element) { return element-mean; } );
    // obtain the (x-E[x])^2 by using diff*diff (inner_product)
    const float sq_sum = std::inner_product(diff.begin(),diff.end(),diff.begin(),0.0);
    return std::sqrt(sq_sum/static_cast<float>(diff.size()));
}

