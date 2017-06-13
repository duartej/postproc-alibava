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
#include "TF1.h"
#include "TCanvas.h"

// System headers
#include <fstream>
#include <map>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

AlibavaPostProcessor::AlibavaPostProcessor():
    _pedestal_subtracted(false),
    _postproc_treename("postproc_Events")
{
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
            output[chip].first.push_back( (gausfunc[chip][ichan])->GetParameter(1) );
            output[chip].second.push_back( (gausfunc[chip][ichan])->GetParameter(2) );
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
    // Deleting memory
    delete cc;
    cc = nullptr;
    for(int chip=0; chip < ALIBAVA::NOOFCHIPS; ++chip)
    {
        for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
        {
            if( histos[chip][ichan] != nullptr )
            {
                delete histos[chip][ichan];
                histos[chip][ichan] = nullptr;
            }
            if( gausfunc[chip][ichan] != nullptr )
            {
                delete gausfunc[chip][ichan];
                gausfunc[chip][ichan] = nullptr;
            }
        }
    }
    // RE-activate all branches (and reset their object addresses)
    pedestal.reset_events_tree();

    return output;
}


void AlibavaPostProcessor::get_pedestal_noise_free(const IOManager & pedestal, 
      const PedestalNoiseBeetleMap & mean_ped_map)
{
    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { "data_beetle1", "data_beetle2"};
    pedestal.set_events_tree_access(data_names);
    
    // Helper map 
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    pedestal.set_events_tree_branch_address("data_beetle1",&(thedata[0]));
    pedestal.set_events_tree_branch_address("data_beetle2",&(thedata[1]));

    // The new tree
    TTree * t = new TTree(this->_postproc_treename.c_str(),"Post-processed Events");
    // Create the new branches
    // Helper map 
    std::map<int,std::vector<float>*> postproc_thedata = { {0,nullptr}, {1,nullptr} };
    t->Branch("postproc_data_beetle1",&(postproc_thedata[0]));
    t->Branch("postproc_data_beetle2",&(postproc_thedata[1]));
    
    // loop over the tree to fill the histograms
    const int nentries = pedestal.get_events_number_entries();
    for(int k = 0; k < nentries; ++k)
    {
        if( postproc_thedata[0] != nullptr )
        {
            postproc_thedata[0]->clear();
        }
        if( postproc_thedata[1] != nullptr )
        {
            postproc_thedata[1]->clear();
        }

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
                nullsignal.push_back( (*(thedata[chip]))[ichan]-mean_ped_map.at(chip).first[ichan] );
            }
            // Now calculate the common noise over the null-signal ADCs
            std::pair<float,float> cmmd_cerr= this->calculate_common_noise(nullsignal);
            // And Subtract the common noise to the pedestals: 
            // pedestals without the common noise
            for(int ichan=0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
            {
                postproc_thedata[chip]->push_back((*(thedata[chip]))[ichan]-cmmd_cerr.first);
            }
        }
        t->Fill();
    }
    // To store it with it
    pedestal.get_events_tree()->AddFriend(this->_postproc_treename.c_str());
    
    // RE-activate all branches (and reset their object addresses)
    pedestal.reset_events_tree();

    // Let the instance know the pedestal was subtracted, therefore
    // the postprocEvents tree exits
    this->_pedestal_subtracted = true;
}

std::pair<float,float> AlibavaPostProcessor::calculate_common_noise(const std::vector<float> & nullsignal)
{
    /*if(!this->_pedestal_subtracted)
    {
        throw std::logic_error("\033[1;31mAlibavaPostProcessor ERROR\033[1;m "\
                "Invalid call to calculate_common_noise, the function must "\
                "be called after AlibavabaPostProcessor::subtract_pedestal");
    }
    
    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { "data_beetle1", "data_beetle2"};
    pedestal.set_events_tree_access(data_names);
    
    // Helper map (Note that the branches should exist in the postproc event
    // created with the pedestal subtraction function)
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    pedestal.set_events_tree_branch_address("postproc_data_beetle1",&(thedata[0]));
    pedestal.set_events_tree_branch_address("postproc_data_beetle2",&(thedata[1]));*/

    // The criteria to discard signal ADC entries
    const float _NoiseDeviation = 2.5;

    // 1. Obtain the mean and the standard deviation for the signal ADCs 
    // //    (pedestal substracted) (i.e. the common noise)
    float mean_signal  = this->get_mean(nullsignal);
    float stddev_signal= this->get_std_dev(nullsignal,mean_signal);
    // 2. Use the mean and standard deviation to exclude real
    //    signal presence. Criteria: anything above 2.5 sigmas from the mean
    //    is considered signal
    // Map to keep only non-signal ADCs 
    // ichannel: ADC
    std::map<int,float> non_signal_map;
    // XXX: Mask them ?
    // std::vector<int> masked_channels;
    // initialize the map
    for(int i=0; i < static_cast<int>(nullsignal.size()); ++i)
    {
        non_signal_map.emplace(i,nullsignal[i]);
    }
    unsigned int last_vec_size = nullsignal.size();
    while( non_signal_map.size() != last_vec_size )
    {
        for(auto it = non_signal_map.cbegin(); it != non_signal_map.cend(); /* no increment*/)
        {
            // Remove the signal channels 
            if( std::abs(it->second-mean_signal)/stddev_signal > _NoiseDeviation )
            {
                non_signal_map.erase(it);
                //masked_channels.push_back(it->first);
            }
            else
            {
                // or keep it
                ++it;
            }
        }
        // 3. Recalculate the mean with the excluded signal channels
        mean_signal  = this->get_mean(non_signal_map);
        stddev_signal= this->get_std_dev(non_signal_map,mean_signal);
        // and the new vector size
        last_vec_size = non_signal_map.size();
    }
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
    return this->get_mean(this->convert_map_in_vector(m));
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
    return this->get_std_dev(this->convert_map_in_vector(m),mean);
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
