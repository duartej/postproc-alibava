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

PedestalNoiseBeetleMap AlibavaPostProcessor::calculate_pedestal_noise(const IOManager & pedestal)
{
    // Speed up access, just using the branches we want:
    const std::vector<std::string> data_names = { "data_beetle1", "data_beetle2"};
    pedestal.set_events_tree_access(data_names);
    
    // Helper map 
    std::map<int,std::vector<float>*> thedata = { {0,nullptr}, {1,nullptr} };
    pedestal.set_events_tree_branch_address("data_beetle1",&(thedata[0]));
    pedestal.set_events_tree_branch_address("data_beetle2",&(thedata[1]));
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

    // RE-activate all branches (and reset them)
    pedestal.reset_events_tree();

    return output;
}
