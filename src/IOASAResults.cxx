// *************************************************
// 
// Storage of results from the AlibavaSensorAnalysis
// class
//
// *************************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "IOASAResults.h"

// The alibava related classes
#include "AlibavaPostProcessor.h"
#include "ALIBAVA.h"

// ROOT 
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TTree.h"

// System headers
#include <algorithm>
#include <cmath>

// FIXME Promote to a generic place (it is used also by IOManager)
namespace auxtree
{
    void set_tree_access(TTree * tree,const std::vector<std::string> & branch_list) 
    {
        // XXX: No mechanism to check the validity of the branches (id they exist)
        // Speed up access, just using the branches we want
        tree->SetBranchStatus("*",0);
        for(const auto & branch: branch_list)
        {
            // Be sure the string contains anything
            if(branch.empty())
            {
                continue;
            }
            tree->SetBranchStatus(branch.c_str(),1);
        }
    }
    
    void reset_tree(TTree * tree) 
    {
        tree->ResetBranchAddresses();
        tree->SetBranchStatus("*",1);
    }
}

IOASAResults::IOASAResults(const std::string & filename):
    _file(nullptr),
    _tree(nullptr)
{
    // XXX: Initialize the file 
    _file = new TFile(filename.c_str(),"RECREATE");
    _tree = new TTree("AlibavaAnalysis","Alibava cluster analysis");

}

IOASAResults::~IOASAResults()
{
    for(auto & h: _histos)
    {
        if(h.second != nullptr)
        {
            delete h.second;
            h.second = nullptr;
        }
    }
}

std::map<std::string,std::string> IOASAResults::get_branch_names(bool was_pedfile_processed)
{
    // Note the definition of the placeholder "@" in order to be substituted for the chip number
    //
    // Initialize with the independent (of the pedestal file processsing) branches
    // Plotname: { {x-branchname,y-branchname}, cutstring }
    std::map<std::string,std::string> branches_map = 
    {
        {"temperature","temperature"},
        {"tdc","eventTime"}
    };

    const std::string chipnum = std::to_string(_chip_number);

    if(was_pedfile_processed)
    {
        branches_map["pedestal"]   = "pedestal_cmmd_beetle"+chipnum;
        branches_map["noise"]      = "noise_cmmd_beetle"+chipnum;
        branches_map["commonmode"] = "postproc_cmmd_beetle"+chipnum;
        branches_map["noiseevent"] = "postproc_cmmd_noise_beetle"+chipnum;
        branches_map["signal"]     = "postproc_data_beetle"+chipnum;
        branches_map["timeprofile"]= "eventTime";
    }
    else
    {
        branches_map["pedestal"]   = "header_pedestal";
        branches_map["noise"]      = "header_noise";
        branches_map["signal"]     = "data_beetle"+chipnum;
        branches_map["timeprofile"]= "eventTime";
    }

    return branches_map;
}

void IOASAResults::book_plots()
{

    // The monitoring plot initialization
    // ----------------------------------
    // Note that most of the histograms are defined here
    // with a dummy number of bins and/or axis ranges, 
    // they are going to be fixed dynamically afterwards
    
    TH1::AddDirectory(0);

    // 1. Calibration (needs a calibration plot)
    const std::string suffix_name("chip-"+std::to_string(_chip_number));
    _histos["calibration"] = new TH2F(std::string("histo_Calibration_"+suffix_name).c_str(),
            std::string("Calibration curves "+suffix_name+";Channel number;Injected pulse [e] x10^{3};Signal average [ADC]").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1,
            1,0,1);
    // 2. Pedestal/Noise (needs pedestal)
    _histos["pedestal"] = new TGraph();
    static_cast<TGraph*>(_histos["pedestal"])->SetName(std::string("histo_Pedestals_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["pedestal"])->SetTitle(std::string("Pedestals "+suffix_name+";Channel number; Pedestal [ADC]").c_str());
    _histos["noise"] = new TGraph();
    static_cast<TGraph*>(_histos["noise"])->SetName(std::string("histo_Noise_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["noise"])->SetTitle(std::string("Noise "+suffix_name+";Channel number; Noise [ADC]").c_str());
    
    // 3. Temperature
    _histos["temperature"] = new TGraph();
    static_cast<TGraph*>(_histos["temperature"])->SetName(std::string("histo_Temperature_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["temperature"])->SetTitle(std::string("Temperature per event "+suffix_name+";Event number; Temperature [^{o}C]").c_str());

    // 4. TDC
    _histos["tdc"] = new TGraph();
    static_cast<TGraph*>(_histos["tdc"])->SetName(std::string("histo_TDC_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["tdc"])->SetTitle(std::string("TDC per event "+suffix_name+";Event number; Time [ns]").c_str());

    // 5. Signal (needs pedestal)
    _histos["signal"] = new TH1F(std::string("histo_Signal_"+suffix_name).c_str(),
            std::string("Signal "+suffix_name+";Signal [ADC]; Entries").c_str(),
            1200,-600,600);
    // 6. Hits   (needs pedestal)
    _histos["hits"] = new TH1F(std::string("histo_Hits_"+suffix_name).c_str(),
            std::string("Hits "+suffix_name+";Channel number; Entries").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1);
    // 7. Time profile (needs pedestal)
    _histos["timeprofile"] = new TProfile(std::string("histo_TimeProfile_"+suffix_name).c_str(),
            std::string("Time profile "+suffix_name+";Time [ns]; Signal average [ADC]").c_str(),
            50,0,100);
    // 8. Noise per event (needs pedestal)
    _histos["noiseevent"] = new TGraph();
    static_cast<TGraph*>(_histos["noiseevent"])->SetName(std::string("histo_NoiseEvent_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["noiseevent"])->SetTitle(std::string("Noise per event "+suffix_name+";Event number; Noise [ADC]").c_str());
    
    _histos["commonmode"] = new TGraph();
    static_cast<TGraph*>(_histos["commonmode"])->SetName(std::string("histo_NoiseEvent_"+suffix_name).c_str());
    static_cast<TGraph*>(_histos["commonmode"])->SetTitle(std::string("Common mode"+suffix_name+";Event number; Noise [ADC]").c_str());

    TH1::AddDirectory(1);
}

void IOASAResults::book_plot(const std::string & name, const TObject * theplot)
{
    // Consistency check
    if(_histos.find(name) != _histos.end())
    {
        std::cerr << "[IOASAResults::book_plot] The plot object '"
            << name << "' already exists or it is already booked a different object "
            << " with the same name" << std::endl;
        // Throw exception?
        return;
    }

    // Get the object and clone it without appending it in the gDirectory of TROOT
    const std::string suffix_name("chip-"+std::to_string(this->_chip_number));
    _histos[name] = theplot->Clone(std::string("externalbook_histo_"+name+"_"+suffix_name).c_str());
    // Only de-attach from current ROOT-directory if it was
    if(std::string(_histos[name]->ClassName()).find("TGraph") != 0)
    {
        dynamic_cast<TH1*>(_histos[name])->SetDirectory(0);
    }
}

void IOASAResults::fill_diagnostic_plots(TTree * event_tree, TTree * header_tree)
{
    // The workhorse method of this class: the histograms are filled by running over
    // the relevant trees. 
    
    // Check if the pedestal file was processed
    bool was_pedfile_proc = false;
    if(std::string(event_tree->GetName()).find("postproc_") == 0)
    {
        was_pedfile_proc = true;
    }

    // -- Get the list of branches needed for each plot  
    //auto br_map = this->get_needed_branches(was_pedfile_proc); --> Draw option
    // -- Get the equivalent branch name for the plot
    auto br_map = this->get_branch_names(was_pedfile_proc);
    // -- Attach the trees in order to access them consistently
    event_tree->AddFriend(header_tree);

    
    // 1A. Process it
    // -- Activate only those branches which are going to be used
    std::vector<std::string> branch_names;
    std::for_each(br_map.begin(),br_map.end(), 
            [&branch_names] (const std::pair<std::string,std::string> & pl_br) { branch_names.push_back(pl_br.second); });
    auxtree::set_tree_access(event_tree,branch_names);

    // Attach the objects to access them
    std::map<std::string,float> float_obj = { {"eventTime",-9999.9}, {"temperature",-9999.9} };
    std::map<std::string,std::vector<float>* > vec_float_obj;
    // -- Note that those that are not float, are vector<float*>
    for(auto & brname: branch_names)
    {
        // It is already attached in the float_obj 
        if(float_obj.find(brname) != float_obj.end())
        {
            continue;
        }
        // Also these branches are floats, update the map as well
        if(brname.find("postproc_cmmd_beetle") == 0 
                || brname.find("postproc_cmmd_noise_beetle") == 0)
        {
            float_obj[brname] = -9999.9;
            continue;
        }
        // Otherwise, it is a vector of floats
        vec_float_obj[brname] = nullptr;
    }
    // Now attach the floats to the tree (avoid to do it in the previous loop
    // in order to keep controlled the memory re-allocation)
    for(auto & name_f: float_obj)
    {
        event_tree->SetBranchAddress(name_f.first.c_str(),&(name_f.second));
    }
    // And attach the vectors (avoid to do it in the previous loop
    // in order to keep controlled the memory re-allocation)
    for(auto & name_vect: vec_float_obj)
    {
        event_tree->SetBranchAddress(name_vect.first.c_str(),&(name_vect.second));
    }
    
    // The signal data (to be corrected by the pedestal 
    // and  the common mode if needed
    std::vector<float>* sg_pedestal_free = nullptr;
    std::pair<float,float> cmmd;
    // The noise and pedestal vector are needed to differenciate the cases between pedestal file
    // processed or not. In the case of the not processed pedestal file, the vector obtained 
    // from the header contains all the values for all the chips
    std::vector<float> pedestal(ALIBAVA::NOOFCHANNELS);
    std::vector<float> noise(ALIBAVA::NOOFCHANNELS);
    
    // Go back to the previous position, move up 2 line, and set 
    // the 80 column
    std::cout << std::endl;
    std::cout << "\033[2A\033[45C["<<this->_chip_number<<"][000%]" << std::flush;
    const int nentries= event_tree->GetEntries();
    float point = float(nentries)/(100.0);
    // the event loop
    for(int k=0; k < nentries; ++k)
    {
        std::cout << "\r\033[45C["<<this->_chip_number<<"][" << std::setfill('0') << std::setw(3) 
            << std::round(float(k)/point) << "%]" << std::flush;

        event_tree->GetEntry(k);
        // Calculate pedestal and noise (just using first entry as they belong 
        // to the runHeader tree)
        if(k == 0)
        {
            // Extract the pedestal and noise
            if(!was_pedfile_proc)
            {
                // Note here the pedestal,noise are extracted from the header, and
                // they contain the two chips together in the same vector
                // WARNING: the chip number convention used in this class 
                //      (IOASAResults chips start at 1) is different from
                //      the used along the package (chips start at 0)
                unsigned int el = static_cast<unsigned int>((this->_chip_number-1)*ALIBAVA::NOOFCHANNELS);
                for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
                {
                    pedestal[ichan] = (*vec_float_obj[br_map["pedestal"]])[el];
                    noise[ichan] = (*vec_float_obj[br_map["noise"]])[el];
                    ++el;
                }

            }
            else
            {
                for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ++ichan)
                {
                    pedestal[ichan] = (*vec_float_obj[br_map["pedestal"]])[ichan];
                    noise[ichan] = (*vec_float_obj[br_map["noise"]])[ichan];
                }
            }
            // And fill the corresponding graphs
            for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
            {
                static_cast<TGraph*>(_histos["pedestal"])->SetPoint(ich,ich,pedestal[ich]);
                static_cast<TGraph*>(_histos["noise"])->SetPoint(ich,ich,noise[ich]);
            }
        }
        // The data signal (if there was no pedestal file processed,
        // we must calculate it now
        sg_pedestal_free = vec_float_obj[br_map["signal"]];
        if(!was_pedfile_proc)
        {
            for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
            {
                (*sg_pedestal_free)[ich] -= pedestal[ich];
            }
            cmmd = AlibavaPostProcessor::calculate_common_noise(*sg_pedestal_free);
        }
        else
        {
            // Get them from the branches
            cmmd = std::pair<float,float>(float_obj[br_map["commonmode"]],float_obj[br_map["noiseevent"]]);
        }

        // Now filling the histograms
        // Per event plots
        // -- temperature and tdc
        static_cast<TGraph*>(_histos["temperature"])->SetPoint(k,k,float_obj["temperature"]);
        static_cast<TGraph*>(_histos["tdc"])->SetPoint(k,k,float_obj["eventTime"]);
        // -- common mode and noise
        static_cast<TGraph*>(_histos["commonmode"])->SetPoint(k,k,cmmd.first);
        static_cast<TGraph*>(_histos["noiseevent"])->SetPoint(k,k,cmmd.second);
        // -- per channel plots
        for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
        {
            // subtract the common mode if needed
            if(!was_pedfile_proc)
            {
                (*sg_pedestal_free)[ich] -= cmmd.first;
            }
            
            // Evaluate if this channel is defining the seed of a cluster
            // (a rough estimation of the hits --see description in the method--)
            if(this->is_seed_cluster((*sg_pedestal_free)[ich],noise[ich]))
            {
                // Store the signal (pedestal and common-free)
                // WARNING: The stored signal is refered only to the seed channel,
                //          therefore a cut in ~ |S| > 5 <Noise> will be present 
                //          in the signal distribution
                static_cast<TH1F*>(_histos["signal"])->Fill((*sg_pedestal_free)[ich]);
                static_cast<TH1F*>(_histos["hits"])->Fill(ich);
                static_cast<TProfile*>(_histos["timeprofile"])->Fill(float_obj["eventTime"],(*sg_pedestal_free)[ich]);
            }
        }
    }
    std::cout << std::endl;

    auxtree::reset_tree(header_tree);
    auxtree::reset_tree(event_tree);
}

const std::vector<TObject*> IOASAResults::get_calibration_plots() const
{
    // Get the calibration plots booked externally
    // XXX FIXME - the name should be defined as a global (in the ALIBAVA.h probably)
    std::vector<TObject*> theobjects;
    std::string calname("calibration_profile_");
    for(int ichannel = 0; ichannel < ALIBAVA::NOOFCHANNELS; ++ichannel)
    {
        // -- first check whether they are here or not
        const std::string name(calname+std::to_string(ichannel));
        if(_histos.find(name) == _histos.end())
        {
            std::cerr << "[IOASAResults::get_calibration_plots] The plot object '"
                << name << "' was not booked! Inconsistent use of this method...'" << std::endl;
            // XXX throw exception ?
            continue;
        }
        theobjects.push_back(_histos.at(name));
    }
    return theobjects;
}

void IOASAResults::set_calibration_plot(const std::vector<TObject*> & curves)
{
    // First check at least there is something
    if(curves.size() == 0)
    {
        std::cerr << "[IOASAResults::set_calibration_plots] No calibration plots "
            << "present in the argument. Incoherent use of this class... '" << std::endl;
        return;
    }
    // Get the info from the objects in order to re-bin the histogram
    TProfile * auxobj = static_cast<TProfile*>(curves[0]);
    const int bins_y = auxobj->GetNbinsX();
    const float ymin = auxobj->GetXaxis()->GetBinLowEdge(1);
    const float ymax = auxobj->GetXaxis()->GetBinUpEdge(bins_y);
    
    // a handler to work with the histogram
    TH2F * thehisto = static_cast<TH2F*>(_histos["calibration"]);
    
    // Rebin the histogram depending the profiles bins
    thehisto->SetBins(thehisto->GetNbinsX(),thehisto->GetXaxis()->GetBinLowEdge(1),thehisto->GetXaxis()->GetBinUpEdge(thehisto->GetNbinsX()),
            bins_y,ymin,ymax);
    // And fill the histogram (remember was it fill in channel order)
    for(unsigned int i=0; i < curves.size(); ++i)
    {
        // The vector elements orden corresponds to the channel number
        const int ich = static_cast<int>(i);
        const TProfile* prf = static_cast<TProfile*>(curves[i]);
        // for each channel, let's obtain the injected charge (x-profile) with
        // its measured ADC counts (y-profile)
        for(int ibin=1; ibin < prf->GetNbinsX()+1; ++ibin)
        {
            const int bin = thehisto->FindBin(ich,prf->GetBinCenter(ibin));
            thehisto->SetBinContent(bin,prf->GetBinContent(ibin));
        }
    }
}

void IOASAResults::deliver_plots()
{
    // Store all external booked plots if any
    // XXX: a directory can be used to create a more ordered structure
    // in that case it is needed either a TFile argument or a gDirectory...
    for(auto & h: _histos)
    {
        if(h.second != nullptr && std::string(h.second->GetName()).find("externalbook") == 0)
        {
            h.second->Write();
        }
    }

    // Plot all the monitoring plots in a canvas and 
    // store that canvas in the file
    // ----------------------------------------------
    // XXX FIXME: Define the Style
    // Define the canvas
    const int wsize = 650;
    const int hsize = 850;
    TCanvas * canvas = new TCanvas(std::string("monitor_plots_beetle"+std::to_string(_chip_number)).c_str(),
            std::string("Monitor plots Beetle "+std::to_string(_chip_number)).c_str(),
            wsize,hsize);
    // and its layout: colummns, row, separation betwwen columns and rows
    canvas->Divide(2,4,0.01,0.01);
    // And plot (first define the order of plotting (except those pads
    // with more than 1 histo per pad)
    // Note that the map is defining the pad where the histogram should go
    // and the Draw option to be used
    const std::map<int,std::pair<std::string,std::string>> plotorder = 
    { 
        {1,{"calibration","SURF3Z"}},
        {3,{"signal",""}},
        {4,{"hits",""}},
        {5,{"timeprofile",""}},
        {6,{"temperature",""}},
        {7,{"tdc",""}}
    };
    // Define the pads and fill them
    for(const auto & padplot: plotorder)
    {
        canvas->cd(padplot.first);
        // Check if the plot is present
        if(_histos[padplot.second.first] == nullptr)
        {
            // This behaviour is also allowed (see calibration or pedestal
            // processing cases)
            delete canvas;
            canvas = nullptr;
            return;
        }
        _histos[padplot.second.first]->Draw(padplot.second.second.c_str());
    }
    // And now the composite pads
    // --------------------------
    // PAD 2: Pedestals - Noise
    auto pad_comp = canvas->cd(2);
    pad_comp->Divide(2,1,0.001,0.001);
    pad_comp->cd(1);
    _histos["pedestal"]->Draw();
    pad_comp->cd(2);
    _histos["noise"]->Draw();
    // PAD 8: Noise per Event 
    auto pad_comp8 = canvas->cd(8);
    pad_comp8->Divide(2,1,0.001,0.001);
    pad_comp8->cd(1);
    _histos["noiseevent"]->Draw();
    pad_comp8->cd(2);
    _histos["commonmode"]->Draw();
    // Save it to the TFile present in the gDirectory 
    // (managed by the client of this class: IOManager)
    canvas->Write();
    // and deallocate memory 
    delete canvas;
    canvas = nullptr;
    // Note that the pads are already deleted by the canvas
}

bool IOASAResults::is_seed_cluster(const float & signal, const float & noise, const float & SoNtimes) const
{
    // This method is just evaluating if the amount of signal is beyond the SoNtimes x noise
    // and in that case, it is assuming that it could be seed of a cluster.
    // WARNING: this is not a fully finding cluster algorithm, as ignores to include the 
    // neighbour channels, it is just giving an indication of the hit profile
    return std::fabs(signal) > SoNtimes*noise;
}
