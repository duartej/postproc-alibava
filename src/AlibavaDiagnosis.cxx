// ***********************************************
// 
// Monitoring and diagnosis tool for the Alibava
// data-taking and (post-)processing 
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "AlibavaDiagnosis.h"

// The alibava run and event headers
//#include "AuxiliaryStructures.h"
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
//#include <fstream>
//#include <map>
//#include <algorithm>
//#include <numeric>
//#include <stdexcept>
//#include <utility>

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

AlibavaDiagnosis::AlibavaDiagnosis(const int & chipnumber):
    _chip_number(chipnumber)
{
}

AlibavaDiagnosis::~AlibavaDiagnosis()
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

std::map<std::string,std::pair<std::vector<std::string>,std::string> > AlibavaDiagnosis::get_needed_branches(bool was_pedfile_processed)
{
    // Note the definition of the placeholder "@" in order to be substituted for the chip number
    //
    // Initialize with the independent (of the pedestal file processsing) branches
    // Plotname: { {x-branchname,y-branchname}, cutstring }
    std::map<std::string,std::pair<std::vector<std::string>,std::string> > branches_map = 
    {
        {"temperature",{{"Entry$","temperature"},""}},
        {"tdc",{{"Entry$","eventTime"},""}}
    };

    const std::string chipnum = std::to_string(_chip_number);

    if(was_pedfile_processed)
    {
        branches_map["pedestal"]   = {{"Iteration$","pedestal_cmmd_beetle"+chipnum},""};
        branches_map["noise"]      = {{"Iteration$","noise_cmmd_beetle"+chipnum},""};
        branches_map["commonmode"] = {{"Entry$","postproc_cmmd_beetle"+chipnum},""};
        branches_map["noiseevent"] = {{"Entry$","postproc_cmmd_noise_beetle"+chipnum},""};
        branches_map["signal"]     = {{"postproc_data_beetle"+chipnum},""};
        branches_map["hits"]       = {{"Iteration$"},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
        branches_map["timeprofile"]= {{"eventTime","postproc_data_beetle"+chipnum},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
    }
    else
    {
        // WAIT!!
        branches_map["pedestal"]   = {{"Iteration$","header_pedestal"},""};
        branches_map["noise"]      = {{"Iteration$","header_noise"},""};
        branches_map["commonmode"] = {{""},""};
        branches_map["noiseevent"] = {{""},""};
        branches_map["signal"]     = {{"data_beetle"+chipnum},""};
        branches_map["hits"]       = {{"Iteration$"},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
        branches_map["timeprofile"]= {{"eventTime","postproc_data_beetle"+chipnum},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
    }

    return branches_map;
}

/*std::map<std::string,std::pair<std::vector<std::string>,std::string> > AlibavaDiagnosis::get_needed_branches(bool was_pedfile_processed)
{
    // Note the definition of the placeholder "@" in order to be substituted for the chip number
    //
    // Initialize with the independent (of the pedestal file processsing) branches
    // Plotname: { {x-branchname,y-branchname}, cutstring }
    std::map<std::string,std::pair<std::vector<std::string>,std::string> > branches_map = 
    {
        {"temperature",{{"Entry$","temperature"},""}},
        {"tdc",{{"Entry$","eventTime"},""}}
    };

    const std::string chipnum = std::to_string(_chip_number);

    if(was_pedfile_processed)
    {
        branches_map["pedestal"]   = {{"Iteration$","pedestal_cmmd_beetle"+chipnum},""};
        branches_map["noise"]      = {{"Iteration$","noise_cmmd_beetle"+chipnum},""};
        branches_map["commonmode"] = {{"Entry$","postproc_cmmd_beetle"+chipnum},""};
        branches_map["noiseevent"] = {{"Entry$","postproc_cmmd_noise_beetle"+chipnum},""};
        branches_map["signal"]     = {{"postproc_cmmd_beetle"+chipnum},""};
        branches_map["hits"]       = {{"Iteration$"},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
        branches_map["timeprofile"]= {{"eventTime","postproc_data_beetle"+chipnum},"abs(postproc_cmmd_beetle"+chipnum+") > 5.0*noise_cmmd_beetle"+chipnum};
    }

    return branches_map;
}*/


void AlibavaDiagnosis::book_plots()
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
    // Should it be assigned to no-where or belong to any file?
    // --> for all the histos, then AddDirectory(0)a
    /*for(auto & h: _histos)
    {
        dynamic_cast<TH1*>(h.second)->SetDirectory(0);
    }*/
}

void AlibavaDiagnosis::book_plot(const std::string & name, const TObject * theplot)
{
    // Consistency check
    if(_histos.find(name) != _histos.end())
    {
        std::cerr << "[AlibavaDiagnosis::book_plot] The plot object '"
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

void AlibavaDiagnosis::fill_diagnostic_plots(TTree * event_tree, TTree * header_tree)
{
    // Check if the pedestal file was processed
    bool was_pedfile_proc = false;
    if(std::string(event_tree->GetName()).find("postproc_") == 0)
    {
        was_pedfile_proc = true;
    }

    // A. Process it (needed anyway if there is no pedestal)
    // B. Using Draw methods
    
    // -- Get the list of branches needed for each plot 
    auto br_map = this->get_needed_branches(was_pedfile_proc);
    // -- Attach the trees in order to access them consistently
    event_tree->AddFriend(header_tree);

    
    // 1A. Process it
    // -- Activate only those branches which are going to be used
    std::vector<std::string> branch_names;
    for(auto & plotname_auxpair: br_map)
    {
        auto varvector= plotname_auxpair.second.first;
        for(auto & brname: varvector)
        {
            if(brname != "Iteration$" && brname != "Entry$")
            {
                branch_names.push_back(brname);
            }
        }
    }
    auxtree::set_tree_access(event_tree,branch_names);

    // Attach the objects to access them
    std::map<std::string,float> float_obj = { {"eventTime",-9999.9}, {"temperature",-9999.9} };
    std::map<std::string,std::vector<float>* > vfloat_obj;
    // -- Note that those that are not float, are vector<float*>
    for(auto & brname: branch_names)
    {
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
        vfloat_obj[brname] = nullptr;
    }
    // Now attach the floats to the tree
    for(auto & name_f: float_obj)
    {
        event_tree->SetBranchAddress(name_f.first.c_str(),&(name_f.second));
    }
    // And attach the vectors (avoid to do it in the previous loop
    // in order to keep controlled the memory re-allocation
    for(auto & name_vect: vfloat_obj)
    {
        event_tree->SetBranchAddress(name_vect.first.c_str(),&(name_vect.second));
    }
    // Some useful variables: the name of the branch containing the 
    // pedestal and noise
    const std::string pedestal_branch =(br_map["pedestal"].first)[1];
    const std::string noise_branch =(br_map["noise"].first)[1];
    // the data (pedestal,noise subtracted or not)
    const std::string signal_branch = (br_map["noise"].first)[0];
    
    // Go back to the previous position, move up 2 line, and set 
    // the 80 column
    std::cout << std::endl;
    std::cout << "\033[2A\033[45C[000%]" << std::flush;
    const int nentries= event_tree->GetEntries();
    float point = float(nentries)/(100.0);
    // the event loop
    for(int k=0; k < nentries; ++k)
    {
        std::cout << "\r\033[45C[" << std::setfill('0') << std::setw(3) 
            << std::round(float(k)/point) << "%]" << std::flush;

        event_tree->GetEntry(k);
        // Calculate pedestal and noise (just using first entry as they belong 
        // to the runHeader tree)
        if(k == 0)
        {
            for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
            {
                static_cast<TGraph*>(_histos["pedestal"])->SetPoint(ich,ich,(*vfloat_obj[pedestal_branch])[ich]);
                static_cast<TGraph*>(_histos["noise"])->SetPoint(ich,ich,(*vfloat_obj[noise_branch])[ich]);
            }
        }
        // HERE OR up
        std::vector<float>* sg_pedestal_free = vfloat_obj[signal_branch];
        std::pair<float,float> cmmd_and_noise;
        if(was_pedfile_proc)
        {
            for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
            {
                (*sg_pedestal_free)[ich] -= (*vfloat_obj[pedestal_branch])[ich];
            }
        }

        // Now filling the histograms
        // Per event plots
        // -- temperature and tdc
        static_cast<TGraph*>(_histos["temperature"])->SetPoint(k,k,float_obj["temperature"]);
        static_cast<TGraph*>(_histos["tdc"])->SetPoint(k,k,float_obj["eventTime"]);
        // -- per channel plots
        //for(int ich = 0; ich < ALIBAVA::NOOFCHANNELS; ++ich)
        //{
        //}
    }
    std::cout << std::endl;

    auxtree::reset_tree(event_tree);
    //
    // 1B. Draw methods
    /*auto drawopts = get_draw_option();
    for(auto & plotname_auxpair: br_map)
    {
        auto plotname = plotname_auxpair.first;
        auto varvector= plotname_auxpair.second.first;
        auto cutstring= plotname_auxpair.second.second;
        // Prepare the var-expression
        std::string varexpression;
        for(auto & var: varvector)
        {
            varexpression.append(var+":");
        }
        // Remove the last ':'
        varexpression.erase(varexpression.end()-1);
        // Plotting using the draw method
        event_tree->Draw(varexpression.c_str(),cutstring.c_str(),drawopts[plotname].c_str());
        // Filling
        if(std::string(_histos[plotname]->ClassName()).find("TGraph") == 0)
        {
            TGraph * thgr = dynamic_cast<TGraph*>(_histos[plotname]);
            auto vx = event_tree->GetV1();
            auto vy = event_tree->GetV2();
            for(unsigned int i=0; i < event_tree->GetSelectedRows() ; ++i)
            {
                thgr->SetPoint(i,vx[i],vy[i]);
            }
        }
        else if(std::string(_histos[plotname]->ClassName()).find("TH1F") == 0)
        {
            TH1F * th1 = dynamic_cast<TH1F*>(_histos[plotname]);
            if(th1 == nullptr)
            {
                continue;
            }
            auto v1 = event_tree->GetV1();
            for(unsigned int i=0; i < event_tree->GetSelectedRows(); ++i)
            {       
                th1->Fill(v1[i]);
            }
        }
    }*/
}

const std::vector<TObject*> AlibavaDiagnosis::get_calibration_plots() const
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
            std::cerr << "[AlibavaDiagnosis::get_calibration_plots] The plot object '"
                << name << "' was not booked! Inconsistent use of this method...'" << std::endl;
            // XXX throw exception ?
            continue;
        }
        theobjects.push_back(_histos.at(name));
    }
    return theobjects;
}

void AlibavaDiagnosis::set_calibration_plot(const std::vector<TObject*> & curves)
{
    // First check at least there is something
    if(curves.size() == 0)
    {
        std::cerr << "[AlibavaDiagnosis::set_calibration_plots] No calibration plots "
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
        

void AlibavaDiagnosis::deliver_plots()
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
    TCanvas * canvas = new TCanvas(std::string("monitor_canvas_"+std::to_string(_chip_number)).c_str(),
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
