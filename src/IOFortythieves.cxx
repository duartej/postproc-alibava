// ****************************************************
// 
// Manager to read the ROOT output for the fortythieve
// exec. https://github.com/duartej/postproc-alibava/blob/master/bin/fortythieves.cc
//
// Entry point for the post-processing and cluster 
// analysis of the ROOT created by fortythieves. 
// 
// XXX DOC IMPLEMENTATION
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "IOFortythieves.h"

//#include "ALIBAVA.h"

// ROOT headers
#include "TFile.h"
#include "TTree.h"

// System headers
//

// The names of the trees
std::list<std::string> IOFortythieves::_TREENAMES( { "runHeader", "postproc_runHeader", 
        "Events", "postproc_Events" }) ;

IOFortythieves::IOFortythieves(const std::string & filename,const int & chip) :
    _chip(chip),
    _entries(-1),
    _file(nullptr),
    _run_number(-1),
    _event_number(-1),
    _temperature(-1),
    _adc_data(nullptr)
{
    // XXX Check the chip number chip < ALIBAVA::NOOFCHIPS
    
    // Create the root file
    _file = TFile::Open(filename.c_str());
    // XXX: Check validity of the pointer
    
    // And obtain the relevant trees
    for(auto & treename: _TREENAMES)
    {
        _trees[treename] = static_cast<TTree*>(_file->Get(treename.c_str()));
        // Check entries
        if(treename == "postproc_Events")
        {
            _entries = _trees[treename]->GetEntries();
        }
        // XXX: Check validity of the pointer
        // De-activate branches (afterwards activate only the used ones)
        _trees[treename]->SetBranchStatus("*",0);
   }
    
    // keep track of the branches which are going to be used
    std::map<std::string,std::vector<std::string> > used_branches;
    used_branches["Events"].push_back("eventTime");
    used_branches["Events"].push_back("temperature");
    used_branches["Events"].push_back("eventNumber");
    used_branches["Events"].push_back("runNumber");

    // Remember the ROOT file notation for the chip is chip+1 (it has beetle 1 and 2)
    const std::string adc_data_str("postproc_data_beetle"+std::to_string(chip+1));
    used_branches["postproc_Events"].push_back(adc_data_str);

    // Read both trees at the same time
    _trees["postproc_Events"]->AddFriend(_trees["Events"]);

    // Activate the branches to be used
    for(const auto & treename_vectorbr: used_branches)
    {
        for(const auto & branch_name: treename_vectorbr.second)
        {
            _trees[treename_vectorbr.first]->SetBranchStatus(std::string(branch_name+"*").c_str(),1);
        }
    }

    // Attach the trees to the corresponding data-member
    // -- TDC time, temperature, run and event numbers
    _trees["Events"]->SetBranchAddress("eventTime",&_event_time);
    _trees["Events"]->SetBranchAddress("temperature",&_temperature);
    _trees["Events"]->SetBranchAddress("eventNumber",&_event_number);
    _trees["Events"]->SetBranchAddress("runNumber",&_run_number);

    // -- ADCs corrected by common noise and pedestals
    _trees["postproc_Events"]->SetBranchAddress(adc_data_str.c_str(),&_adc_data);
}

IOFortythieves::~IOFortythieves()
{
    // delete and de-allocate memory
    if(_file != nullptr)
    {
        _file->Close();
        delete _file;
        _file = nullptr;
    }
}

void IOFortythieves::initialize()
{
    // Obtain the pedestal, noise and calibration from the postproc 
    // runHeader
    
    // 1. Create the auxiliary vectors to set the branch
    //    address to the runHeader
    // Remember the ROOT file notation for the chip is chip+1 (it has beetle 1 and 2)
    const std::string noise_str("noise_cmmd_beetle"+std::to_string(_chip+1));
    const std::string pedes_str("pedestal_cmmd_beetle"+std::to_string(_chip+1));
    const std::string calib_str("electronADC_beetle"+std::to_string(_chip+1));
    // helper map 
    std::map<std::string,std::vector<float>* > br_names = { {noise_str,nullptr}, 
        {pedes_str,nullptr}, {calib_str,nullptr} };
    for(auto & name_vector: br_names)
    {
        _trees["postproc_runHeader"]->SetBranchStatus(std::string(name_vector.first+"*").c_str(),1);
        _trees["postproc_runHeader"]->SetBranchAddress(name_vector.first.c_str(),&(name_vector.second));
    }

    // 2. Get the first entry of the header and fill the related datamembers
    _trees["postproc_runHeader"]->GetEntry(0);
    
    // 3. Pass the filled vectors to the data members
    _noise = *(br_names[noise_str]);
    _pedestal = *(br_names[pedes_str]);
    _calibration = *(br_names[calib_str]);
    
    // 4. De-activate the tree again, not needed anymore
    _trees["postproc_runHeader"]->ResetBranchAddresses();
    _trees["postproc_runHeader"]->SetBranchStatus("*",0);
    for(auto & br_vtr: br_names)
    {
        if( br_vtr.second != nullptr)
        {
            delete br_vtr.second;
            br_vtr.second = nullptr;
        }
    }
}

void IOFortythieves::process(const int & i)
{
    // Recall that the Events treee is friend of this,
    _trees["postproc_Events"]->GetEntry(i);
}

// C-Wrapper for the python use with ctypes
#ifdef __cplusplus
extern "C"
{
#endif
    IOFortythieves * ioft_new(const char * fn, int chip) { return new IOFortythieves(fn,chip); }
    void ioft_delete(IOFortythieves * io_ft) { io_ft->~IOFortythieves(); }
    void ioft_initialize(IOFortythieves * io_ft) { io_ft->initialize(); }
    void ioft_process(IOFortythieves * io_ft, int i) { io_ft->process(i); }
    int ioft_get_entries(IOFortythieves * io_ft) { return io_ft->get_entries(); }
#ifdef __cplusplus
}
#endif
