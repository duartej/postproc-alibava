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
#include "AlibavaSensorAnalysis.h"
#include "IOFortythieves.h"
#include "StripCluster.h"

// ROOT 
#include "TFile.h"
#include "TTree.h"

// System headers
#include <algorithm>
#include <cmath>
#include <memory>

IOASAResults::IOASAResults(const std::string & filename):
    _file(nullptr),
    _tree(nullptr),
    _event_number(-1),
    _polarity(0),
    _cluster_number(-1),
    _event_time(-1),
    _temperature(-1),
    _common_mode(-99999.9),
    _event_noise(-99999.9),
    _cluster_size(nullptr),
    _cluster_seed_channel(nullptr),
    _cluster_charge(nullptr),
    _cluster_seed_charge(nullptr),
    _cluster_eta_seed(nullptr),
    _cluster_eta(nullptr),
    _cluster_channels(nullptr)
{
    // Initialize the file 
    _file = new TFile(filename.c_str(),"RECREATE");

    // Associate all the vectors to some helper containers
    _branches_int.push_back(&_cluster_size);
    _branches_int.push_back(&_cluster_seed_channel);

    _branches_float.push_back(&_cluster_charge);
    _branches_float.push_back(&_cluster_seed_charge);
    _branches_float.push_back(&_cluster_eta_seed);
    _branches_float.push_back(&_cluster_eta);

    // Allocate memory
    for(auto & el: _branches_int)
    {
        *el = new std::vector<int>;
    }
    for(auto & el: _branches_float)
    {
        *el = new std::vector<float>;
    }
}

IOASAResults::~IOASAResults()
{
    // Storing results at the file
    _file->Write(nullptr,TObject::kOverwrite);
    _file->Close();
    // deallocating memory
    /*if(_file != nullptr)
    {
        delete _file;
        _file = nullptr;
    }
    if(_tree != nullptr)
    {
        delete _tree;
        _tree = nullptr;
    }
    
    // Freeing memory
    for(auto & el: _branches_int)
    {
        if(el != nullptr)
        {
            delete el;
            el = nullptr;
        }
    }
    for(auto & el: _branches_float)
    {
        if(el != nullptr)
        {
            delete el;
            el = nullptr;
        }
    }*/
    
    /*for(auto & h: _histos)
    {
        if(h.second != nullptr)
        {
            delete h.second;
            h.second = nullptr;
        }
    }*/
}

void IOASAResults::clear_variables()
{
    // Plain variables
    _event_number = 0;
    _polarity = 0;
    _event_time = -1;
    _temperature = -1;
    _cluster_number = -1;
    _common_mode = -99999.9;
    _event_noise = -99999.9;

    // vectors
    for(auto & el: _branches_int)
    {
        (*el)->clear();
    }
    for(auto & el: _branches_float)
    {
        (*el)->clear();
    }
}

void IOASAResults::book_tree()
{
    _tree = new TTree("alibava_clusters","Alibava cluster analysis");

    _tree->Branch("eventNumber",&(_event_number));
    _tree->Branch("eventTime",&(_event_time));
    _tree->Branch("temperature",&(_temperature));
    _tree->Branch("polarity",&(_polarity));
    _tree->Branch("common_mode",&(_common_mode));
    _tree->Branch("event_noise",&(_event_noise));
    _tree->Branch("cluster_number",&(_cluster_number));
    _tree->Branch("cluster_size",&(_cluster_size));
    _tree->Branch("cluster_charge",&(_cluster_charge));
    _tree->Branch("cluster_seed_channel",&(_cluster_seed_channel));
    _tree->Branch("cluster_seed_charge",&(_cluster_seed_charge));
    _tree->Branch("cluster_eta_seed",&(_cluster_eta_seed));
    _tree->Branch("cluster_eta",&(_cluster_eta));

    clear_variables();
}

void IOASAResults::fill_tree(const IOFortythieves * ioft_inst, AlibavaSensorAnalysis * aa_inst)
{
    // Variables from the AlibavaSensorAnalysis instance and fill 
    // the tree. Note that this method uses the IOFortythieves
    // getters at event-based, i.e., 
    
    // And fill event-based
    _event_number = ioft_inst->event_number();
    _event_time   = ioft_inst->event_time();
    _temperature  = ioft_inst->temperature();
    _polarity     = aa_inst->get_polarity();
    _common_mode  = ioft_inst->common_mode();
    _event_noise  = ioft_inst->event_noise();
    
    // Get the clusters
    const std::vector<std::unique_ptr<StripCluster> > clusters = aa_inst->find_clusters(ioft_inst);

    _cluster_number = clusters.size();

    // The cluster variables
    for(const auto & cl: clusters)
    {
        _cluster_size->push_back(cl->size());
        _cluster_charge->push_back(cl->charge());
        _cluster_seed_channel->push_back(cl->channels(0));
        _cluster_seed_charge->push_back(cl->charge(0));
        _cluster_eta_seed->push_back(cl->eta_seed()); 
        _cluster_eta->push_back(cl->eta()); 
    }
    _tree->Fill();
    
    // And prepare all variables for the next iteration
    clear_variables();
}


// C-Wrapper to be used by the python module
#ifdef __cplusplus
extern "C"
{
#endif
    // Constructor and destructor
    IOASAResults * ioresults_new(const char * filename) { return new IOASAResults(filename); }
    void ioresults_delete(IOASAResults * ioresults_inst) { ioresults_inst->~IOASAResults(); }
    // Book and fill
    void ioresults_book_tree(IOASAResults * ioresults_inst) { ioresults_inst->book_tree(); }
    void ioresults_fill_tree(IOASAResults * iores, IOFortythieves * ioft, AlibavaSensorAnalysis * aa) { iores->fill_tree(ioft,aa); }
#ifdef __cplusplus
}
#endif
