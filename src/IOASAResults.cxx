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
    _event_time(-1),
    _temperature(-1),
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
    _branches_int.push_back(_cluster_size);
    _branches_int.push_back(_cluster_seed_channel);

    _branches_float.push_back(_cluster_charge);
    _branches_float.push_back(_cluster_seed_charge);
    _branches_float.push_back(_cluster_eta_seed);
    _branches_float.push_back(_cluster_eta);

    // Allocate memory
    for(auto & el: _branches_int)
    {
        el = new std::vector<int>;
    }
    for(auto & el: _branches_float)
    {
        el = new std::vector<float>;
    }
}

IOASAResults::~IOASAResults()
{
    // Storing results at the file
    _file->Close();
    delete _file;
    _file = nullptr;
    if(_tree = nullptr)
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
    }
    
    /*for(auto & h: _histos)
    {
        if(h.second != nullptr)
        {
            delete h.second;
            h.second = nullptr;
        }
    }*/
}

void clear_variables()
{
    // Plain variables
    _event_number = 0;
    _polarity = 0;
    _event_time = -1;
    _temperature = -1;

    // vectors
    for(auto & el: _branches_int)
    {
        el->clear();
    }
    for(auto & el: _branches_float)
    {
        el->clear();
    }
}

void IOASAResults::book_tree()
{
    _tree = new TTree("alibava_clusters","Alibava cluster analysis");

    _tree->Branch("eventNumber",&(_events_number));
    _tree->Branch("eventTime",&(_event_time));
    _tree->Branch("temperature",&(_temperature));
    _tree->Branch("polarity",&(_polarity));
    _tree->Branch("cluster_size",&(_cluster_size));
    _tree->Branch("cluster_charge",&(_cluster_charge));
    _tree->Branch("cluster_seed_channel",&(_cluster_seed_channel));
    _tree->Branch("cluster_seed_charge",&(_cluster_seed_charge));
    _tree->Branch("cluster_eta_seed",&(_cluster_eta_seed));
    _tree->Branch("cluster_eta",&(_cluster_eta_seed));
}

void IOASAResults::fill_tree(const AlibavaSensorAnalysis * aa_inst)
{
    // Variables from the AlibavaSensorAnalysis instance and fill 
    // the tree. Note that this method uses the IOFortythieves
    // getters at event-based, i.e., 
    
    // First clear all variables
    clear_variableS();
    // And fill event-based
    _event_number = aa_inst->get_event_number();
    _event_time = aa_inst->get_event_time();
    _temperature = aa_inst->get_temperature();
    _polarity = aa_inst->get_polarity();
    
    // Get the clusters
    auto clusters = aa_inst->find_clusters(_event_number);

    // The cluster variables
    for(const auto & cl: clusters)
    {
        _cluster_size->push_back(cl->size());
        _cluster_charge->push_back(cl->charge());
        _cluster_seed_channel->push_back(cl->channels(0));
        _cluster_seed_charge->push_back(cl->charge(0));
        // XXX
        _cluster_eta_seed->push_back(0); 
        _cluster_eta->push_back(0); 
    }
    _tree->Fill();
}



