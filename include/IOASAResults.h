// *************************************************
// 
// Storage of results from the AlibavaSensorAnalysis
// class
//
// *************************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// *************************************************
#ifndef IOASARESULTS_H
#define IOASARESULTS_H


// System headers
#include<string>
#include<map>
#include<vector>

// forward declarations
class AlibavaSensorAnalysis;
class IOFortythieves;
class TFile;
class TTree;

class IOASAResults
{
    public:
        IOASAResults(const std::string & filename);
        IOASAResults() = delete;
        ~IOASAResults();

        // Tree initialization
        void book_tree();

        // Tree filling
        void fill_tree(const IOFortythieves * ioft_inst, AlibavaSensorAnalysis * aa_inst);

    private:
        TFile *_file; 
        TTree *_tree;

        // The variables to get attached the new branches:
        // -- Simple elements
        int _event_masked;
        int _event_number;
        int _polarity;
        int _cluster_number;
        float _event_time;
        float _temperature;
        float _common_mode;
        float _event_noise;
        // -- int vectors
        std::vector<int> * _cluster_size;
        std::vector<int> * _cluster_seed_channel;
        // -- float vectors
        std::vector<float> * _cluster_charge;
        std::vector<float> * _cluster_seed_charge;
        std::vector<float> * _cluster_eta_seed;
        std::vector<float> * _cluster_eta;
        std::vector<float> * _cluster_position_cog;
        // -- maps for cluster-dependent variable
        //    cluster Id: vector
        std::map<int,std::vector<float> > * _cluster_channels;

        // Helper vector to easy and speed up common and
        // rutinary steps
        // -- for int vectos
        std::vector<std::vector<int> **> _branches_int;
        // -- for float
        std::vector<std::vector<float> **> _branches_float;

        // clear vectors members before fill them
        void clear_variables();
};

#endif
