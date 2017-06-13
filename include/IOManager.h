// ***********************************************
// 
// Output manager to persistify ROOT files
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef IOMANAGER_H
#define IOMANAGER_H

#include "input_options.h"

#include <string>
#include <vector>
#include <map>

// forward declarations
class TFile;
class TTree;

class AlibavaRunHeader;
class AlibavaEvent;

using PedestalNoiseBeetleMap = std::map<int,std::pair<std::vector<float>,std::vector<float>> >;

class IOManager
{
    private:
        // datamembers
        std::string _rootfilename;
        TFile * _file;
        TTree * _tree_header;
        TTree * _tree_events;
        int     _eventsProcessed;

        // The auxiliary functions
        AlibavaRunHeader * _runheader;
        AlibavaEvent * _events;
        
        // store the friends
        void aux_store_friends(TTree * tree);

    public:
        IOManager(const std::string & rootfilename);
        ~IOManager();

        void update(const PedestalNoiseBeetleMap & pednoise_m);

        void book_tree_header();
        void book_tree();
        
        void fill_header(const AlibavaRunHeader * aheader) const;
        void fill_event(const AlibavaEvent * anAlibavaEvent) const;
        
        // Update the Tree (Runheader one) with pedestal/noise info
        int update_runheader(const IOManager & other);
        
        void close();
    
        // Event Tree related functions
        // Events tree getter
        TTree * get_events_tree() const;
        // Speed up access to the tree
        void set_events_tree_access(const std::vector<std::string> & branch_list) const;
        // Re-activate everything (and reset branch addresses)
        void reset_events_tree() const;
        // Probably a template, but just vector of floats so far
        void set_events_tree_branch_address(const std::string & branch, std::vector<float> ** v) const;
        // Number of entries and Get the entry
        int get_events_number_entries() const;
        void get_events_entry(const int & i) const;
};

#endif
