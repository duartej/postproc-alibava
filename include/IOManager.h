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

#include <string>

// forward declarations
class TFile;
class TTree;

class AlibavaRunHeader;
class AlibavaEvent;

class IOManager
{
    private:
        // datamembers
        TFile * _file;
        TTree * _tree_header;
        TTree * _tree_events;
        int     _eventsProcessed;

        // The auxiliary functions
        AlibavaRunHeader * _runheader;
        AlibavaEvent * _events;

    public:
        IOManager(const std::string & rootfilename);
        ~IOManager();

        void book_tree_header();
        void book_tree();
        
        void fill_header(const AlibavaRunHeader * aheader) const;
        void fill_event(const AlibavaEvent * anAlibavaEvent) const;
        
        void close();
};

#endif
