// ***********************************************
// 
// Output manager to persistify ROOT files
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "IOManager.h"

// The alibava run and event headers
#include "AuxiliaryStructures.h"
#include "ALIBAVA.h"

// ROOT 
#include "TFile.h"
#include "TTree.h"

// System headers
#include <fstream>

// Some auxiliary functions ----------------------------
std::string trim_right(const std::string & s)
{
	const std::string b=" \t\n";
	std::string str = s;
	return str.erase(str.find_last_not_of(b) +1);
}

std::string trim_left(const std::string &s)
{
	std::string b=" \t\n";
	std::string str = s;
	return str.erase( 0, str.find_first_not_of(b) );
}

std::string trim_str(const std::string &s)
{
	std::string str = s;
	return trim_left(trim_right(str) );
}

// See alibava user manual
double tdc_time(unsigned int tdcTime)
{
    unsigned short fpart = tdcTime & 0xffff;
    short ipart = (tdcTime & 0xffff0000)>>16;
    if (ipart<0)
    {
        fpart *= -1;
    }
    //double tt = 100.*(1. -(ipart + (fpart/65535.)));
    double tt = 100.0*(ipart + (fpart/65535.));
    return tt;
}

// Coded temperature
double get_temperature(unsigned short temp)
{
    if(temp==0)
    {
        return 9999.;
    }
    else
    {
        return 0.12*temp - 39.8;
    }
}


// --
IOManager::IOManager(const std::string & rootfilename): 
    _runheader(nullptr), 
    _events(nullptr) 
{ 
    _file = new TFile(rootfilename.c_str(),"RECREATE"); 
}

IOManager::~IOManager()
{
    if(_runheader != nullptr)
    {
        delete _runheader;
        _runheader = nullptr;
    }
    if(_events != nullptr)
    {
        delete _events;
        _events = nullptr;
    }
}

void IOManager::book_tree_header()
{
    _tree_header = new TTree("runHeader","run header");

    _runheader   = new AlibavaRunHeader;
    _tree_header->Branch("version",&(_runheader->version));
    _tree_header->Branch("data_type",&(_runheader->data_type));
    _tree_header->Branch("run_number",&(_runheader->run_number));
    //_tree_header->Branch("header",&(_runheader.header),"I");
    //_tree_header->Branch("date_time",&(_runheader.header),"I");
    _tree_header->Branch("chipSelection",&(_runheader->chipSelection));
    _tree_header->Branch("header_pedestal",&(_runheader->header_pedestal));
    _tree_header->Branch("header_noise",&(_runheader->header_noise));
}

void IOManager::book_tree()
{
    _tree_events = new TTree("Events","Alibava events");

    _events   = new AlibavaEvent;
    _tree_events->Branch("type",&(_events->eventType));
    _tree_events->Branch("size",&(_events->eventSize));
    _tree_events->Branch("clock",&(_events->eventClock));
    _tree_events->Branch("runNumber",&(_events->runNumber));
    _tree_events->Branch("eventNumber",&(_events->eventNumber));
    _tree_events->Branch("eventTime",&(_events->eventTime));
    _tree_events->Branch("temperature",&(_events->eventTemp));
    _tree_events->Branch("calibration_charge",&(_events->calCharge));
    _tree_events->Branch("calibration_delay",&(_events->calDelay));
    _tree_events->Branch("chipheader_beetle1",&(_events->beetle1_chipheader));
    _tree_events->Branch("chipheader_beetle2",&(_events->beetle2_chipheader));
    _tree_events->Branch("data_beetle1",&(_events->beetle1_data));
    _tree_events->Branch("data_beetle2",&(_events->beetle2_data));
}


void IOManager::fill_header(const AlibavaRunHeader * aheader) const
{
    *_runheader = *aheader;
    _tree_header->Fill();
}

void IOManager::fill_event(const AlibavaEvent * anEvent) const
{
    *_events = *anEvent;
    _tree_events->Fill();
}

void IOManager::close()
{
    _tree_header->Write("", TTree::kOverwrite);
    _tree_events->Write("", TTree::kOverwrite);
    if(_file != nullptr)
    {
        _file->Close();
        delete _file;
        _file = nullptr;
    }
}

// Getter
TTree * IOManager::get_events_tree() const
{
    return this->_tree_events;
}

void IOManager::set_events_tree_access(const std::vector<std::string> & branch_list) const
{
    // XXX: No mechanism to check the validity of the branches (id they exist)
    // Speed up access, just using the branches we want
    this->_tree_events->SetBranchStatus("*",0);
    for(const auto & branch: branch_list)
    {
        this->_tree_events->SetBranchStatus(branch.c_str(),1);
    }
}

void IOManager::reset_events_tree() const
{
    this->_tree_events->ResetBranchAddresses();
    this->_tree_events->SetBranchStatus("*",1);
}

void IOManager::set_events_tree_branch_address(const std::string & branch_name, std::vector<float> ** v) const
{
    // XXX: No mechanism to check the validity of the branches (id they exist)
    this->_tree_events->SetBranchAddress(branch_name.c_str(),v);
}

int IOManager::get_events_number_entries() const
{
    return this->_tree_events->GetEntries();
}

void IOManager::get_events_entry(const int & i) const
{
    this->_tree_events->GetEntry(i);
}

