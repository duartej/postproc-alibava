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
#include "AlibavaDiagnosis.h"
#include "AuxiliaryStructures.h"
#include "AlibavaPostProcessor.h"
#include "ALIBAVA.h"

// ROOT 
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFriendElement.h"

// System headers
#include <fstream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <algorithm>
#include <ios>

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

// deallocate memory
template <class T> void deallocate_memory(std::map<int,std::vector<T>*> v)
{
    for(auto & i_v: v)
    {
        if(i_v.second != nullptr)
        {
            delete i_v.second;
            i_v.second = nullptr;
        }
    }
}

// clear and reserve memory before filling
template <class T> void clear_and_reserve_channels(std::map<int,std::vector<T>*> v)
{
    for(auto & i_v: v)
    {
        if( i_v.second != nullptr )
        {
            i_v.second->clear();
            i_v.second->reserve(ALIBAVA::NOOFCHANNELS);
        }
    }
}

// --
IOManager::IOManager(const std::string & rootfilename): 
    _rootfilename(rootfilename),
    _file(nullptr),
    _tree_header(nullptr),
    _tree_events(nullptr),
    _eventsProcessed(0),
    _channel_mask({{0,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)},
            {1,std::vector<int>(ALIBAVA::NOOFCHANNELS,1)}}),
    _monitor_plots_booked(false),
    _cal_parameters(nullptr),
    _runheader(nullptr), 
    _events(nullptr)
{ 
    _file = new TFile(_rootfilename.c_str(),"RECREATE");

    _monitor_plots = { {1,new AlibavaDiagnosis(1)}, {2,new AlibavaDiagnosis(2)} };
}


IOManager::~IOManager()
{
    // Store the monitor plots classes, for that, just open
    // the file first.
    // Check the file is closed, otherwise don't do anything
    if( _file != nullptr )
    {
        std::cerr << "[IOManager::~IOManager WARNING] Trying to update"
            << " an still open file. Please close it first. " << std::endl;
        std::cerr << "[IOManager::~IOManager WARNING] The monitor plots"
            << " are not going to be saved! " << std::endl;
    }
    else
    {
        _file = new TFile(_rootfilename.c_str(),"UPDATE"); 
        for(auto & mon: _monitor_plots)
        {
            if(mon.second != nullptr)
            {
                // Create and store the 
                mon.second->deliver_plots();
                delete mon.second;
                mon.second = nullptr;
            }
        }
    }
    // Closing the file
    this->close();
    
    // And deallocate the other members
    if(_cal_parameters != nullptr)
    {
        delete _cal_parameters;
        _cal_parameters = nullptr;
    }
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

void IOManager::set_calibration_parameters(const std::vector<int> & calparam_v)
{
    // DEPRECATED... NOT NEEDED ANYMORE
    _cal_parameters = new CalibrationParameters;
    _cal_parameters->nPulses        = calparam_v[0];
    _cal_parameters->initialCharge  = calparam_v[1];
    _cal_parameters->finalCharge    = calparam_v[2];
    _cal_parameters->deltaCharge    = calparam_v[3];
    // Note that the missing nSamplesPerPulse should be
    // updated once all events have been read
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

void IOManager::book_monitor_plots()
{
    // Initialize the monitor instances
    _monitor_plots_booked=true;
    for(const auto & mon: _monitor_plots)
    {
        mon.second->book_plots();
    }
}


void IOManager::fill_header(const AlibavaRunHeader * aheader) const
{
    *_runheader = *aheader;
    _tree_header->Fill();
}

void IOManager::fill_event(const AlibavaEvent * anEvent) 
{
    *_events = *anEvent;
    _tree_events->Fill();
    ++_eventsProcessed;
}

void IOManager::aux_store_friends(TTree * tree)
{
    // store the friends trees as well
    TIter next(tree->GetListOfFriends());
    TFriendElement * obj= nullptr;
    while((obj = dynamic_cast<TFriendElement*>(next())))
    {
        obj->GetTree()->Write("", TTree::kOverwrite);
    }
}

void IOManager::close()
{
    if( _tree_header != nullptr )
    {
        _tree_header->Write("", TTree::kOverwrite);
        this->aux_store_friends(_tree_header);
    }
    if( _tree_events != nullptr )
    {
        _tree_events->Write("", TTree::kOverwrite);
        this->aux_store_friends(_tree_events);
    }
    if(_file != nullptr)
    {
        _file->Close();
        delete _file;
        _file = nullptr;
        _tree_header = nullptr;
        _tree_events = nullptr;
    }
}

// Getter
TTree * IOManager::get_events_tree() const
{
    return this->_tree_events;
}

// Re-open the event trees
void IOManager::resurrect_events_tree()
{
    if( this->_tree_events != nullptr )
    {
        std::cerr << "[IOManager::resurrect_events_tree WARNING] Trying to"
            << " resurrect an still live TTree. Ignoring the order." << std::endl;
        return;
    }
    this->_tree_events = static_cast<TTree*>(this->_file->Get("Events"));
}

void IOManager::resurrect_header_tree()
{
    if( this->_tree_header != nullptr )
    {
        std::cerr << "[IOManager::resurrect_header_tree WARNING] Trying to"
            << " resurrect an still live TTree. Ignoring the order." << std::endl;
        return;
    }
    this->_tree_header = static_cast<TTree*>(this->_file->Get("runHeader"));
}

void IOManager::set_events_tree_access(const std::vector<std::string> & branch_list) const
{
    // XXX: No mechanism to check the validity of the branches (id they exist)
    // Speed up access, just using the branches we want
    this->_tree_events->SetBranchStatus("*",0);
    for(const auto & branch: branch_list)
    {
        // Be sure the string contains anything
        if(branch.empty())
        {
            continue;
        }
        this->_tree_events->SetBranchStatus(branch.c_str(),1);
    }
}

void IOManager::reset_events_tree() const
{
    this->_tree_events->ResetBranchAddresses();
    this->_tree_events->SetBranchStatus("*",1);
}

void IOManager::set_events_tree_branch_address(const std::string & branch_name, float * v) const
{
    // XXX: No mechanism to check the validity of the branches (id they exist)
    this->_tree_events->SetBranchAddress(branch_name.c_str(),v);
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

void IOManager::update(const CalibrateBeetleMap & cal_m)
{
    // Check the file is closed, otherwise don't do anything
    if( _file != nullptr )
    {
        std::cerr << "[IOManager::update WARNING] Trying to update"
            << " an still open file. Please close it first. " << std::endl;
        return;
    }

    _file = new TFile(_rootfilename.c_str(),"UPDATE"); 
    // Create a new Tree header 
    TTree * t = new TTree("postproc_runHeader","post-processed pedestals and common noise");
    
    std::map<int,std::vector<int>*> elec_adc = { {0, new std::vector<int>}, { 1, new std::vector<int>} }; 
    t->Branch("electronADC_beetle1",&(elec_adc[0]));
    t->Branch("electronADC_beetle2",&(elec_adc[1]));
    // The loop (first clean and reserve memory to fill the chip
    // vectors)
    clear_and_reserve_channels<int>(elec_adc);
    for(const auto & chip_p: cal_m)
    {
        for(int i = 0; i < static_cast<int>(chip_p.second.size()); ++i)
        {
            elec_adc[chip_p.first]->push_back(chip_p.second[i]);
        }
    }
    t->Fill();
    // Write it
    t->Write("", TTree::kOverwrite);    
    // Deallocate ttree (is already in the file)
    delete t;
    t = nullptr;
    
    // Deallocate memory
    deallocate_memory<int>(elec_adc);

    // Note that close will take care of closing the appropiate 
    // trees (events, runHeader)
    this->close();
}


void IOManager::update(const PedestalNoiseBeetleMap & pednoise_m)
{
    // Check the file is closed, otherwise don't do anything
    if( _file != nullptr )
    {
        std::cerr << "[IOManager::update WARNING] Trying to update"
            << " an still open file. Please close it first. " << std::endl;
        return;
    }
    
    _file = new TFile(_rootfilename.c_str(),"UPDATE"); 
    // Check if the tree already exist
    // in that case an extra branches will be place also in the Events
    bool is_calibrated = true;
    TTree *t = dynamic_cast<TTree*>(_file->Get("postproc_runHeader"));
    if( t == nullptr )
    {
        // Create a new Tree header and the Tree events
        t = new TTree("postproc_runHeader","post-processed pedestals and common noise");
        is_calibrated=false;
    }
    
    std::map<int,std::vector<float>*> pedestal = { {0, nullptr}, { 1, nullptr} }; 
    std::map<int,std::vector<float>*> noise = { {0,nullptr}, {1,nullptr} }; 
    std::map<int,std::vector<int>*> channel_mask = { {0,nullptr}, {1,nullptr} }; 
    // Filling it as branches: 
    // XXX: BE CAREFUL though. If storage problems show up maybe change the approach to
    // create a new tree
    // See. https://root.cern.ch/doc/master/classTTree.html (Adding a Branch to a Existing Tree)
    std::vector<TBranch*> newbranches;
    newbranches.push_back(t->Branch("pedestal_cmmd_beetle1",&(pedestal[0])));
    newbranches.push_back(t->Branch("noise_cmmd_beetle1",&(noise[0])));
    newbranches.push_back(t->Branch("channel_mask_beetle1",&(channel_mask[0])));
    newbranches.push_back(t->Branch("pedestal_cmmd_beetle2",&(pedestal[1])));
    newbranches.push_back(t->Branch("noise_cmmd_beetle2",&(noise[1])));
    newbranches.push_back(t->Branch("channel_mask_beetle2",&(channel_mask[1])));

    // The loop (clear the vectors, and reserve memory)
    clear_and_reserve_channels<float>(pedestal);
    clear_and_reserve_channels<float>(noise);
    clear_and_reserve_channels<int>(channel_mask);
    for(const auto & chip_p: pednoise_m)
    {
        for(int i = 0; i < static_cast<int>(chip_p.second.first.size()); ++i)
        {
            pedestal[chip_p.first]->push_back(chip_p.second.first[i]);
            noise[chip_p.first]->push_back(chip_p.second.second[i]);
            // if pedestal and noise are 0 --> the channel is masked
            if(chip_p.second.first[i] < 1e-9 && chip_p.second.second[i] < 1e-9)
            {
                channel_mask[chip_p.first]->push_back(0);
            }
            else
            {
                channel_mask[chip_p.first]->push_back(1);
            }
        }
    }
    // fill all branches
    std::for_each(newbranches.begin(),newbranches.end(), [] (TBranch* br) { br->Fill(); });
    // Write it
    t->Write("", TTree::kOverwrite);
    
    // obtain the calibrated vector if exist
    std::map<int,std::vector<float>*> * elec_per_adc = nullptr;
    if(is_calibrated)
    {
        std::map<int,std::vector<float>*> cal_m;
        t->SetBranchAddress("electronADC_beetle1",&(cal_m[0]));
        t->SetBranchAddress("electronADC_beetle2",&(cal_m[1]));
        t->GetEntry(0);
        // and copy 
        elec_per_adc= new std::map<int,std::vector<float>*>;
        (*elec_per_adc)[0] = cal_m[0];
        (*elec_per_adc)[1] = cal_m[1];
    }
    // Tree not needed anymore, 
    delete t;
    t = nullptr;

    
    // And the events tree: I need to ressurrect the Events tree
    TTree * tevt = new TTree("postproc_Events","post-processed signals");
    std::map<int,std::vector<float>* > data = { {0, new std::vector<float>}, { 1, new std::vector<float>} }; 
    tevt->Branch("postproc_data_beetle1",&(data[0]));
    tevt->Branch("postproc_data_beetle2",&(data[1]));
    std::map<int,float> common_mode = { {0,-9999.9}, { 1, -9999.9} }; 
    tevt->Branch("postproc_cmmd_beetle1",&(common_mode[0]));
    tevt->Branch("postproc_cmmd_beetle2",&(common_mode[1]));
    std::map<int,float> noise_cmmd = { {0, -9999.0}, { 1, -9999.0} }; 
    tevt->Branch("postproc_cmmd_noise_beetle1",&(noise_cmmd[0]));
    tevt->Branch("postproc_cmmd_noise_beetle2",&(noise_cmmd[1]));
    // The calibrated if needed
    //
    std::map<int,std::vector<float>* > out_cal = { {0,new std::vector<float>}, {1,new std::vector<float>} };
    // And define the setter to that collection (a dummy functor if no calibration file)
    std::function<void(const int &, const int &)> fill_calibrated_branches = [&] (const int & /*chip*/, const int & /*istrip*/) { return; };
    if(is_calibrated)
    {
        tevt->Branch("postproc_cal_data_beetle1",&(out_cal[0]));
        tevt->Branch("postproc_cal_data_beetle2",&(out_cal[1]));
        fill_calibrated_branches = [&] (const int & chip, const int & istrip) 
                { out_cal[chip]->push_back( (*(*elec_per_adc)[chip])[istrip] * (*data[chip])[istrip] ) ; };
    }

    // Resurrect the Events tree to extract the data signal, in order 
    // to correct it
    this->resurrect_events_tree();
    // Activate the needed branches
    std::map<int,std::vector<float>* > original_data = { {0, nullptr}, { 1, nullptr} };
    std::map<int,std::string> chip_dataname_map = { {0,"data_beetle1"}, {1,"data_beetle2"} };
    
    set_events_tree_access( {chip_dataname_map[0], chip_dataname_map[1]} );

    // attach the vector to the branches 
    for(auto & i_v: original_data)
    {
        set_events_tree_branch_address(chip_dataname_map[i_v.first],&i_v.second);
    }
    
    // Go back to the previous position, move up 1 line, and set 
    // the 80 column
    std::cout << "\033[1A\033[45C[000%]" << std::flush;
    const int nentries= this->get_events_number_entries();
    float point = float(nentries)/100.0;
    // the event loop
    for(int k=0; k < this->get_events_number_entries(); ++k)
    {
        std::cout << std::fixed << std::setprecision(0) << "\r\033[45C[" << std::setfill('0') << std::setw(3) 
            << std::round(float(k)/point) << "%]" << std::flush;
        clear_and_reserve_channels(data);
        clear_and_reserve_channels(out_cal);
        
        this->get_events_entry(k);

        std::vector<float> pedestal_free;
        pedestal_free.reserve(ALIBAVA::NOOFCHANNELS);
        // And update using the pedestal and noise corrections
        for(const auto & chip_rawdata: original_data)
        {
            pedestal_free.clear();
            pedestal_free.reserve(ALIBAVA::NOOFCHANNELS);
            // Get the common mode for that chip, it needs the pedestal-free
            for(int ch = 0; ch < static_cast<int>(chip_rawdata.second->size()); ++ch)
            {
                pedestal_free.push_back( (*chip_rawdata.second)[ch]-(*pedestal[chip_rawdata.first])[ch] );
            }
            const auto cmmd_and_noise = AlibavaPostProcessor::calculate_common_noise(pedestal_free,*(channel_mask[chip_rawdata.first]));
            // And also store it
            common_mode[chip_rawdata.first]=cmmd_and_noise.first;
            noise_cmmd[chip_rawdata.first] =cmmd_and_noise.second;
            // Subtract noise and pedestals to the raw-data
            for(int ichan = 0 ; ichan < static_cast<int>(chip_rawdata.second->size()); ++ichan)
            {
                // Signal
                if((*channel_mask[chip_rawdata.first])[ichan] == 0)
                {
                    data[chip_rawdata.first]->push_back( 0.0 );
                    //fill_calibrated_branches(chip_rawdata.first,ichan);
                }
                else
                {
                    data[chip_rawdata.first]->push_back( pedestal_free[ichan] - cmmd_and_noise.first );
                    fill_calibrated_branches(chip_rawdata.first,ichan);
                }
            }
        }
        tevt->Fill();
    }
    std::cout << std::endl; 
    // Write it
    tevt->Write("", TTree::kOverwrite);    
    // Deallocate
    delete tevt;
    tevt = nullptr;

    // Deallocate memory
    deallocate_memory<float>(data);
    deallocate_memory<float>(out_cal);
    // the pedestal and noise
    deallocate_memory<float>(pedestal);
    deallocate_memory<float>(noise);

    // Everything back to the original state
    this->reset_events_tree();

    // Note that close will take care of closing the appropiate 
    // trees (events, runHeader)
    this->close();
}


void IOManager::book_monitor_plot(const std::string & plotname, const TObject * theplot, const int & chip)
{
    if(_monitor_plots.find(chip) == _monitor_plots.end())
    {
        std::cerr << "[IOManager::book_monitor_plot ERROR] Invalid chip"
            << " number [" << chip << "] " << std::endl;
        // Exception??
        return;
    }
    // Booking it in the proper monitor
    _monitor_plots[chip]->book_plot(plotname,theplot);
}

const std::vector<TObject*> IOManager::get_calibration_objects(const int & chipnumber) const
{
    return _monitor_plots.at(chipnumber)->get_calibration_plots();
}

void IOManager::set_calibration_plot(const IOManager & cal_manager)
{
    for(auto & chipmon: _monitor_plots)
    {
        // Choose the right monitor manager and send the plots
        chipmon.second->set_calibration_plot(cal_manager.get_calibration_objects(chipmon.first));
    }
}

void IOManager::fill_diagnostic_plots()
{
    bool file_was_closed = false;
    // Check the file is closed to close it afterwards
    if( _file == nullptr )
    {
        file_was_closed = true;
        _file = new TFile(_rootfilename.c_str(),"UPDATE"); 
        this->resurrect_events_tree();
        this->resurrect_header_tree();
    }
    // Check if the pedefile was processed, if yes it should exist 
    // the postproc_Events tree
    TTree *event_tree = dynamic_cast<TTree*>(_file->Get("postproc_Events"));
    TTree * header_tree = nullptr;
    if( event_tree == nullptr )
    {
        // Not postproc available, ressurrect the event tree
        // and assign as the event tree
        event_tree   = _tree_events;
        header_tree = _tree_header;
    }
    else
    {
        // Then it must exist the postproc header as well
        header_tree = dynamic_cast<TTree*>(_file->Get("postproc_runHeader"));
        // Add it as friend the regular Event tree as well
        event_tree->AddFriend(_tree_events);
    }
    // Add it as friend in order to perform whatever manipulations
    // needed to do
    event_tree->AddFriend(header_tree);

    // Call the plotters
    for(auto & chipmon: _monitor_plots)
    {
        chipmon.second->fill_diagnostic_plots(event_tree,header_tree);
    }

    // Everything back to the original state
    this->reset_events_tree();
    
    if(file_was_closed)
    {
        // Note that close will take care of closing the appropiate     
        // trees as well (events, runHeader)
        this->close();
    }
}
