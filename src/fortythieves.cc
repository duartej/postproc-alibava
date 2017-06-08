// ***********************************************
// 
// Alibava binary raw data convertor into ROOT
//
// See the section 'The Alibava data format' at 
// https://www.alibavasystems.com/images/Catalogo/alibava-usermanual.pdf
// 
// An alternative to extract the binary data can
// be found at the AsciiRoot.cc code inside the
// root_macros folder of the alibava code:
// https://www.alibavasystems.com/downloads-alibava-systems/alibava-classic-downloads.html
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// (based on the Marlin processor, AlibavaCreator, 
// created by Eda Yildirim
// https://github.com/eutelescope/eutelescope/blob/v1.0-tag/include/alibava/AlibavaConverter.h )
// ****************************************************************************

#include "ALIBAVA.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

// system header
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <memory>

// Define alibava header and events (not using ROOT dictionaries
// although it can be think to use it later)
struct AlibavaRunHeader
{
    AlibavaRunHeader() 
    {
        version = -1;
        data_type = -1;
        run_number = -1;
        header = "";
        date_time = "";
        chipSelection.push_back(0);
        chipSelection.push_back(1);
    }
    int version;
    int data_type;
    int run_number;
    std::string header;
    std::string date_time;
    std::vector<int> chipSelection;
    std::vector<float> header_pedestal;
    std::vector<float> header_noise;
};

struct AlibavaEvent
{
    AlibavaEvent() 
    { 
        runNumber = -1;
        eventNumber = -1;
        eventType   = 9999;
        eventSize   = 9999;
        eventClock  = 9999;
        eventTime = -1.0;
        eventTemp = -1.0;
        calCharge = -9999;
        calDelay  = -9999;

        beetle1_chipheader.reserve(ALIBAVA::NOOFCHANNELS);
        beetle2_chipheader.reserve(ALIBAVA::NOOFCHANNELS);

        beetle1_data.reserve(ALIBAVA::NOOFCHANNELS);
        beetle2_data.reserve(ALIBAVA::NOOFCHANNELS);
    }
    unsigned int eventType ;
    unsigned int eventSize ;
    unsigned int eventClock;
    int runNumber;
    int eventNumber;
    float eventTime ;
    float eventTemp ;
    float calCharge ;
    float calDelay  ;

    std::vector<float> beetle1_chipheader;
    std::vector<float> beetle2_chipheader;
    
    std::vector<float> beetle1_data;
    std::vector<float> beetle2_data;
};

// ROOT dedicated classes
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
        IOManager(const std::string & rootfilename): 
            _runheader(nullptr), 
            _events(nullptr) { _file = new TFile(rootfilename.c_str(),"RECREATE"); }
        ~IOManager();

        void book_tree_header();
        void book_tree();
        
        void fill_header(const AlibavaRunHeader * aheader) const;
        void fill_event(const AlibavaEvent * anAlibavaEvent) const;
        
        void close();
};

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
//---------

// Just summary of error codes
void print_error_codes()
{
    std::cout << "*************** Error codes ***************" << std::endl;
    std::cout << "0: Succesful execution" << std::endl;
    std::cout << "1: Invalid number of input arguments" << std::endl;
    std::cout << "2: Undefined input option" << std::endl;
    std::cout << "3: Invalid raw input binary file" << std::endl;
    std::cout << "4: Invalid firmware version from raw data" << std::endl;
    std::cout << "5: Invalid data type (user)" << std::endl;
}

// Helper structure to fill all the input options
struct input_options
{
    input_options() 
    {
        // Defaults
        outputFilename = "fortythieves.root";
        storeHeaderPedestalNoise = false;
        runNumber = -1;
        cmndfile  = "";
        startEventNum = -1;
        stopEventNum  = -1;
    }
    std::string cmndfile;
    std::string outputFilename;
    int runNumber;
    int startEventNum;
    int stopEventNum;
    bool storeHeaderPedestalNoise;
};

// Some auxiliary functions
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



// Input parser functions
void display_usage()
{
    std::cout << "\033[37musage:\033[m fortythieves [OPTIONS] alibava_data.raw"
        << std::endl;
    std::cout << std::endl;
    std::cout << "Convert a raw binary data from the ALIBAVA DAQ into a ROOT file" 
        << std::endl;
    std::cout << std::endl;
    std::cout << "[OPTIONS]\n -o name of the ROOT output file [fortythieves.root]\n"
      << " -p flag to store pedestal and noise header [false]\n"
      << " -r run number [-1]\n"
      /*<< " -m mis-identification probability, a value different from 0 will force '-t pions_kaons'"
      << " regardless of the user input [default: 0.0]\n"
      */
      << " -h show this help" << std::endl;
}

// Declaration
int readData(const input_options & opt, const IOManager & iomanager);

int main(int argc, char* argv[]) 
{
    // Check that correct number of command-line arguments
    if(argc < 2 && std::string(argv[1]) != "-h") 
    {
        std::cerr << " Unexpected number of command-line arguments. \n You are"
            << " expected to provide one input file name. \n"
            << " Program stopped! " << std::endl;
        return 1;
    }

    // The options 
    input_options opt;

    // get options
    for(int i = 1; i < argc; ++i)
    {
        if( strcmp(argv[i],"-h") == 0 )
        {
            display_usage();
            return 0;
        }
        else if( strcmp(argv[i],"-o") == 0 )
        {
            opt.outputFilename = argv[i+1];
            ++i;
        }
        else if( strcmp(argv[i],"-p") == 0 )
        {
            opt.storeHeaderPedestalNoise = true;
        }
        else if( strcmp(argv[i],"-r") == 0 )
        {
            opt.runNumber = std::stoi(argv[i+1]);
            ++i;
        }
        /*
        else if( strcmp(argv[i],"-m") == 0 )
        {
            misid_ratio = std::stof(argv[i+1]);
            ++i;
        }
        */
        else
        {
            // Check that the provided input name corresponds to an existing file.
            std::ifstream is(argv[i]);
            if(!is && std::string(argv[i]) != "-h") 
            {
                std::cerr << " Command-line file '" << argv[i] << "' was not found. \n"
                    << " Program stopped! " << std::endl;
                return 2;
            }
            opt.cmndfile = argv[i];
        }
    }

    // Initialize and prepare the ROOT output
    IOManager iomanager(opt.outputFilename);
    iomanager.book_tree_header();
    iomanager.book_tree();

    // process the file
    int status = readData(opt,iomanager);

    // close output ROOT file
    iomanager.close();

    return status;
}


int readData(const input_options & opt, const IOManager & iomanager) 
{
    // this event counter is used to stop the processing when it is
    // greater than numEvents.
    int eventCounter = 0;
    
    // print out a debug message
    std::cout << "fortythieves: Reading " << opt.cmndfile << ", run number: " << opt.runNumber << std::endl;
    
    ////////////////
    // Open File  //
    ////////////////
    std::ifstream infile;
    infile.open(opt.cmndfile);
    if(!infile.is_open()) 
    {
        std::cerr << "'fortythieves' could not read the file "
            <<opt.cmndfile<<" correctly. Please check the path and file names that" 
            << " have been input" << std::endl;
        return 3;
    }
    else
    {
        std::cout<<"Input file "<<opt.cmndfile<<" is now opened"<<std::endl;
    }
    /////////////////
    // Read Header //
    /////////////////
    // --- Time of start of the run 
    unsigned int date = 0;
    infile.read(reinterpret_cast<char*>(&date), sizeof(unsigned int));
    // --- Run type: 1:Calibration; 2:Laser Syn.; 3:Laser; 4: Rad. source; 5: Pedestal
    int type = -1;
    infile.read(reinterpret_cast<char*>(&type), sizeof(int));
    // --- Header lenght
    unsigned int lheader = 0; 
    infile.read(reinterpret_cast<char*>(&lheader), sizeof(unsigned int));
    // read the next field until reach the lenght of the header
    std::string header("");
    header.clear();
    for(unsigned int ic=0; ic<lheader; ++ic)
    {
        char tmp_c;
        infile.read(&tmp_c, sizeof(char));
        header.append(1, tmp_c);
    }
    header = trim_str(header);
    
    // Firmware version
    int version(-1);
    if (header[0]!='V' && header[0]!='v')
    {
        version = 0;
    }
    else
    {
        version = int(header[1]-'0');
        header = header.substr(5);
    }
    /// XXX: Different treatment when calibration run?
    ////////////////////////////////////
    // Read header pedestal and noise //
    // /////////////////////////////////
	
    // Alibava stores a pedestal and noise set in the run header. 
    // These values are not used in te rest of the analysis, so it is 
    // optional to store it. By default it will not be stored, but it you
    // want you can set _storeHeaderPedestalNoise variable to true.
    float tmp_float=-11;
    std::vector<float> headerPedestal;
    // first pedestal
    for(int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ++ichan) 
    {
        infile.read(reinterpret_cast<char*> (&tmp_float), sizeof(double));
        headerPedestal.push_back(tmp_float);
    }	
    // now noise
    std::vector<float> headerNoise;
    for(int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ichan++) 
    {
        infile.read(reinterpret_cast<char*>(&tmp_float), sizeof(double));
        headerNoise.push_back(tmp_float);
    }
    ////////////////////
    // Process Header //
    ////////////////////
    AlibavaRunHeader * runHeader = new AlibavaRunHeader;
    runHeader->header = header;
    runHeader->version= version;
    runHeader->data_type =type;
    runHeader->date_time = std::string(ctime(reinterpret_cast<time_t*>(&date)));
    if(opt.storeHeaderPedestalNoise) 
    {
        runHeader->header_pedestal = headerPedestal;
        runHeader->header_noise = headerNoise;
    }
    runHeader->run_number = opt.runNumber;
    
    // Store the run header
    iomanager.fill_header(runHeader);
    delete runHeader;
    runHeader = nullptr;

    ////////////////
    // Read Event //
    ////////////////
    // Expected different streams depending the firmware version
    // used to store the data
    if(version<2) 
    {
        // this code is not written for version<=1.
        std::cerr <<" Not supported data version found (version="
            <<version<<"<2). Data is not saved."<<std::endl;
	return 4;
    }
    
    // Event loop
    while( !(infile.bad() || infile.eof()) )
    {
        if( eventCounter % 1000 == 0 )
        {
            std::cout << "\r Processing  "<< eventCounter << " from run " 
                << opt.runNumber << std::flush;
        }
	
        unsigned int headerCode = 0;
        unsigned int eventTypeCode=0;
        do
        {
            infile.read(reinterpret_cast< char *> (&headerCode), sizeof(unsigned int));
            if(infile.bad() || infile.eof())
            {
                std::cout << std::endl;
                return 0;
            }
	    eventTypeCode = (headerCode>>16) & 0xFFFF;
        } while( eventTypeCode != 0xCAFE );
	
        eventTypeCode = headerCode & 0x0FFF;

        unsigned int userEventTypeCode = headerCode & 0x1000;
        if(userEventTypeCode)
        {	
            std::cout << std::endl;
            std::cout<<"Unexpected data type (user type). Data is not saved"<<std::endl;
            return 5;
        }
        unsigned int eventSize = 0;
        infile.read(reinterpret_cast< char *> (&eventSize), sizeof(unsigned int));
        
        double value = -1;
        infile.read(reinterpret_cast< char *> (&value), sizeof(double));
		
	//see AlibavaGUI.cc
        double charge = int(value) & 0xFF;
        double delay = int(value) >> 16;
        charge = charge * 1024;
	
        // The timestamp
        unsigned int clock = 0;
        // Thomas 13.05.2015: Firmware 3 introduces the clock to the header!
        // for now this is not stored...
        if(version==3)
        {
            infile.read(reinterpret_cast< char *> (&clock), sizeof(unsigned int));
        }
        // Time digital converter 
        unsigned int tdcTime = 0;
        infile.read(reinterpret_cast< char *> (&tdcTime), sizeof(unsigned int));
        // Temperature measured on daughter board
	unsigned short temp = 0;      
        infile.read(reinterpret_cast< char *> (&temp), sizeof(unsigned short));
        
        // vector for data
        std::map<int,std::vector<float> > beetles_data; // = { {0,{}}, {1,{}} };
        beetles_data[0].reserve(ALIBAVA::NOOFCHANNELS);
        beetles_data[1].reserve(ALIBAVA::NOOFCHANNELS);
        // vector for chip header
        std::map<int,std::vector<float> > beetles_chipheaders; // = { {0.{}}, {1,{}} };
        beetles_chipheaders[0].reserve(ALIBAVA::CHIPHEADERLENGTH);
        beetles_chipheaders[1].reserve(ALIBAVA::CHIPHEADERLENGTH);
        
        // An auxiliary variable (contains the data before pushing in the 
        // previous vector
        short tmp_short;
        // Chip header
        unsigned short chipHeader[2][ALIBAVA::CHIPHEADERLENGTH];
        // iterate over number of chips and store the chip headers 
        // and data
        for(int ichip=0; ichip < ALIBAVA::NOOFCHIPS; ++ichip)
        {
            infile.read(reinterpret_cast< char *>(chipHeader[ichip]), ALIBAVA::CHIPHEADERLENGTH*sizeof(unsigned short));
            // store chip header in all_chipheaders vector
            for(int j = 0; j < ALIBAVA::CHIPHEADERLENGTH; ++j)
            {
                beetles_chipheaders[ichip].push_back(float(chipHeader[ichip][j]));
            }
            // store data in all_data vector
	    for(int ichan=0; ichan < ALIBAVA::NOOFCHANNELS; ichan++) 
            {
                infile.read(reinterpret_cast< char *> (&tmp_short), sizeof(unsigned short));
                beetles_data[ichip].push_back(float(tmp_short));
            }
        }
        // -- so far the event is read
        ++eventCounter;
        
        // Now store it in the new format
        //
        // Skip event 
        if(opt.startEventNum!=-1 && eventCounter == opt.startEventNum) 
        {
            std::cout << std::endl;
        }
        if(opt.startEventNum!=-1 && eventCounter < opt.startEventNum) 
        {
            std::cout << "\r Skipping event "<< eventCounter <<". StartEventNum is set to "
                <<opt.startEventNum << std::flush;
            if(eventCounter == opt.startEventNum)
            {
                std::cout << std::endl;
            }
            continue;
        }
        
        if(opt.stopEventNum!=-1 && eventCounter >= opt.stopEventNum) 
        {
            std::cout << std::endl;
            std::cout <<"fortythieves: Reached StopEventNum '"
                << opt.stopEventNum <<"'. Last saved event number is "
                << eventCounter-1 << std::endl;
            break;
        }
        // Store it
        AlibavaEvent* anEvent = new AlibavaEvent();
        anEvent->runNumber = opt.runNumber;
        anEvent->eventNumber = eventCounter;
        anEvent->eventType   = eventTypeCode;
        anEvent->eventSize   = eventSize;
        if(version==3)
        {
            anEvent->eventClock = clock;
        }
        anEvent->eventTime = static_cast<float>(tdc_time(tdcTime));
	anEvent->eventTemp = static_cast<float>(get_temperature(temp));
	anEvent->calCharge = charge;
	anEvent->calDelay  = delay;
        
        anEvent->beetle1_chipheader = beetles_chipheaders[0];
        anEvent->beetle2_chipheader = beetles_chipheaders[1];

        anEvent->beetle1_data = beetles_data[0];
        anEvent->beetle2_data = beetles_data[1];
        
        // and store it 
        iomanager.fill_event(anEvent);

        // Free memory
        delete anEvent;
    }
    std::cout << std::endl;
    
    infile.close();
    if(opt.stopEventNum!=-1 && eventCounter<opt.stopEventNum)
    {
        std::cout <<"fortythieves WARNING: Stopped before reaching"
            << " StopEventNum: " <<opt.stopEventNum<<". The file has "
            << eventCounter << " events."<<std::endl;
    }
    return 0;
}

