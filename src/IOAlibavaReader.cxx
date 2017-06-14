// ************************************************
// 
// Manager to read binary raw data from the ALIBAVA
// DAQs
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "IOAlibavaReader.h"

// The alibava run and event headers
#include "AuxiliaryStructures.h"
#include "ALIBAVA.h"

// System headers
#include <fstream>
#include <map>

// Some auxiliary functions ----------------------------
namespace auxfunc
{
    std::string path_filename(const std::string & fullpath)
    {
        size_t pos = fullpath.find_last_of("/\\");
        return fullpath.substr(pos+1);
    }
    
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
}

// See alibava user manual
double IOAlibavaReader::tdc_time(const unsigned int & tdcTime)
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
double IOAlibavaReader::get_temperature(const unsigned short & temp)
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


int IOAlibavaReader::read_data(const input_options & opt,const IOManager & iomanager) 
{
    // this event counter is used to stop the processing when it is
    // greater than numEvents.
    int eventCounter = 0;
    
    // print out a debug message
    std::cout << "IOAlibavaReader.read_data: Reading '" 
        <<  auxfunc::path_filename(opt.cmndfile) << "', run number: " << opt.runNumber << std::endl;
    
    ////////////////
    // Open File  //
    ////////////////
    std::ifstream infile;
    infile.open(opt.cmndfile);
    if(!infile.is_open()) 
    {
        std::cerr << "\033[1;31mIOAlibavaReader.read_data\033[1;m could not read the file "
            << auxfunc::path_filename(opt.cmndfile)<<" correctly. Please check the path and file names that" 
            << " have been input" << std::endl;
        return 3;
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
    header = auxfunc::trim_str(header);
    
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
            std::cout << "\rProcessing  "<< eventCounter << " from run " 
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
        anEvent->eventTime = static_cast<float>(IOAlibavaReader::tdc_time(tdcTime));
	anEvent->eventTemp = static_cast<float>(IOAlibavaReader::get_temperature(temp));
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


