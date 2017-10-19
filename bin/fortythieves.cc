// ***************************************************************************
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

// system header
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <memory>
#include <map>

// Define alibava header and events (not using ROOT dictionaries
// although it can be think to use it later)
#include "AuxiliaryStructures.h"

// Alibava reader and post-processing
#include "IOAlibavaReader.h"
#include "AlibavaPostProcessor.h"

// ROOT dedicated classes
#include "IOManager.h"

// The input options
#include "input_options.h"

// Just summary of error codes
void print_error_codes()
{
    std::cout << "*************** Error codes ***************" << std::endl;
    std::cout << "0: Succesful execution" << std::endl;
    std::cout << "-1: Invalid number of input arguments" << std::endl;
    std::cout << "-2: Undefined input option" << std::endl;
    std::cout << "-3: Invalid raw input binary file" << std::endl;
    std::cout << "-4: Invalid firmware version from raw data" << std::endl;
    std::cout << "-5: Invalid data type (user)" << std::endl;
    std::cout << "-6: Inconsistent Run Header in calibration file" << std::endl;
}

// Helper function to extract and parse the content of the -u option
int parse_use_channels_opt(std::string raw_str,std::map<int,std::vector<int> > & parsed_d)
{
    // Add a delimeter to facilitate the parsing
    raw_str.append(",");
    // Expected format is chip:r1-r2,chip:r3-r4,...
    const std::string delim_set(",");
    const std::string delim_chip(":");
    const std::string delim_range("-");
    // The loop to find the sets of range1,range2,..
    size_t pos = 0;
    while( (pos=raw_str.find(delim_set)) != std::string::npos )
    {
        /*if( pos != raw_str.length()-1 && raw_str[pos] != delim_set[0] )  
        {
            continue;
        }*/
        std::string element(raw_str.substr(0,pos));
        raw_str.erase(0,pos+1);
        // Now evaluate the string "chip:r1-r2"
        const size_t pos_el = element.find(delim_chip);
        const int chip = std::stoi(element.substr(0,pos));
        if( parsed_d.count(chip) != 1 )
        {
            std::cerr << "\033[1;31mfortythieves ERROR\033[1;m "  
                << "Invalid chip number: " << chip << std::endl;
            return -7;
        }
        const size_t pos_ranges = element.find(delim_range);
        const int rlow  = std::stoi(element.substr(pos_el+1,pos_ranges));
        const int rhigh = std::stoi(element.substr(pos_ranges+1));
        if( rlow >= rhigh )
        {
            std::cerr << "\033[1;31mfortythieves ERROR\033[1;m "
                << "Invalid channel ranges: " << rlow << "-" << rhigh << std::endl;
            return -8;
        }
        // XXX error control!!
        // Now create an xrange equivalen
        for(int i=rlow; i <= rhigh;++i)
        {
            parsed_d[chip].push_back(i);
        }
    }
    return 0;
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
      << " -p alibava raw-binary containing the pedestal run\n"
      << " -c alibava raw-binary containing the calibration run\n"
      << " -r run number [-1]\n"
      << " -u use channels introduced as <chip>:ch0-ch1,<chip>:ch2-ch3,...\n "
      << "   Those channels not in the list are masked\n"
      << " -n no automask noisy channels\n"
      << " -h show this help" << std::endl;
}

int main(int argc, char* argv[]) 
{
    // Check that correct number of command-line arguments
    if(argc < 2 && std::string(argv[1]) != "-h") 
    {
        std::cerr << " Unexpected number of command-line arguments. \n You are"
            << " expected to provide one input file name. \n"
            << " Program stopped! " << std::endl;
        return -1;
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
            opt.pedestal_file = argv[i+1];
            opt.storeHeaderPedestalNoise = true;
            ++i;
        }
        else if( strcmp(argv[i],"-c") == 0 )
        {
            opt.calibration_file = argv[i+1];
            ++i;
        }
        else if( strcmp(argv[i],"-r") == 0 )
        {
            opt.runNumber = std::stoi(argv[i+1]);
            ++i;
        }
        else if( strcmp(argv[i],"-u") == 0 )
        {
            const int return_value=parse_use_channels_opt(argv[i+1],opt.use_channels);
            if(return_value != 0)
            {
                return return_value;
            }
            ++i;
        }
        else if( strcmp(argv[i],"-n") == 0 )
        {
            opt.no_automasking = true;
        }
        else
        {
            // Check that the provided input name corresponds to an existing file.
            std::ifstream is(argv[i]);
            if(!is && std::string(argv[i]) != "-h") 
            {
                std::cerr << " Command-line file '" << argv[i] << "' was not found. \n"
                    << " Program stopped! " << std::endl;
                return -2;
            }
            opt.cmndfile = argv[i];
        }
    }

    // Initialize and prepare the ROOT output
    IOManager iomanager(opt.outputFilename);
    iomanager.book_tree_header();
    iomanager.book_tree();
    // Monitor plots, booking
    iomanager.book_monitor_plots();

    // process the file
    const int events = IOAlibavaReader::read_data(opt,&iomanager);
    std::cout << "Processed " << events << std::endl;
    int status = 0; 
    if(events < 0)
    {
        // A event < 0 is actually an error code, something went 
        // wrong
        status = events;
    }
    iomanager.close();
    
    // process calibration file
    if(opt.calibration_file != "")
    {
        std::cout << "\033[1;34mfortythieves INFO\033[1;m: " 
            << "Calibration actived" << std::endl;
        // The name of the output ROOT file
        const auto dotpos = opt.outputFilename.find(".root");
        std::string calfile(opt.outputFilename.replace(dotpos,5,"_cal.root"));
        
        // Initialize and prepare the output
        IOManager iomanager_cal(calfile);
        iomanager_cal.book_tree_header();
        iomanager_cal.book_tree();
        // process the pedestal file
        // change the name of the input file
        input_options opt_cal(opt);
        opt_cal.cmndfile = opt.calibration_file;
        const int events_cal = IOAlibavaReader::read_data(opt_cal,&iomanager_cal);
        if(events_cal < 0)
        {
            status += events_cal;
        }
        // calibrating
        // XXX: Create an unique postproc
        AlibavaPostProcessor postproc;
        std::cout << " - Calibrating" << std::endl;
        CalibrateBeetleMap cal_map = postproc.calibrate(iomanager_cal);
        // Extract the summary plot for calibration and set it in the main manager
        iomanager.set_calibration_plot(iomanager_cal);
        // close the calibration file
        iomanager_cal.close();
        // And update the beam file with the calibration vector
        // included in the runHeader postproc 
        iomanager.update(cal_map);
    }

    // process pedestal file
    if(opt.storeHeaderPedestalNoise)
    {
        std::cout << "\033[1;34mfortythieves INFO\033[1;m: " 
            << "Pedestal and noise calculations actived" << std::endl;
        // The name of the output ROOT file
        const auto dotpos = opt.outputFilename.find(".root");
        std::string pedfile(opt.outputFilename.replace(dotpos,5,"_ped.root"));
        // Initialize and prepare the output
        IOManager iomanager_ped(pedfile);
        iomanager_ped.book_tree_header();
        iomanager_ped.book_tree();
        // process the pedestal file
        // change the name of the input file
        input_options opt_ped(opt);
        opt_ped.cmndfile = opt.pedestal_file;
        const int events_ped= IOAlibavaReader::read_data(opt_ped,&iomanager_ped);
        if(events_ped < 0)
        {
            status += events_ped;
        }
    
        // calculate pedestals and common noise: { bettle: { channels ,,, } }
        // XXX: Create an unique postproc
        AlibavaPostProcessor postproc(opt.use_channels);
        std::cout << " - Calculating <pedestals>" << std::endl;
        PedestalNoiseBeetleMap pednoise_cmmdnot = postproc.calculate_pedestal_noise(iomanager_ped);
        postproc.print_channels_mask();
        
        std::cout << " - Re-evaluating pedestals (noise extracted) and common mode" << std::endl;
        postproc.get_pedestal_noise_free(iomanager_ped,pednoise_cmmdnot);
        
        std::cout << " - Re-calculating pedestals and noise per channel" << std::endl;
        PedestalNoiseBeetleMap pednoise_cmmd = postproc.calculate_pedestal_noise(iomanager_ped);
        // close the pedestal noise file
        iomanager_ped.close();
        // And update the beam file with the pedestals and common noise values
        // included in the runHeader postproc 
        std::cout << " - Updating the tree: create signal branches" << std::endl;
        iomanager.update(pednoise_cmmd);
    }
    // process the diagnostic plots
    std::cout << "\033[1;34mfortythieves INFO\033[1;m: " 
        << "Filling monitor plots CHIP" << std::endl;
    iomanager.fill_diagnostic_plots();

    return status;
}


