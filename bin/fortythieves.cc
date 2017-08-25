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
      /*<< " -m mis-identification probability, a value different from 0 will force '-t pions_kaons'"
      << " regardless of the user input [default: 0.0]\n"
      */
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
        /*
        else if( strcmp(argv[i],"-m") == 0 )
        {
            misid_ratio = std::stof(argv[i+1]);
            ++i;
        }
        else if( strcmp(argv[i],"-m") == 0 )
        {
            boolean_one = true;
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
                return -2;
            }
            opt.cmndfile = argv[i];
        }
    }

    // Initialize and prepare the ROOT output
    IOManager iomanager(opt.outputFilename);
    iomanager.book_tree_header();
    iomanager.book_tree();
    // Monitor
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
    // process the diagnostic plots
    iomanager.diagnostic_plots();
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
        // --> Not needed, just need the calibration iomanager_cal.book_monitor_plots();
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
        //XXX Fill the diagnostic plots for the calibration
        // iomanager.diagnostic_plots(cal_map);
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
        AlibavaPostProcessor postproc;
        std::cout << " - Calculating <pedestals>" << std::endl;
        PedestalNoiseBeetleMap pednoise_cmmdnot = postproc.calculate_pedestal_noise(iomanager_ped);
        
        std::cout << " - Re-evaluating pedestals (noise extracted)" << std::endl;
        postproc.get_pedestal_noise_free(iomanager_ped,pednoise_cmmdnot);
        
        std::cout << " - Re-calculating pedestals and common noise" << std::endl;
        PedestalNoiseBeetleMap pednoise_cmmd = postproc.calculate_pedestal_noise(iomanager_ped);
        /*for(auto & chip_m: pednoise_cmmdnot)
        {
            std::cout << "CHIP: " << chip_m.first << std::endl;
            std::string up;
            std::string down;
            for(unsigned int i=0; i < chip_m.second.first.size(); ++i)
            {
                up += std::to_string(chip_m.second.first[i])+" ("+
                        std::to_string(pednoise_cmmd[chip_m.first].first[i])+") " ;
                down += std::to_string(chip_m.second.second[i])+" ("+
                        std::to_string(pednoise_cmmd[chip_m.first].second[i])+") ";
            }
            std::cout << up << std::endl;
            std::cout << down << std::endl;
        }*/
        // close the pedestal noise file
        iomanager_ped.close();
        // Fill the diagnostic plots for the pedestal and noise
        iomanager.diagnostic_plots(pednoise_cmmd);
        // And update the beam file with the pedestals and common noise values
        // included in the runHeader postproc 
        iomanager.update(pednoise_cmmd);
    }
    std::cout << std::endl;

    return status;
}


