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
    std::cout << "1: Invalid number of input arguments" << std::endl;
    std::cout << "2: Undefined input option" << std::endl;
    std::cout << "3: Invalid raw input binary file" << std::endl;
    std::cout << "4: Invalid firmware version from raw data" << std::endl;
    std::cout << "5: Invalid data type (user)" << std::endl;
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
      << " -p alibava raw-binary containing the pedestal runs\n"
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
    int status = IOAlibavaReader::read_data(opt,iomanager);
    iomanager.close();

    // process pedestal file (Calibration ??)
    if(opt.storeHeaderPedestalNoise)
    {
        std::cout << "\033[1;34mfortythieves INFO\033[1;m: " 
            << "Processing pedestal file" << std::endl;
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
        status += IOAlibavaReader::read_data(opt_ped,iomanager_ped);
    
        // calculate pedestals and common noise: { bettle: { channels ,,, } }
        AlibavaPostProcessor postproc;
        std::map<int,std::pair<std::vector<float>,std::vector<float>> > pedestals_cmmdnot = postproc.calculate_pedestal_noise(iomanager_ped);
        //std::map<int,std::vector<float> > noise_cmmd = AlibavaPostProcessor::calculate_commonnoise();
        //std::map<int,std::vector<float> > pedestals_cmmd = AlibavaPostProcessor::calculate_pedestals(pedestals_cmmdnot,noise_cmmd);
        // and update the main beam file
        //iomanager.update(pedestals_cmmd,noise_cmmd);
        // close the pedestal noise file
        iomanager_ped.close();
    }

    // close output ROOT file
    // iomanager.update();
    //iomanager.close();

    return status;
}


