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
      << " -p flag to store pedestal and noise header [false]\n"
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
    int status = iomanager.read_data(opt);

    // close output ROOT file
    iomanager.close();

    return status;
}


