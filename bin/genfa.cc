// ***************************************************************************
// genfa: Get the Event Number From an Alibava file
//
// Obtain the number of events available at a alibava binary raw data 
//
// See the section 'The Alibava data format' at 
// https://www.alibavasystems.com/images/Catalogo/alibava-usermanual.pdf
// 
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ****************************************************************************

// system header
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

// Alibava reader
#include "IOAlibavaReader.h"

// The input options
#include "input_options.h"


// Input parser functions
void display_usage()
{
    std::cout << "\033[37musage:\033[m genfaol [OPTIONS] alibava_data.raw"
        << std::endl;
    std::cout << std::endl;
    std::cout << "Extract the number of events from ALIBAVA raw binary data" 
        << std::endl;
    std::cout << std::endl;
    std::cout << "[OPTIONS]\n"
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

    // Decide if the file is alibava or LCIO
    
    // process the filea
    const int events = IOAlibavaReader::read_data(opt,nullptr);
    if(events < 0)
    {
        return events;
    }
    std::cout << events << std::endl;
    return 0;
}
