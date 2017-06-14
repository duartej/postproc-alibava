// ***********************************************
// 
// Pool of functions to post-process alibava data
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef ALIBAVAPOSTPROCESSOR_H
#define ALIBAVAPOSTPROCESSOR_H

#include "IOManager.h"

// System headers
#include <string>
#include <map>
#include <vector>


// improving readability
using PedestalNoiseBeetleMap = std::map<int,std::pair<std::vector<float>,std::vector<float> > >;

class AlibavaPostProcessor
{
    public:
        AlibavaPostProcessor();
        ~AlibavaPostProcessor() { ; }

        // Fit a gaussian per channel, the mean is the <pedestal> and the 
        // sigma the <noise> (note that to evaluate the common noise you 
        // need to provide a null signal (or pedestal subtracted) gaussians
        PedestalNoiseBeetleMap calculate_pedestal_noise(const IOManager & pedestal);
        // Get the pedestals with the common noise subtracted
        void get_pedestal_noise_free(const IOManager & pedestal, const PedestalNoiseBeetleMap & mean_ped);

        // Auxiliary functions: (XXX to be gathered in a pool)
        // The mean of a list
        float get_mean(const std::vector<float> & v);
        // Overloaded method for a map
        float get_mean(const std::map<int,float> & m);
        // The standard deviation
        float get_std_dev(const std::vector<float> & v, const float & mean);
        // Overload function for a map
        float get_std_dev(const std::map<int,float> & v, const float & mean);

    private:
        bool _pedestal_subtracted;
        std::string _postproc_treename;

        // Evaluate the common noise, it must be called previously the subtract_pedestal
        // class
        std::pair<float,float> calculate_common_noise(const std::vector<float> & nullsignal);
        // Function to convert a map mimicking a vector list into a 
        // pure vector
        std::vector<float> convert_map_in_vector(const std::map<int,float> & m);
};

#endif
