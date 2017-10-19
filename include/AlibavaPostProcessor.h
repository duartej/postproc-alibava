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

// forward declarations


class AlibavaPostProcessor
{
    public:
        AlibavaPostProcessor();
        AlibavaPostProcessor(const std::map<int,std::vector<int> > & use_channel);
        ~AlibavaPostProcessor() { ; }

        // Check if a channel is masked
        inline bool is_channel_masked(const int & chip, const int & channel) const 
            {  return (_channel_mask.at(chip)[channel] == 0); }
        void print_channels_mask() const;

     
        // Obtain the equivalent electrons per ADC counts
        CalibrateBeetleMap calibrate(IOManager & gauge);

        // Fit a gaussian per channel, the mean is the <pedestal> and the 
        // sigma the <noise> (note that to evaluate the common noise you 
        // need to provide a null signal (or pedestal subtracted) gaussians
        PedestalNoiseBeetleMap calculate_pedestal_noise(const IOManager & pedestal);
        // Get the pedestals with the common noise subtracted (if modify_file_active is 
        // set to true, a new branch will be created using the beam data subtracted
        // the pedestals and noise)
        void get_pedestal_noise_free(const IOManager & pedestal, const PedestalNoiseBeetleMap & mean_ped);

        // Auxiliary functions: static methods useful along the package
        // The mean of a list
        static float get_mean(const std::vector<float> & v);
        // Overloaded method for a map
        static float get_mean(const std::map<int,float> & m);
        // The standard deviation
        static float get_std_dev(const std::vector<float> & v, const float & mean);
        // Overload function for a map
        static float get_std_dev(const std::map<int,float> & v, const float & mean);
        
        // Evaluate the common noise of a set of signals
        static std::pair<float,float> calculate_common_noise(const std::vector<float> & signal, const std::vector<int> & channel_mask);
        // Overloaded version with no channels masked
        static std::pair<float,float> calculate_common_noise(const std::vector<float> & signal);

    private:
        bool _pedestal_subtracted;
        std::string _postproc_treename;
        // The channel mask : 0-masked 1-use it
        std::map<int,std::vector<int> > _channel_mask;
        // Whether or not to automask noisy channels
        bool _automasking;
        // Criteria of automasking
        float _mask_criterium;

        // Function to convert a map mimicking a vector list into a 
        // pure vector
        static std::vector<float> convert_map_in_vector(const std::map<int,float> & m);
        // -- Auto-mask noisy channels
        void mask_channels(const PedestalNoiseBeetleMap & ped_noise);
};

#endif
