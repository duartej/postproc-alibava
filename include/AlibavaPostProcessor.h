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
#include <map>
#include <vector>


// improving readability
using PedestalNoiseBeetleMap = std::map<int,std::pair<std::vector<float>,std::vector<float> > >;

class AlibavaPostProcessor
{
    public:
        AlibavaPostProcessor() { ; }
        ~AlibavaPostProcessor() { ; }

        // Fit a gaussian per channel, the mean is the <pedestal> and the 
        // sigma the <noise> (note that to evaluate the common noise you 
        // need to provide a null signal (or pedestal subtracted) gaussians
        PedestalNoiseBeetleMap calculate_pedestal_noise(const IOManager & pedestal);
    //private:

};

#endif
