// ***********************************************
// 
// Cluster definition  (heavly inspired from the 
// AlibavaCluster class at EUTelescope)
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "StripCluster.h"

// ROOT

// System headers

StripCluster::StripCluster() :
    _signal_polarity(0),
    _sensitive_direction(-1),
    _eta_seed(-1),
    _eta_cluster(-1)
{
}
StripCluster::~StripCluster()
{
}

void StripCluster::add(const int & channel, const float & adc)
{
    _signal_map[channel]=adc;
    _channels.push_back(channel); 
}


