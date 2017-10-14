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
#include <algorithm>

StripCluster::StripCluster():
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
    _element_index[_channels.size()-1] = channel;
}

float StripCluster::charge() const
{
    return std::accumulate(_signal_map.begin(),_signal_map.end(),0.0, 
            [] (int previous, const std::map<int,float>::value_type & p) { return previous+p.second; }) ;
}

