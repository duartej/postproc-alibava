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
#include <map>

StripCluster::StripCluster():
    _signal_polarity(0),
    _sensitive_direction(-1),
    _eta_calculated(false),
    _eta_seed(-1),
    _eta_cluster(-1)
{
}
StripCluster::~StripCluster()
{
}

void StripCluster::add(const int & channel, const float & adc)
{
    // Re-evaluate the eta when ask for it
    _eta_calculated = false;
    // add the strip
    _signal_map[channel]=adc;
    _channels.push_back(channel); 
    _element_index[_channels.size()-1] = channel;
}

int StripCluster::left_channel() const
{
    return *std::min_element(_channels.begin(),_channels.end());
}

int StripCluster::right_channel() const
{
    return *std::max_element(_channels.begin(),_channels.end());
}

float StripCluster::charge() const
{
    // -- Memorize it as well?
    return std::accumulate(_signal_map.begin(),_signal_map.end(),0.0, 
            [this] (const float & previous, const std::map<int,float>::value_type & p) { return previous+_signal_polarity*p.second; }) ;
}

float StripCluster::calibrated_charge(const std::vector<float> & calibration) 
{
    float total_charge = 0.0;
    for(int ch=0; ch < static_cast<int>(calibration.size()); ++ch)
    {
        if(std::find(_channels.begin(),_channels.end(),ch) != _channels.end())
        {
            total_charge += _signal_polarity*_signal_map[ch]*calibration[ch];
        }
        
    }
    return total_charge;
}

float StripCluster::snr(const std::vector<float> & noise)
{
    float cluster_noise = 0.0;
    for(int ch=0; ch < static_cast<int>(noise.size()); ++ch)
    {
        if(std::find(_channels.begin(),_channels.end(),ch) != _channels.end())
        {
            cluster_noise += noise[ch];
        }
        
    }
    return this->charge()/cluster_noise;
}

float StripCluster::eta()
{
    // First check if the data was already calculated, not
    // doing it again
    if(_eta_calculated)
    {
        return _eta_cluster;
    }

    // [1] It would define the position by subtracting the
    // right strip  x = X_R - p f(eta)
    // In case of clusters > 2, it will use the two higher
    // signal channels
    
    // Calculate the charge sharing for clusters > 2, otherwise
    // using the eta_seed
    if(this->size() == 1)
    {
        _eta_cluster = _eta_seed;
        _eta_calculated = true;
        return _eta_cluster;
    }

    // Note that the seed is always the first element of the cluster
    std::multimap<float,int> channels_ordered;
    // Store the multimap to use the highest signal channels
    for(const auto & channel_signal: _signal_map)
    {
        channels_ordered.insert(std::pair<float,int>(_signal_polarity*channel_signal.second,channel_signal.first));
    }
    // Use only the two last elements (the highest signals) of the map
    auto highest_signal = channels_ordered.rbegin();
    auto next_to_highest_signal = std::next(channels_ordered.rbegin());
    
    // A guess
    auto left_signal = highest_signal;
    auto right_signal = next_to_highest_signal;
    // Check if the guess was right, i.e. the highest signal channel
    // is the leftest channel (i.e, the lower channel), otherwise, change the order
    if( left_signal->second > right_signal->second)
    {
        left_signal = next_to_highest_signal;
        right_signal= highest_signal;
    }
    // Update the data-member and memorize 
    //_eta_cluster = left_signal->first/(left_signal->first+right_signal->first); --> DEPRECATED
    _eta_cluster = right_signal->first/(left_signal->first+right_signal->first);
    _eta_calculated = true;

    return _eta_cluster;
}

float StripCluster::position_cog()
{
    if(this->size() == 1)
    {
        // The unique member, therefore its position  
        return _channels[0];
    }

    float total_signal = 0;
    float weighted_sum = 0;
    for(int i=0; i < static_cast<int>(this->size()) ; ++i) 
    {
        total_signal += this->charge(i);
        weighted_sum += this->charge(i)*_channels[i];
    }
	
    return weighted_sum/total_signal;
}

