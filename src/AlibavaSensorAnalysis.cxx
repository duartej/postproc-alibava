// ***********************************************
// 
// Pool of functions to post-process and analyse 
// the ROOT data created by the fortythieves exec
// (https://github.com/duartej/postproc-alibava/blob/master/bin/fortythieves.cc)
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "AlibavaSensorAnalysis.h"
#include "IOFortythieves.h"

// ROOT headers

// System headers
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

// XXX: Better to have a centralized pool of functions, but if
//      I use the AlibavaPostProcessor versions, I have to link
//      the IOManagersAlibava (to obtain their implementations...)
//      Maybe create a library with useful functions....
//      By the moment,just copy paste the AlibavaPostProcessor funcs.a
namespace auxfunc
{
    std::vector<float> convert_map_in_vector(const std::map<int,float> & m)
    {
        std::vector<float> v;
        v.reserve(m.size());
        for(const auto & _f: m)
        {
            v.push_back(_f.second);
        }
        return v;
    }
    
    float get_mean(const std::vector<float> & v)
    {
        if(v.size() == 0)
        {   
            return 0.0;
        }
        return std::accumulate(v.begin(),v.end(),0.0)/float(v.size());
    }
    
    float get_mean(const std::map<int,float> & m)
    {
        return get_mean(convert_map_in_vector(m));
    }
    
    float get_std_dev(const std::vector<float> & v,const float & mean)
    {
        if(v.size()==0)
        {
            return 0.0;
        }
        // standard deviation = sqrt( E[(x-E[x])^2] )
        std::vector<float> diff(v.size());
        // Get x-E[x] in `diff`
        std::transform(v.begin(),v.end(),diff.begin(),
                [mean](const float & element) { return element-mean; } );
        // obtain the (x-E[x])^2 by using diff*diff (inner_product)
        const float sq_sum = std::inner_product(diff.begin(),diff.end(),diff.begin(),0.0);
        return std::sqrt(sq_sum/static_cast<float>(diff.size()));
    }
    
    float get_std_dev(const std::map<int,float> & m,const float & mean)
    {
        return get_std_dev(convert_map_in_vector(m),mean);
    }
}


AlibavaSensorAnalysis::AlibavaSensorAnalysis(IOFortythieves * io_ft):
    _ft(io_ft),
    _polarity(-1),
    _masked_channels(nullptr),
    _mask_criterium(2.5),
    _snr_seed(5),
    _snr_neighbour(3)
{
    // Fix the size of the TDC cut vector
    _tdc_cut.resize(2);
}

AlibavaSensorAnalysis::~AlibavaSensorAnalysis()
{
    if( _masked_channels != nullptr )
    {
        delete _masked_channels;
        _masked_channels = nullptr;
    }
}

bool AlibavaSensorAnalysis::is_channel_masked(const int & ich)
{
    return (std::find(_masked_channels->begin(),_masked_channels->end(),ich) != _masked_channels->end());
}

void AlibavaSensorAnalysis::mask_channels()
{
    // Asumming configuration XXX: Maybe a state flag stating that?
    if(_masked_channels)
    {
        // The user decide to put this masking 
        // XXX A warning message saying so?
        // Do nothing, the list of channels are already present
        return;
    }
    // Automatically mask algorithm  
    
    // 0 -- Initialize the masked_channels vector
    _masked_channels = new std::vector<int>;

    // 1. -- Define a map to keep only non-noisy channels
    // ichannel: noise
    std::map<int,float> non_noisy_map;
    
    // Get the noise vector and check the noisy channel condition
    // which is given by:
    //    i-strip is noisy if |N_{i} - <Noise> > Xsigma
    for(int ch=0; ch < static_cast<int>(_ft->noise().size()); ++ch)
    {
        // Not using already masked channels
        if(this->is_channel_masked(ch))
        {
            continue;
        }
        non_noisy_map.emplace(ch,_ft->noise()[ch]);
    }

    // 1. Obtain the mean and the standard deviation for the noise
    float mean_noise = auxfunc::get_mean(non_noisy_map);
    float stddev_noise= auxfunc::get_std_dev(non_noisy_map,mean_noise);
    
    // 2. Use the mean and std. dev to find noise channels
    //    while using criteria defined by the user, by
    //    removing iteratively noisy channels until converge, i.e
    //    the non_noisy_map will not loose any element anymore
    unsigned int last_vec_size = _ft->noise().size();
    do
    {
        // after the check from the 'while' statement, update
        // the last know vector size
        last_vec_size = non_noisy_map.size();
        for(auto it = non_noisy_map.cbegin(); it != non_noisy_map.cend(); /* no increment*/)
        {
            if(this->is_channel_masked(it->first))
            {
                ++it;
                continue;
            }
            // Remove the noisy channels and tag it as masked channel
            if( std::abs(it->second-mean_noise)/stddev_noise > _mask_criterium )
            {
                std::cout << "\033[1;34mINFO [AUTO-MASKING]\033[1;m noisy channel: " 
                    << std::setw(3) << it->first << " with noise: " << std::setprecision(3) 
                    << std::setw(5) << it->second << " (<noise>_{all}: " 
                    << std::setprecision(3) << std::setw(5) << mean_noise 
                    << ", st. dev.: " << std::setprecision(3) << std::setw(5) << stddev_noise << ")" 
                    << std::endl;
                _masked_channels->push_back(it->first);
                non_noisy_map.erase(it++);
            }
            else
            {
                // or keep it
                ++it;
            }
        }
        // 3. Recalculate the mean with the excluded signal channels
        mean_noise = auxfunc::get_mean(non_noisy_map);
        stddev_noise= auxfunc::get_std_dev(non_noisy_map,mean_noise);
    } while( non_noisy_map.size() != last_vec_size );
}


StripCluster AlibavaSensorAnalysis::find_clusters(const std::vector<float> & /*adc_corr*/)
{
    return StripCluster();
}


// C-Wrapper to be used by the python module
#ifdef __cplusplus
extern "C"
{
#endif
    // Constructor and destructor
    AlibavaSensorAnalysis * aa_new(IOFortythieves * ioft) { return new AlibavaSensorAnalysis(ioft); }
    void aa_delete(AlibavaSensorAnalysis * aa_inst) { aa_inst->~AlibavaSensorAnalysis(); }
    // Configuration: setters and Getters
    // -- Signal polarity
    void aa_configure_polarity(AlibavaSensorAnalysis * aa_inst, int polarity) { aa_inst->configure_polarity(polarity); }
    int aa_polarity_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_polarity(); }
    // --- Time cut
    void aa_configure_time_cut(AlibavaSensorAnalysis * aa_inst, float t0, float t1) { aa_inst->configure_time_cut(t0,t1); }
    const float * aa_time_cut_getter(AlibavaSensorAnalysis * aa_inst) 
    { 
        // Get the pair and convert it to a float*
        return &(aa_inst->get_time_cut())[0];
    }
    // --- Manually masked channels
    void aa_configure_masked_channels(AlibavaSensorAnalysis * aa_inst, int * mch, int len) 
    { 
        // Convert it in std::vector
        std::vector<int> mch_prov(len);
        mch_prov.assign(mch,mch+len);
        aa_inst->configure_masked_channels(mch_prov); 
    }
    // -- extra function to extrac the size of the masked channesl vector
    int aa_number_masked_channels(AlibavaSensorAnalysis * aa_inst) 
    { 
        if(aa_inst->get_masked_channels() != nullptr)
        {
            return aa_inst->get_masked_channels()->size();
        }
        return 0;
    }
    const int * aa_masked_channels_getter(AlibavaSensorAnalysis * aa_inst)
    {
        // Get the channels vector pointer
        const std::vector<int> * ch_p = aa_inst->get_masked_channels();
        // And check it is not a nullptr
        if(ch_p == nullptr)
        {
            return nullptr;
        }
        // convert the vector to a pointer of ints
        return &(*ch_p)[0];
    }
    void aa_configure_mask_criterium(AlibavaSensorAnalysis * aa_inst, float sigma) { aa_inst->configure_mask_criterium(sigma); }
    float aa_mask_criterium_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_mask_criterium(); }
    // --- The minimum SNR for the seeds of a clusters
    void aa_configure_snr_seed(AlibavaSensorAnalysis * aa_inst, float snr) { aa_inst->configure_snr_seed(snr); }
    float aa_snr_seed_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_snr_seed(); }
    // --- The minimum SNR for a strip candidate of a clusters
    void aa_configure_snr_neighbour(AlibavaSensorAnalysis * aa_inst, float snr) { aa_inst->configure_snr_neighbour(snr); }
    float aa_snr_neighbour_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_snr_neighbour(); }

    // ---------------------------------------------------------
    // Modifier functions
    // ---- Masking channels (part of the initialization of the
    // algorithm)
    void aa_mask_channels(AlibavaSensorAnalysis * aa_inst) { aa_inst->mask_channels(); }
#ifdef __cplusplus
}
#endif



