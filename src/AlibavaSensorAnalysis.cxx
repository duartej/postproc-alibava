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
//

AlibavaSensorAnalysis::AlibavaSensorAnalysis(IOFortythieves * io_ft):
    _ft(io_ft),
    _polarity(-1),
    _masked_channels(nullptr),
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
    // --- The minimum SNR for the seeds of a clusters
    void aa_configure_snr_seed(AlibavaSensorAnalysis * aa_inst, float snr) { aa_inst->configure_snr_seed(snr); }
    float aa_snr_seed_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_snr_seed(); }
    // --- The minimum SNR for a strip candidate of a clusters
    void aa_configure_snr_neighbour(AlibavaSensorAnalysis * aa_inst, float snr) { aa_inst->configure_snr_neighbour(snr); }
    float aa_snr_neighbour_getter(AlibavaSensorAnalysis * aa_inst) { return aa_inst->get_snr_neighbour(); }
#ifdef __cplusplus
}
#endif



