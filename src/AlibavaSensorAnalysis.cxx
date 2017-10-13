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
    _tdc_t0(0.0),
    _tdc_t1(100.0),
    _masked_channels(nullptr),
    _snr_seed(0),
    _snr_neighbour(0)
{
}

StripCluster AlibavaSensorAnalysis::find_clusters(const std::vector<float> & /*adc_corr*/)
{
    return StripCluster();
}

#include <iostream>

// C-Wrapper to be used by the python extension
#ifdef __cplusplus
extern "C"
{
#endif
    // Constructor and destructor
    AlibavaSensorAnalysis * aa_new(IOFortythieves * ioft) { return new AlibavaSensorAnalysis(ioft); }
    void aa_delete(AlibavaSensorAnalysis * aa_inst) { aa_inst->~AlibavaSensorAnalysis(); }
    // Useful functions
    void aa_configure_time_cut(AlibavaSensorAnalysis * aa_inst, float t0, float t1) { aa_inst->configure_time_cut(t0,t1); }
    void aa_configure_masked_channels(AlibavaSensorAnalysis * aa_inst, int * mch, int len) 
    { 
        // Convert it in std::vector
        std::vector<int> mch_prov(len);
        mch_prov.assign(mch,mch+len);
        aa_inst->configure_masked_channels(mch_prov); 
    }
#ifdef __cplusplus
}
#endif



