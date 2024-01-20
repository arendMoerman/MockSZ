#include "Stats.h"
#include "Signal.h"

#ifdef _WIN32
#   define MOCKSZ_DLL __declspec(dllexport)
#else
#   define MOCKSZ_DLL
#endif

#ifndef __InterfaceCPU_h
#define __InterfaceCPU_h

extern "C"
{
    MOCKSZ_DLL void MockSZ_getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu=500);

    MOCKSZ_DLL void MockSZ_getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output); 
    
    MOCKSZ_DLL void MockSZ_getMultiScatteringMJ(double *s_arr, int n_s, double Te, double *output, int n_beta=500);
    
    MOCKSZ_DLL void MockSZ_getSignal_tSZ(double *nu, int n_nu, double Te, double *output, int n_s=500, int n_beta=500, bool no_CMB=false);

    MOCKSZ_DLL void MockSZ_getSignal_kSZ(double *nu, int n_nu, double v_pec, double *output, int n_mu=500);
}

#endif
