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
    MOCKSZ_DLL void MockSZ_getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu);

    MOCKSZ_DLL void MockSZ_getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output); 
    
    MOCKSZ_DLL void MockSZ_getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output);
    
    MOCKSZ_DLL void MockSZ_getMultiScatteringMJ(double *s_arr, int n_s, double Te, double *output, int n_beta);
    
    MOCKSZ_DLL void MockSZ_getMultiScatteringPL(double *s_arr, int n_s, double alpha, double *output, int n_beta);
    
    MOCKSZ_DLL void MockSZ_getSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);
    
    MOCKSZ_DLL void MockSZ_getSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);

    MOCKSZ_DLL void MockSZ_getSignal_kSZ(double *nu, int n_nu, double v_pec, double tau_e, double *output, int n_mu);

    MOCKSZ_DLL void MockSZ_getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid);

    MOCKSZ_DLL void MockSZ_getCMB(double *nu, int n_nu, double *output);
}

#endif
