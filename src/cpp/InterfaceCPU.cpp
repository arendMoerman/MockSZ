/*! \file InterfaceCPU.cpp
    \brief Implementations of ctypes bindings.
*/

#include "InterfaceCPU.h"

MOCKSZ_DLL void MockSZ_getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu) {
    getThomsonScatter(s_arr, n_s, beta, output, num_mu);
}

MOCKSZ_DLL void MockSZ_getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output) {
    getMaxwellJuttner(beta_arr, n_beta, Te, output);
} 

MOCKSZ_DLL void MockSZ_getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output) {
    getPowerlaw(beta_arr, n_beta, alpha, output);
}

MOCKSZ_DLL void MockSZ_getMultiScatteringMJ(double *s_arr, int n_s, double Te, double *output, int n_beta, int nThreads) {
    getMultiScatteringMJ(s_arr, n_s, n_beta, Te, output, nThreads);
}

MOCKSZ_DLL void MockSZ_getMultiScatteringPL(double *s_arr, int n_s, double alpha, double *output, int n_beta, int nThreads) {
    getMultiScatteringPL(s_arr, n_s, n_beta, alpha, output, nThreads);
}

MOCKSZ_DLL void MockSZ_getSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, int n_s, int n_beta, bool no_CMB, int nThreads) {
    calcSignal_tSZ(nu, n_nu, Te, tau_e, output, n_s, n_beta, no_CMB, nThreads);
}

MOCKSZ_DLL void MockSZ_getSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, int n_s, int n_beta, bool no_CMB, int nThreads) {
    calcSignal_ntSZ(nu, n_nu, alpha, tau_e, output, n_s, n_beta, no_CMB, nThreads);
}

MOCKSZ_DLL void MockSZ_getSignal_kSZ(double *nu, int n_nu, double v_pec, double tau_e, double *output, int n_mu) {
    calcSignal_kSZ(nu, n_nu, v_pec, tau_e, output, n_mu);
}

MOCKSZ_DLL void MockSZ_getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid) {
    getIsoBeta(Az, El, n_Az, n_El, ibeta, ne0, thetac, Da, output, grid);
}
    
MOCKSZ_DLL void MockSZ_getCMB(double *nu, int n_nu, double *output) {
    for(int i=0; i<n_nu; i++) {
        output[i] = get_CMB(nu[i]);
    }
}
