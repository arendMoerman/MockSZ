/*! \file InterfaceCPU.cpp
    \brief Implementations of ctypes bindings.
*/

#include "InterfaceCPU.h"

MOCKSZ_DLL void MockSZ_getThomsonScatter(double *s_arr, int n_s, double beta, double *output, double acc) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &getThomsonScatter;

    
    double mu1, mu2;
    double err;
    for(int i=0; i<n_s; i++) {
        get_lims_mu(s_arr[i], beta, mu1, mu2);
        
        struct thom_params params = { s_arr[i], beta };

        F.params = &params;    
        gsl_integration_qag(&F, mu1, mu2, acc, acc, NW_INT, GQMODE, w, &(output[i]), &err);
        if(output[i] < 0) {output[i] = 0;}
    }
    gsl_integration_workspace_free (w);
}

MOCKSZ_DLL void MockSZ_getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output, double acc) {
    for(int i=0; i<n_beta; i++) {
        output[i] = getMaxwellJuttner(beta_arr[i], Te);
    }
} 

MOCKSZ_DLL void MockSZ_getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output, double acc) {
    double A;
    double gamma2 = beta_gamma(1 - DBL_EPSILON);
    double gamma1 = 1.;
    
    getNormPL(gamma1, gamma2, alpha, A);
    
    for(int i=0; i<n_beta; i++) {
        output[i] = getPowerlaw(beta_arr[i], alpha, A);
    }
}

MOCKSZ_DLL void MockSZ_getMultiScatteringMJ(double *s_arr, int n_s, double Te, double *output, double acc) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &getMultiScatteringMJ;
    double beta0, err;
    
    for(int i=0; i<n_s; i++) {
        beta0 = (exp(abs(s_arr[i])) - 1) / (exp(abs(s_arr[i])) + 1) + DBL_EPSILON;
        
        struct MS_params ms_params = { s_arr[i], Te };

        F.params = &ms_params;    
        gsl_integration_qag(&F, beta0, BETA1, acc, acc, NW_INT, GQMODE, w, &(output[i]), &err);
    }
    gsl_integration_workspace_free (w);
}

MOCKSZ_DLL void MockSZ_getMultiScatteringPL(double *s_arr, int n_s, double alpha, double *output, double acc) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &getMultiScatteringPL;
    double beta0, err;
    
    for(int i=0; i<n_s; i++) {
        beta0 = (exp(abs(s_arr[i])) - 1) / (exp(abs(s_arr[i])) + 1) + DBL_EPSILON;
        
        struct MS_params ms_params = { s_arr[i], alpha };

        F.params = &ms_params;    
        gsl_integration_qag(&F, beta0, BETA1, acc, acc, NW_INT, GQMODE, w, &(output[i]), &err);
    }
    gsl_integration_workspace_free (w);
}

MOCKSZ_DLL void MockSZ_getSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, double acc) { 
    double s0 = -3;
    double s1 = 3;
    
    double *func_evals = new double[MEVALS * MEVALS];
    int n_eval = romberg_write(&get_n_eval, &getMultiScatteringMJ, s0, s1, Te, func_evals, MEVALS, acc); 
    for(int i=0; i<n_nu; i++) {
        output[i] = 0.;
        double args[2] = {nu[i], tau_e};
        
        output[i] += romberg_read(&conv_CMB_scatt, s0, s1, args, func_evals, n_eval);
        output[i] -= tau_e * get_CMB(nu[i]);
    }
    delete[] func_evals;
}

MOCKSZ_DLL void MockSZ_getSignal_tSZ_beta2(double *nu, int n_nu, double Te, double tau_e, double *output, double beta) {
    for(int i=0; i<n_nu; i++) {
        output[i] = tau_e * calcSignal_tSZ_beta2(nu[i], Te, beta);
    }
}

MOCKSZ_DLL void MockSZ_getSignal_kSZ_betatheta(double *nu, int n_nu, double Te, double tau_e, double *output, double prefac) {
    for(int i=0; i<n_nu; i++) {
        output[i] = tau_e * calcSignal_kSZ_betatheta(nu[i], Te, prefac);
    }
}

MOCKSZ_DLL void MockSZ_getSignal_kSZ_betat2heta(double *nu, int n_nu, double Te, double tau_e, double *output, double prefac) {
    for(int i=0; i<n_nu; i++) {
        output[i] = tau_e * calcSignal_kSZ_betat2heta(nu[i], Te, prefac);
    }
}

MOCKSZ_DLL void MockSZ_getSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, double acc) {
    double s0 = -9;
    double s1 = 18;
    
    double *func_evals = new double[MEVALS * MEVALS];
    int n_eval = romberg_write(&get_n_eval, &getMultiScatteringPL, s0, s1, alpha, func_evals, MEVALS, acc); 

    for(int i=0; i<n_nu; i++) {
        output[i] = 0.;
        double args[2] = {nu[i], tau_e};
        
        output[i] += romberg_read(&conv_CMB_scatt, s0, s1, args, func_evals, n_eval);
        output[i] -= tau_e * get_CMB(nu[i]);
    }
    delete[] func_evals;
}

MOCKSZ_DLL void MockSZ_getSignal_kSZ(double *nu, int n_nu, double v_pec, double tau_e, double *output, double acc) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &calcSignal_kSZ;
    
    double mu0 = -1.;
    double mu1 = 1.;
    double err;
    
    for(int i=0; i<n_nu; i++) {
        struct kSZ_params ksz_params = { nu[i], v_pec, tau_e };
        F.params = &ksz_params;    
        
        gsl_integration_qag(&F, mu0, mu1, acc, acc, NW_INT, GQMODE, w, &(output[i]), &err);
    }
    gsl_integration_workspace_free (w);
}

MOCKSZ_DLL void MockSZ_getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid) {
    getIsoBeta(Az, El, n_Az, n_El, ibeta, ne0, thetac, Da, output, grid);
}
    
MOCKSZ_DLL void MockSZ_getCMB(double *nu, int n_nu, double *output) {
    for(int i=0; i<n_nu; i++) {
        output[i] = get_CMB(nu[i]);
    }
}
