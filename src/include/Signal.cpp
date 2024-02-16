/*! \file Signal.cpp
    \brief Implementations of the single-pointing SZ signals in Signal.h.
*/

#include "Signal.h"

double calcSignal_kSZ(double mu, void *args) {
    struct kSZ_params *ksz_params = (struct kSZ_params *)args;
    double nu = (ksz_params->nu);
    double v_pec = (ksz_params->v_pec);
    double tau_e = (ksz_params->tau_e);
    
    double beta_pec = v_beta(v_pec*1e3);
    double gamma_pec = v_gamma(v_pec);

    double I_CMB;

    double x1, x2;

    x1 = CH * nu / KB / TCMB;
    I_CMB = get_CMB(nu);

    x2 = x1 * gamma_pec*gamma_pec * (1 + beta_pec) * (1 - beta_pec * mu);
    return tau_e * I_CMB * 3./8. * (1 + mu*mu) * ((exp(x1) - 1) / (exp(x2) - 1) - 1);
}

double get_CMB(double nu) {
    double prefac = 2 * CH * nu*nu*nu / (CL*CL);
    double distri = 1 / (exp(CH * nu / (KB * TCMB)) - 1);

    return prefac*distri;
}

double get_n_eval(double (*func)(double, void*), double s, double arg) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = func;

    double beta0, err, output;
    
    beta0 = (exp(abs(s)) - 1) / (exp(abs(s)) + 1) + DBL_EPSILON;
    
    struct MS_params ms_params = { s, arg };

    F.params = &ms_params;    
    gsl_integration_qag(&F, beta0, BETA1, 1e-6, 1e-6, NW_INT, GQMODE, w, &output, &err);
    gsl_integration_workspace_free (w);

    return output;
}

double conv_CMB_scatt(double s, double *args) {
    double nu = args[0];
    double tau_e = args[1];
    return tau_e * get_CMB(nu*exp(-s));
}
