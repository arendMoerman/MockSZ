/*! \file Signal.cpp
    \brief Implementations of the single-pointing SZ signals in Signal.h.
*/

#include "Signal.h"

double calcSignal_kSZ(double mu, void *args) {
    struct kSZ_params *ksz_params = (struct kSZ_params *)args;
    double nu = (ksz_params->nu);
    double tau_e = (ksz_params->tau_e);
    
    double beta_pec = (ksz_params->v_pec);
    double gamma_pec = beta_gamma(beta_pec);

    double I_CMB;

    double x1, x2;

    x1 = CH * nu / KB / TCMB;
    I_CMB = get_CMB(nu);

    x2 = x1 * gamma_pec*gamma_pec * (1 + beta_pec) * (1 - beta_pec * mu);
    return tau_e * I_CMB * 3./8. * (1 + mu*mu) * ((exp(x1) - 1) / (exp(x2) - 1) - 1);
}

double calcSignal_tSZ_beta2(double nu, double Te, double beta) {
    double X = nu_x(nu);
    double theta = Te_theta(keV_Temp(Te));
    double eX = exp(X);
    double Xt = X * cosh(X/2) / sinh(X/2);
    double St = X / sinh(X/2);

    double Y0 = Xt - 4;
    double Y1 = -10 + Xt*(47./2 - Xt*(42./5 - Xt*7./10)) + St*(-21./5 + Xt*7./5);

    double conv = 2*CH * nu*nu*nu / CL/ CL * eX / (eX-1) / (eX-1);

    return conv * beta*beta * (Y0/3 + theta*(5.*Y0/6 + 2.*Y1/3));
}

double calcSignal_kSZ_betatheta(double nu, double Te, double prefac) {
    double X = nu_x(nu);
    double theta = Te_theta(keV_Temp(Te));
    double eX = exp(X);
    double Xt = X * cosh(X/2) / sinh(X/2);
    double St = X / sinh(X/2);

    double C1 = 10 - Xt*(47./5 - Xt*7./5) + 7./10*St*St;
    double C2 = 25 + Xt*(-111.7 + Xt*(84.7 + Xt*(-18.3 + 11./10*Xt))) + 
            St*St*(84.7/2 + Xt*(-183./5 + 12.1/2*Xt) + 11./10*St*St);

    double conv = 2*CH * nu*nu*nu / CL/ CL * X * eX / (eX-1) / (eX-1);

    return -conv * prefac * theta*(C1 + theta*C2);
}

double calcSignal_kSZ_betat2heta(double nu, double Te, double prefac) {
    double X = nu_x(nu);
    double theta = Te_theta(keV_Temp(Te));
    double eX = exp(X);
    double Xt = X * cosh(X/2) / sinh(X/2);
    double St = X / sinh(X/2);
    
    double D0 = -2./3 + 11./30*Xt;
    double D1 = -4 + Xt*(12 + Xt*(-6 + 19./30*Xt)) + St*St*(-3 + 19./15 * Xt);

    double conv = 2*CH * nu*nu*nu / CL/ CL * X * eX / (eX-1) / (eX-1);

    return -conv * prefac * (D0 + theta * D1);
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
