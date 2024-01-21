#include "Signal.h"

void calcSignal_tSZ(double *nu, int n_nu, double Te, double *output, int n_s, int n_beta, bool no_CMB) {
    double s0 = -1.2;
    double s1 = 2.;
    
    double ds = (s1 - s0) / n_s;

    double *s_arr = new double[n_s];
    double *P1_arr = new double[n_s];

    double I_CMB;

    for(int i=0; i<n_s; i++) {
        s_arr[i] = s0 + (i + 0.5)*ds;
    }

    getMultiScatteringMJ(s_arr, n_s, n_beta, Te, P1_arr);    

    for(int i=0; i<n_nu; i++) {
        output[i] = 0.;
        for(int j=0; j<n_s; j++) {
            output[i] += get_CMB(nu[i]*exp(-s_arr[j])) * ds * P1_arr[j];
        }
    }

    if(no_CMB) {
        for(int i=0; i<n_nu; i++) {
            output[i] -= get_CMB(nu[i]);
        }
    }

    delete[] s_arr;
    delete[] P1_arr;
}

void calcSignal_ntSZ(double *nu, int n_nu, double alpha, double *output, int n_s, int n_beta, bool no_CMB) {
    double s0 = -1.2;
    double s1 = 10.;
    
    double ds = (s1 - s0) / n_s;

    double *s_arr = new double[n_s];
    double *P1_arr = new double[n_s];

    double I_CMB;

    for(int i=0; i<n_s; i++) {
        s_arr[i] = s0 + (i + 0.5)*ds;
    }

    getMultiScatteringPL(s_arr, n_s, n_beta, alpha, P1_arr);    

    for(int i=0; i<n_nu; i++) {
        output[i] = 0.;
        for(int j=0; j<n_s; j++) {
            output[i] += get_CMB(nu[i]*exp(-s_arr[j])) * ds * P1_arr[j];
        }
    }

    if(no_CMB) {
        for(int i=0; i<n_nu; i++) {
            output[i] -= get_CMB(nu[i]);
        }
    }

    delete[] s_arr;
    delete[] P1_arr;
}


void calcSignal_kSZ(double *nu, int n_nu, double v_pec, double *output, int n_mu) {
    double beta_pec = v_beta(v_pec*1e3);
    double gamma_pec = v_gamma(v_pec);

    double mu0 = -1.;
    double mu1 = 1.;
    double mu;
    double I_CMB;

    double dmu = (mu1 - mu0) / n_mu;

    double x1, x2;

    for(int i=0; i<n_nu; i++) {
        x1 = CH * nu[i] / KB / TCMB;
        I_CMB = get_CMB(nu[i]);
        output[i] = 0.;

        for(int j=0; j<n_mu; j++) {
            mu = mu0 + (j + 0.5) * dmu;
            x2 = x1 * gamma_pec*gamma_pec * (1 + beta_pec) * (1 - beta_pec * mu);
            output[i] += I_CMB * 3/8 * (1 + mu*mu) * ((exp(x1) - 1) / (exp(x2) - 1) - 1) * dmu;
        }
    }
}
