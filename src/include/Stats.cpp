/*! \file Stats.cpp
    \brief Implementations of methods in Stats.h.
*/

#include "Stats.h"

void get_lims_mu(double s, double beta, double &mu1, double &mu2) {
    if (s < 0) {
        mu1 = -1.;
        mu2 = (1 - exp(-s)*(1 - beta)) / beta;
    } 
    else {
        mu1 = (1 - exp(-s)*(1 + beta)) / beta;
        mu2 = 1.;
    }
}

double getThomsonScatter(double mu, void *args) {
    struct thom_params *params = (struct thom_params *)args;
    double s = (params->s);
    double beta = (params->beta);

    double mu_prime;

    double output;

    double gamma = beta_gamma(beta);
    double prefac = 3 / (16 * gamma*gamma*gamma*gamma * beta); 

    mu_prime = (exp(s) * (1 - beta*mu) - 1) / beta;
    output = (1 + beta*mu_prime) * (1 + mu*mu * mu_prime*mu_prime +
            0.5 * (1 - mu*mu) * (1 - mu_prime*mu_prime)) / 
            ((1 - beta * mu)*(1 - beta * mu)*(1 - beta * mu)) *
            prefac;

    return output;
}

double getMaxwellJuttner(double beta, double Te) {
    double theta = Te_theta(keV_Temp(Te));
    double gamma;
    double nominator, denominator;

    gamma = beta_gamma(beta);
    nominator = gamma*gamma*gamma*gamma*gamma * beta*beta * exp(-gamma / theta);
    denominator = theta * gsl_sf_bessel_Kn(2, 1/theta);
    return nominator / denominator;

}

double getPowerlaw(double beta, double alpha, double A) {
    double gamma = beta_gamma(beta);
    return A * pow(gamma, -alpha) * beta * pow(1 - beta*beta, -1.5);
}

double getMultiScatteringMJ(double beta, void *args) {
    struct MS_params *ms_params = (struct MS_params *)args;
    double s = (ms_params->s);
    double Te = (ms_params->param);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &getThomsonScatter;

    double pmu, pe, err;
    
    double mu1, mu2;
    get_lims_mu(s, beta, mu1, mu2);
    
    struct thom_params th_params = { s, beta };
    F.params = &th_params;    

    gsl_integration_qag(&F, mu1, mu2, 0, 1e-6, NW_INT, GQMODE, w, &pmu, &err);
    pe = getMaxwellJuttner(beta, Te);
    gsl_integration_workspace_free (w);

    return pmu * pe;
}

double getMultiScatteringPL(double beta, void *args) {
    struct MS_params *ms_params = (struct MS_params *)args;
    double s = (ms_params->s);
    double alpha = (ms_params->param);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (NW_INT);
    gsl_function F;
    F.function = &getThomsonScatter;

    double A;
    double gamma2 = beta_gamma(1 - DBL_EPSILON);
    double gamma1 = 1.;
    getNormPL(gamma1, gamma2, alpha, A);

    double pmu, pe, err;
    
    // Should be calculated from ntSZ
    double mu1, mu2;
    get_lims_mu(s, beta, mu1, mu2);

    struct thom_params th_params = { s, beta };
    F.params = &th_params;    

    gsl_integration_qag(&F, mu1, mu2, 0, 1e-6, NW_INT, GQMODE, w, &pmu, &err);
    pe = getPowerlaw(beta, alpha, A);
    gsl_integration_workspace_free (w);

    return pmu * pe;
}

void getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid) {
    double Da_si = pc_m(Da * 1e6); 
    double theta_c_si = thetac / 3600 / 180 * PI;
    double rc = theta_c_si * Da_si;

    double te0 = ne0*1e6 * ST * rc * sqrt(PI) * gsl_sf_gamma(3/2*ibeta - 0.5) / gsl_sf_gamma(3/2*ibeta);
    
    double theta2;

    if(grid) {
        for(int i=0; i<n_Az; i++) {
            for(int j=0; j<n_El; j++) {
                theta2 = Az[i]*Az[i] + El[j]*El[j];
                output[i*n_El + j] = te0*pow(1 + theta2/(thetac*thetac), 0.5-1.5*ibeta);
            }
        }
    }

    else {
        for(int i=0; i<n_Az; i++) {
            theta2 = Az[i]*Az[i] + El[i]*El[i];
            output[i] = te0*pow(1 + theta2/(thetac*thetac), 0.5-1.5*ibeta);
        }
    }
}

void getNormPL(double &gamma1, double &gamma2, double &alpha, double &A) {
    if(alpha < 0) {
        A = log10(gamma2/gamma1);
        alpha = 1.;
    }
    else {
        A = (1 - alpha) / (pow(gamma2, 1-alpha) - pow(gamma1, 1-alpha));
    }
}
