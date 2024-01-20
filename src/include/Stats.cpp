#include "Stats.h"

void getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu) {
    double output_scalar;

    for(int i=0; i<n_s; i++) {
        getThomsonScatter(s_arr[i], beta, output_scalar, num_mu);
        output[i] = output_scalar;
    }
}

void getThomsonScatter(double s, double beta, double &output, int num_mu) {
    double mu1;
    double mu2;
    double mu;
    double mu_prime;
    double dmu;

    double integrand = 0;

    double gamma = beta_gamma(beta);
    double prefac = 3 / (16 * gamma*gamma*gamma*gamma * beta); 

    if (s < 0) {
        mu1 = -1.;
        mu2 = (1 - exp(-s)*(1 - beta)) / beta;
    } 

    else {
        mu1 = (1 - exp(-s)*(1 + beta)) / beta;
        mu2 = 1.;
    }
    
    dmu = (mu2 - mu1) / num_mu;

    for(int j=0; j<num_mu; j++) {
        mu = mu1 + (j + 0.5) * dmu;
        mu_prime = (exp(s) * (1 - beta*mu) - 1) / beta;

        integrand += (1 + beta*mu_prime) * (1 + mu*mu * mu_prime*mu_prime +
                0.5 * (1 - mu*mu) * (1 - mu_prime*mu_prime)) *
                dmu / ((1 - beta * mu)*(1 - beta * mu)*(1 - beta * mu));
    }
    if (integrand < 0) {
        output = 0.;
    }
    else {
        output = prefac * integrand;
    }
    integrand = 0.;
}

void getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output) {
    double beta;
    double output_scalar;

    for(int i=0; i<n_beta; i++) {
        beta = beta_arr[i];
        getMaxwellJuttner(beta, Te, output_scalar);
        output[i] = output_scalar;
    }
}

void getMaxwellJuttner(double beta, double Te, double &output) {
    double theta = Te_theta(KeV_Temp(Te));
    double gamma;
    double nominator, denominator;

    gamma = beta_gamma(beta);
    nominator = gamma*gamma*gamma*gamma*gamma * beta*beta * exp(-gamma / theta);
    denominator = theta * gsl_sf_bessel_Kn(2, 1/theta);
    output = nominator / denominator;
}

void getMultiScatteringMJ(double *s_arr, int n_s, int n_beta, double Te, double *output) {
    double beta0;
    double dbeta;
    double beta;
    double beta1 = 1 - DBL_EPSILON;

    double output_scalar1, output_scalar2;
    for(int i=0; i<n_s; i++) {
        beta0 = (exp(abs(s_arr[i])) - 1) / (exp(abs(s_arr[i])) + 1); 
        dbeta = (beta1 - beta0) / n_beta;
        
        output[i] = 0.;
        for(int j=0; j<n_beta; j++) {
            beta = beta0 + (j + 0.5) * dbeta;
            getThomsonScatter(s_arr[i], beta, output_scalar1);
            getMaxwellJuttner(beta, Te, output_scalar2);

            output[i] += output_scalar1 * output_scalar2 * dbeta;
        }
    }
}

