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

void getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output) {
    double A;
    double gamma2 = beta_gamma(1 - DBL_EPSILON);
    double gamma1 = beta_gamma(beta_arr[0]);
    double output_scalar;

    if(alpha < 0) {
        A = log10(gamma2/gamma1);
        alpha = 1.;
    }
    else {
        A = (1 - alpha) / (pow(gamma2, 1-alpha) - pow(gamma1, 1-alpha));
    }

    double beta;
    for(int i=0; i<n_beta; i++) {
        beta = beta_arr[i];
        getPowerlaw(beta, alpha, A, output_scalar);
        output[i] = output_scalar;
    }
}

void getPowerlaw(double beta, double alpha, double A, double &output) {
    double gamma = beta_gamma(beta);
    output = A * pow(gamma, -alpha) * beta * pow(1 - beta*beta, -1.5);
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

void getMultiScatteringPL(double *s_arr, int n_s, int n_beta, double alpha, double *output) {
    double beta0;
    double dbeta;
    double beta;
    
    double A;
    double gamma2 = beta_gamma(1 - DBL_EPSILON);

    double *beta0_arr = new double[n_s];
    for(int i=0; i<n_s; i++) {
        beta0_arr[i] = (exp(abs(s_arr[i])) - 1) / (exp(abs(s_arr[i])) + 1);
    }
    
    double min_beta0 = get_min(beta0_arr, n_s);
    double gamma1 = beta_gamma(min_beta0);

    double output_scalar1, output_scalar2;

    if(alpha < 0) {
        A = log10(gamma2/gamma1);
        alpha = 1.;
    }
    else {
        A = (1 - alpha) / (pow(gamma2, 1-alpha) - pow(gamma1, 1-alpha));
    }

    for(int i=0; i<n_s; i++) {
        dbeta = (1 - DBL_EPSILON - beta0_arr[i]) / n_beta;
        output[i] = 0.;

        for(int j=0; j<n_beta; j++) {
            beta = beta0_arr[i] + (j + 0.5) * dbeta;
            getThomsonScatter(s_arr[i], beta, output_scalar1);
            getPowerlaw(beta, alpha, A, output_scalar2);

            output[i] += output_scalar1 * output_scalar2 * dbeta;
        }
    }
    delete[] beta0_arr;
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
