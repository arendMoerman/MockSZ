/*! \file Romberg.cpp
   \brief Implementations of functions used for romberg integrations.
*/

#include <stdio.h>
#include <math.h>

#include "Romberg.h"

int romberg_write(double (*f)(double (*g)(double, void*), double, double), double (*gg)(double, void*), double a, double b, double arg, double *write_arr, size_t max_steps, double acc) {
    double R1[max_steps], R2[max_steps]; // buffers
    double *Rp = &R1[0], *Rc = &R2[0]; // Rp is previous row, Rc is current row
    double h = b-a; //step size
    write_arr[0] = f(gg, a, arg);
    write_arr[1] = f(gg, b, arg);

    Rp[0] = (write_arr[0] + write_arr[1])*h*0.5; // first trapezoidal step
  
    int n_eval = 2;

    for (size_t i = 1; i < max_steps; ++i) {
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1); //2^(n-1)
        for (size_t j = 1; j <= ep; ++j) {
            write_arr[n_eval] = f(gg, a + (2*j-1) * h, arg);
            c += write_arr[n_eval];

            n_eval++;
        }
        Rc[0] = h*c + .5*Rp[0]; // R(i,0)

        for (size_t j = 1; j <= i; ++j) {
            double n_k = pow(4, j);
            Rc[j] = (n_k*Rc[j-1] - Rp[j-1]) / (n_k-1); // compute R(i,j)
        }

        if (i > 1 && fabs(Rp[i-1]-Rc[i]) < acc) {
            //printf("%.3e %d \n", write_arr[n_eval-1], n_eval);
            return i;
        }

        // swap Rn and Rc as we only need the last row
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
    }
    return max_steps; // return our best guess
}

double romberg_read(double (*f)(double, double*), double a, double b, double *args, double *read_arr, size_t n_eval) {
    double R1[n_eval], R2[n_eval]; // buffers
    double *Rp = &R1[0], *Rc = &R2[0]; // Rp is previous row, Rc is current row
    double h = b-a; //step size

    Rp[0] = (read_arr[0]*f(a, args) + read_arr[1]*f(b, args))*h*0.5; // first trapezoidal step
  
    int n = 2;

    for (size_t i = 1; i <= n_eval; ++i) {
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1); //2^(n-1)
        for (size_t j = 1; j <= ep; ++j) {
            c += read_arr[n] * f(a + (2*j-1) * h, args);

            n++;
        }
        Rc[0] = h*c + .5*Rp[0]; // R(i,0)

        for (size_t j = 1; j <= i; ++j) {
            double n_k = pow(4, j);
            Rc[j] = (n_k*Rc[j-1] - Rp[j-1]) / (n_k-1); // compute R(i,j)
        }
        // swap Rn and Rc as we only need the last row
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
    }
    return Rc[n_eval-1];
}
