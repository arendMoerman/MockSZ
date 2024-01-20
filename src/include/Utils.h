#include "math.h"
#include "Constants.h"

#ifndef __Utils_h
#define __Utils_h

double get_min(double *arr, int n_arr);
double get_CMB(double nu);


inline double get_min(double *arr, int n_arr) {
    double min = arr[0];
    for(int i=1; i<n_arr; i++) {
        if(arr[i] < min) {min = arr[i];}
    }

    return min;
}

inline double get_CMB(double nu) {
    double prefac = 2 * CH * nu*nu*nu / (CL*CL);
    double distri = 1 / (exp(CH * nu / (KB * TCMB)) - 1);

    return prefac*distri;
}
#endif
