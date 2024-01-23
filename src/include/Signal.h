#include "Constants.h"
#include "Conversions.h"
#include "Utils.h"
#include "Stats.h"

#ifndef __Signal_h
#define __Signal_h

void calcSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);

void calcSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);

void calcSignal_kSZ(double *nu, int n_nu, double v_pec, double tau_e, double *output, int n_mu);

#endif
