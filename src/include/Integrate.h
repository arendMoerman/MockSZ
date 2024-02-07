/*! \file Integrate.h
    \brief Declarations of the Romberg integrator for MockSZ.
*/

#include "Constants.h"
#include "Conversions.h"
#include "Utils.h"
#include "Stats.h"

#ifndef __Integrate_h
#define __Integrate_h

template <class T>
class Romberg {
    int k;
    int j;

    T &func;

    Romberg(T &func);
    

} 

#endif

