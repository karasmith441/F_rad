#ifndef WINDPARAMS_HPP
#define WINDPARAMS_HPP

#include <iostream>
#include <math.h>
#include <vector>

#include "consts.hpp"
#include "integrator.hpp"

#define ldouble long double

/*==============================================================================
This class handles storage of wind parameters and calculates a sight line
dependent opacity as well as the directional gradient of the velocity. Like
the Disk class, it requires it's observer to be set, but it takes both the
spherical polar and cartesian coordinates as well as the trig functions of
theta. This class should really only be used in the context of the Disk.
==============================================================================*/
class WindParams
{
public:

  // parameters to describe the wind
  ldouble M;                           // Central object mass
  ldouble k;                           // CAK k parameter
  ldouble alpha;                       // CAK alpha parameter
  ldouble v_th;                        // Thermal velocity
  ldouble rho;                         // density
  ldouble N_0;                         // Line distribution normalization factor
  ldouble sigma_ref;                   // reference opacity
  ldouble eta_max;                     // maximum opacity
  ldouble sz;                          // power law for wind velocity field

  // positional parameters
  ldouble r0, theta0, x0, z0;
  ldouble ct, st, tt;

  WindParams() {};
  WindParams(ldouble M, ldouble k, ldouble alpha, ldouble v_th, ldouble rho, ldouble N_0, ldouble sigma_ref, ldouble eta_max, ldouble sz);

  // Sets the position of the observer in the wind
  void setObserver(ldouble r0, ldouble theta0, ldouble x0, ldouble z0, ldouble ct, ldouble st, ldouble tt);

  ldouble get_t(ldouble mu, ldouble phi_1);        // returns the dimensionless optical depth parameter, t
  ldouble dv_vertical(ldouble mu, ldouble phi_1);  // returns the directional gradient of velocity for the vertical velocity field.

private:

};












#endif //WINDPARAMS_HPP
