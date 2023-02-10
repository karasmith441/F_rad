#include "windparams.hpp"

WindParams::WindParams(ldouble M, ldouble k, ldouble alpha, ldouble v_th, ldouble rho, ldouble N_0, ldouble sigma_ref, ldouble eta_max, ldouble sz)
{
  // set parameters
  this->alpha = alpha;
  this-> v_th = v_th;
  this->rho = rho;
  this->N_0 = N_0;
  this->k = k;
  this->sigma_ref = sigma_ref;
  this->eta_max = eta_max;
  this->sz = sz;
}

/*
- Parameters
-   (r0, theta0) : polar coordinates of observer
-   (x0, z0)     : cartesian coordinates of the observer
-   ct, st, tt   : cosine, sine, and tangent of theta0
- Sets all of the positional variables for the point being observed, take these
- parameters from the disk.
*/
void WindParams::setObserver(ldouble r0, ldouble theta0, ldouble x0, ldouble z0, ldouble ct, ldouble st, ldouble tt)
{
  // set position variables
  this->r0 = r0;
  this->theta0 = theta0;
  this->x0 = x0;
  this->z0 = z0;
  this->ct = ct;
  this->st = st;
  this->tt = tt;
}

/*
Parameters:
- (mu, phi_1) : describes the line of sight.
- Returns the optical depth parameter t for a given line of sight
*/
ldouble WindParams::get_t(ldouble mu, ldouble phi_1)
{
  return sigma_ref * v_th * rho / fabs(dv_vertical(mu, phi_1));
}

/*
Parameters:
- (mu, phi_1) : describes the line of sight.
- Returns the gradient of velocity along the line of sight
- given line of sight unit vector n; return n dot grad(n dot v)
- This is computed for a vertical velocity flow v = v_0 * (z/r_0)^sz
*/
ldouble WindParams::dv_vertical(ldouble mu, ldouble phi_1)
{
  // components of n_hat (light of sight unit vector)
  ldouble nr = mu;                                                              // radial component
  ldouble nt = sqrt(1 - mu*mu) * cos(phi_1);                                    // theta component
  ldouble np = sqrt(1 - mu*mu) * sin(phi_1);                                    // phi component

  // The equation has be rearranged so that all radial dependence is
  // pulled out front
  ldouble rad = sqrt(G*M / (r0 * st)) * ct / pow(tt, sz) * 1 / (2*r0);

  // The geometric factor of dv is a sum over each pair of n_hat components
  // which come from n dot grad(n dot v)
  ldouble geo = nr*nr                                                           // radial component squared
              + nt*nt*(1 - 2*sz / (ct*ct))                                      // theta component squared
                                                                                // the factor for the phi squared component is 0 so it is left out
              + nr*nt*((1 + 2*sz) / (ct * st))                                  // radial/theta cross component
              + nr*np*(pow(tt, sz) / ct)                                        // radial/phi cross component
              + nt*np*(pow(tt, sz) / st);                                       // theta/phi cross component

  // and for the full result we just take rad * geo
  return rad * geo;
}
