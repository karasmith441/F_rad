#ifndef DISK_HPP
#define DISK_HPP

#include <iostream>
#include <vector>
#include <stdexcept>

#include "consts.hpp"
#include "integrator.hpp"
#include "windparams.hpp"

#define ldouble long double
#define Array1d std::vector<ldouble>
#define Array2d std::vector<std::vector<ldouble> >
#define Array3d std::vector<std::vector<std::vector<ldouble> > >

// ID's of functions that can be evaluated in the integration functions
#define J_NU 0      // Mean intensity
#define F_R 1       // Radial component of Flux/Radiation Force
#define F_THETA 2   // Theta component of Flux/Radiation Force
//#define F_PHI 4     // Phi component of Flux/Radiation Force (commented out as it always goes to 0 for an axisymmetric disk)
#define F_DIST 8    // Flux as seen be a distant observer scaled to the set observer.

/*==============================================================================
This class handles storage of disk parameters and calculates a number of
spectral quantities of interest. It does this using a mix of Gaussian and
Gauss-Laguerre quadrature. When doing any calculations, you must set the
observation point before running the calculations as near a disk, many spectral
variables become extremely dependent on position.
==============================================================================*/
class Disk
{
public:

    // Input Parameters
    ldouble M, Gamma, epsilon, R_star, R_in, R_out, x, beta, nu_X_low, nu_X_high;
    // Computed Parameters
    ldouble M_dot, L_D, A_0;
    // position
    ldouble r0, theta0, x0, z0;                                                 // position in polar (r, theta) and cartesian (x, z)
    ldouble ct, st, tt;                                                         // cos, sin, and tan of theta0

    // Quadrature Memory                                                        // Numbers of quadrature points when integrating over ...
    int nu_pts;                                                                 // ... all frequencies
    int mu_pts;                                                                 // ... line of sight angle mu
    int phi_pts;                                                                // ... lign of sight angle phi_1
    int band_pts;                                                               // ... frequencies in a band
    int dist_pts;                                                               // ... radii for the distant observer approximation
    int star_pts;                                                               // ... frequencies in the corona spectrum

    Disk(ldouble M, ldouble Gamma, ldouble epsilon, ldouble R_star, ldouble R_in, ldouble R_out, ldouble x, ldouble beta, ldouble nu_X_low, ldouble nu_X_high);

    // void setQuadrature(int qID, int num_pts);                                // sets the number of quadrature points
    void setPolarObserver(ldouble r0, ldouble theta0);                          // sets the observation location using polar coords
    void setCartesianObserver(ldouble x0, ldouble z0);                          // sets the observation location using cartesian coords

    ldouble dOmega_integral(int funcID, ldouble nu);                            // Integrates over solid angle to compute moments of specific intensity
    Array1d dOmega_integral(int funcID, Array1d nu);                            // Same as above but for a whole spectrum
    ldouble dnu_integral(int funcID);                                           // Integrates moments of specific intensity over frequency to get moments of intensity
    ldouble power_integral(int funcID, ldouble p);                              // Integrates nu^p * moments of specific intensity
    ldouble band_integral(int funcID, ldouble a);                               // Integrates moments of specific intensity over nu in [a, inf)
    ldouble band_integral(int funcID, ldouble a, ldouble b);                    // Integrates moments of specific intensity over nu in [a, b]
    ldouble mean_photon_energy(int funcID);                                     // Computes mean photon energy weighted by a moment of specific intensity
    ldouble F_rad(int funcID, WindParams wind);                                 // Computes the force due to line absoption for a statistical distribution of lines

private:

  // region handling
  int region, int_n;                                                            // int_n is the number of integrals in a region
  ldouble mu_i, mu_o1, mu_o2, mu_star, mu_disk;                                 // Important mu bounds we can precompute

  // Quadrature Memory
  ldouble scale;                                                                // scales nu so that Laguerre integration converges nicely
  Array1d nu_x, nu_w;                                                           // abscissus and weight arrays for nu
  Array1d mu_x, mu_w;                                                           // ... for mu
  Array2d phi_x, phi_w;                                                         // ... for phi_1
  Array1d eta_x, eta_w;                                                         // ... for eta
  Array1d band_x, band_w;                                                       // ... for nu in a frequency band
  Array1d dist_x, dist_w;                                                       // ... for disk radii
  Array1d star_x, star_w;                                                       // ... for corona nu

  // Pre-computations arrays for integral sums
  Array2d R_grid;                                                               // array to store r_D for all mu and phi_1
  Array2d T_grid;                                                               // array to store temperature for all r_D
  Array2d t_grid;                                                               // array to store optical depth parameter t for all mu and phi_1
  Array2d Mt_grid;                                                              // array to store M(t) for all mu and phi_1
  Array3d I_grid;                                                               // array to store intensity for all mu, phi_1, and nu in quadrature
  Array3d I_spec_grid;                                                          // array to store intensity for all mu, phi_1, and nu in spectrum
  Array3d Opt_grid;                                                             // array to store (1 - e^(-eta t)) / t for all eta, mu, and phi_1

  Integrator quad;                                                              // Numerical integrator

  ldouble R_D(ldouble mu, ldouble phi_1);                                       // Computes the radius at which the line of sight hits the disk
  ldouble T_D(ldouble r_D);                                                     // Computes the temperature of the disk at some radius
  ldouble _B_(ldouble nu, ldouble T_r);                                         // Blackbody specific intensity at temperature T_r
  ldouble I_D(ldouble nu, ldouble T_r);                                         // Specific intensity of the disk from a point with temperature T_r

  ldouble L_star(ldouble nu);                                                   // Specific luminosity of the corona
  ldouble I_star(ldouble nu);                                                   // Speficic intensity of the corona

  void setOmegaGrid(int m);                                                     // Sets the quadrature grid for mu and phi_1
  ldouble OmegaGridHelper(int phiID, ldouble mu);                               // Computers the mu and phi_1 integration bounds

  void setRGrid(int m);                                                         // Precomputes r_D for every (mu, phi_1) pair
  void setTGrid(int m);                                                         // Precomputes T_D for every r_D
  void setSpecIntegrandGrid(int m, Array1d nu);                                 // Precomputes intensity for every (nu, mu, phi_1) triple
  void setIntegrandGrid(int m, ldouble p);                                      // Precomputes modified intensity for every (nu, mu, phi_1) triple
  void settGrid(int m, WindParams wind);                                        // Precomputes the optical depth parameter for every (mu, phi_1) pair

  void setRegion();                                                             // sets the region to keep track of what equations/integrals we use
  ldouble* getMuBounds(int m);                                                  // gets the mu boundaries
  int* getPhiID(int m);                                                         // gets the ID for computing phi boundaries

  ldouble F_distant(ldouble nu);                                                // Helper function for distant observer flux
  ldouble F_distant_gl(ldouble nu);                                             // Helper function for distant observer flux

};





#endif //DISK_HPP
