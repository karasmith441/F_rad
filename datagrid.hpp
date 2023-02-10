#ifndef DATAGRID_HPP
#define DATAGRID_HPP

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

#include "consts.hpp"
#include "integrator.hpp"
#include "disk.hpp"
#include "windparams.hpp"

#define ldouble long double
#define Array1d std::vector<ldouble>
#define Array2d std::vector<std::vector<ldouble> >
#define Array3d std::vector<std::vector<std::vector<ldouble> > >

#ifndef SSTR
#define SSTR( x ) static_cast< std::ostringstream & > ( std::ostringstream() << std::dec << x ).str()
#endif //SSTR

// ID's for coordinates
#define SPHERICAL 0
#define CARTESIAN 1

// ID's of disk functions to sample
#define SPECTRUM 0
#define OMEGA 1
#define NU 2
#define BAND 4
#define POWER 8
#define MEAN_PHOTON_ENERGY 16
#define F_RAD 32

// ID's of functions that can be evaluated in the integration functions
#define J_NU 0      // Mean intensity
#define F_R 1       // Radial component of Flux/Radiation Force
#define F_THETA 2   // Theta component of Flux/Radiation Force
//#define F_PHI 4     // Phi component of Flux/Radiation Force (commented out as it always goes to 0 for an axisymmetric disk)
#define F_DIST 8    // Flux as seen be a distant observer scaled to the set observer.

/*==============================================================================
DataGrid is a 2d array that stores a row for each point (r, theta) that is
sampled. It can be exported as a csv with a header containing all the
parameters, followed by a row with names for the variables calculated, followed
by the data. Variables are added in columns using the various add[...]
functions, and then the DataGrid is exported using Export.
==============================================================================*/
class DataGrid
{
public:

  Disk disk;                      // disk object handles disk calculations
  WindParams wind;                // stores wind parameters and computes optical depth
  int coordID;                    // determines if we are in polar or cartesian coords
  Array1d R, Theta, X, Z;         // Lists of positions

  ldouble param_nu, param_p, param_a, param_b;
  bool verbose;

  std::vector<std::string> Head;  // List of column names
  Array2d Res;                    // Array to store all results before export

  // Constructor
  DataGrid(Disk disk, WindParams wind, int coordID, Array1d POS_1, Array1d POS_2);

  void addOmegaIntegral(int funcID, ldouble nu) {addOmegaIntegral("geometric integral", funcID, nu); }
  void addFrequencyIntegral(int funcID) {addFrequencyIntegral("frequency integral", funcID); }
  void addPowerIntegral(int funcID, ldouble p) {addPowerIntegral("power integral", funcID, p); }
  void addBandIntegral(int funcID, ldouble a) {addBandIntegral("band integral", funcID, a); }
  void addBandIntegral(int funcID, ldouble a, ldouble b) {addBandIntegral("band integral", funcID, a, b); }
  void addMeanPhotonEnergy(int funcID) {addMeanPhotonEnergy("mean photon energy", funcID); }
  void addRadiationForce(int funcID) {addRadiationForce("radiation force", funcID); }

  void addOmegaIntegral(std::string name, int funcID, ldouble nu);              // Adds a column for integral over solid angle to compute moments of specific intensity
  void addFrequencyIntegral(std::string name, int funcID);                      // Adds a column for integral moments of specific intensity over frequency to get moments of intensity
  void addPowerIntegral(std::string name, int funcID, ldouble p);               // Adds a column for integral nu^p * moments of specific intensity
  void addBandIntegral(std::string name, int funcID, ldouble a);                // Adds a column for integral moments of specific intensity over nu in [a, inf)
  void addBandIntegral(std::string name, int funcID, ldouble a, ldouble b);     // Adds a column for integral moments of specific intensity over nu in [a, b]
  void addMeanPhotonEnergy(std::string name, int funcID);                       // Adds a column for mean photon energy weighted by a moment of specific intensity
  void addRadiationForce(std::string name, int funcID);                         // Adds a column for the force due to line absoption for a statistical distribution of lines

  void addSpectraAndExport(std::string ExportPath, int funcID, Array1d nu);     // Integrates over solid angle to compute moments of specific intensity

  // Exports to csv file
  // This function will create a file but will NOT create a directory, so make
  // sure that all folders in path already exist.
  void Export(std::string ExportPath);
  void clear();
  void clearparams();

private:

  // Does the actual loop over all R and Theta values to add a column
  Array1d Sample(int SampID, int funcID);

};








#endif //DATAGRID_HPP
