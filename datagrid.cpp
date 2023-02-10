#include "datagrid.hpp"

// =============================================================================
// Constructors and Destructor
// =============================================================================

DataGrid::DataGrid(Disk disk, WindParams wind, int coordID, Array1d POS_1, Array1d POS_2) : disk(disk), wind(wind)
{

  this->coordID = coordID;

  param_nu = 0;
  param_p = 0;
  param_a = 0;
  param_b = 0;
  verbose = false;

  // set parameters
  if(coordID == SPHERICAL)
  {
    this->R = POS_1;
    this->Theta = POS_2;
    Head.push_back("r");
    Head.push_back("theta");
  }
  else if(coordID == CARTESIAN)
  {
    this->X = POS_1;
    this->Z = POS_2;
    Head.push_back("x");
    Head.push_back("z");
  }

  Array1d col1, col2;

  for(int i = 0; i < POS_1.size(); ++i)
  {
    for(int j = 0; j < POS_2.size(); ++j)
    {
      col1.push_back(POS_1[i]/disk.R_in);
      if(coordID == SPHERICAL)
      {
        col2.push_back(POS_2[j] * 180/PI);
      }
      else if(coordID == CARTESIAN)
      {
        col2.push_back(POS_2[j]/disk.R_in);
      }
    }
  }
  Res.push_back(col1);
  Res.push_back(col2);
}

// =============================================================================
// Sampling functions
// =============================================================================

Array1d DataGrid::Sample(int SampID, int funcID)
{
  Array1d col;
  if(coordID == SPHERICAL)
  {
    //loop over all (r, theta) observer points
    for(int i = 0; i < R.size(); ++i)
    {
      for(int j = 0; j < Theta.size(); ++j)
      {
        disk.setPolarObserver(R[i], Theta[j]);
        if(verbose)
          std::cout << R[i] / disk.R_in << " " << Theta[j] * 180/PI << std::endl;

        switch(SampID)
        {
        case OMEGA:
          col.push_back(disk.dOmega_integral(funcID, param_nu));
          break;
        case NU:
          col.push_back(disk.dnu_integral(funcID));
          break;
        case POWER:
          col.push_back(disk.power_integral(funcID, param_p));
          break;
        case BAND:
          if(param_b == 0)
          {
            col.push_back(disk.band_integral(funcID, param_a));
          }
          else
          {
            col.push_back(disk.band_integral(funcID, param_a, param_b));
          }
          break;
        case MEAN_PHOTON_ENERGY:
          col.push_back(disk.mean_photon_energy(funcID));
          break;
        case F_RAD:
          col.push_back(disk.F_rad(funcID, wind));
          break;
        }
      }
    }
  }
  if(coordID == CARTESIAN)
  {

    for(int i = 0; i < X.size(); ++i)
    {
      for(int j = 0; j < Z.size(); ++j)
      {
        disk.setCartesianObserver(X[i], Z[j]);
        if(verbose)
          std::cout << X[i] / disk.R_in << " " << Z[j] * disk.R_in << std::endl;

        switch(SampID)
        {
        case OMEGA:
          col.push_back(disk.dOmega_integral(funcID, param_nu));
          break;
        case NU:
          col.push_back(disk.dnu_integral(funcID));
          break;
        case POWER:
          col.push_back(disk.power_integral(funcID, param_p));
          break;
        case BAND:
          if(param_b == 0)
          {
            col.push_back(disk.band_integral(funcID, param_a));
          }
          else
          {
            col.push_back(disk.band_integral(funcID, param_a, param_b));
          }
          break;
        case MEAN_PHOTON_ENERGY:
          col.push_back(disk.mean_photon_energy(funcID));
          break;
        case F_RAD:
          col.push_back(disk.F_rad(funcID, wind));
          break;
        }
      }
    }
  }

  clearparams();

  return col;
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   nu     : array of nu values for the spectra
- Pushes the values for the spectral energy distribution for all r and theta.
- Will create a column for each frequency in the nu array passed in.
*/
void DataGrid::addSpectraAndExport(std::string ExportPath, int funcID, Array1d nu)
{
  int ResSize = Res.size();
  for(int i = 0; i < nu.size(); ++i)
  {
    std::string nu_str;
    std::ostringstream convert;
    convert << nu[i];
    Head.push_back(convert.str());
    Res.push_back(Array1d());
  }

  Array1d row;
  if(coordID == SPHERICAL)
  {
    for(int i = 0; i < R.size(); ++i)
    {
      for(int j = 0; j < Theta.size(); ++j)
      {
        disk.setPolarObserver(R[i], Theta[j]);
        if(verbose)
          std::cout << R[i] / disk.R_in << " " << Theta[j] * 180/PI << std::endl;

        row = disk.dOmega_integral(funcID, nu);
        for(int k = ResSize; k < ResSize + nu.size(); ++k)
        {
          Res[k].push_back(row[k - ResSize]);
        }
      }
    }
  }
  else if(coordID == CARTESIAN)
  {
    for(int i = 0; i < X.size(); ++i)
    {
      for(int j = 0; j < Z.size(); ++j)
      {
        disk.setPolarObserver(X[i], Z[j]);
        if(verbose)
          std::cout << X[i] / disk.R_in << " " << Z[j] * disk.R_in << std::endl;

        row = disk.dOmega_integral(funcID, nu);
        for(int k = ResSize; k < ResSize + nu.size(); ++k)
        {
          Res[k].push_back(row[k - ResSize]);
        }
      }
    }
  }

  Export(ExportPath);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   nu     : the value of nu to evaluate at
- Pushes the values for the geometric integral given nu for all r and theta
- into a results column
*/
void DataGrid::addOmegaIntegral(std::string name, int funcID, ldouble nu)
{
    if(verbose)
      std::cout << name << std::endl;
    param_nu = nu;
    Res.push_back(Sample(OMEGA, funcID));
    Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
- Pushes the values for the frequency integral for all r and theta
- into a results column
*/
void DataGrid::addFrequencyIntegral(std::string name, int funcID)
{
  if(verbose)
    std::cout << name << std::endl;
  Res.push_back(Sample(NU, funcID));
  Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   p      : exponent of nu for the power integral
- Pushes the values for the (nu^p)*f(nu) integral for all r and theta
- into a results column
*/
void DataGrid::addPowerIntegral(std::string name, int funcID, ldouble p)
{
  if(verbose)
    std::cout << name << std::endl;
  param_p = p;
  Res.push_back(Sample(POWER, funcID));
  Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   a      : lower bound of the frequency band
- Pushes the values for the frequency band integral for all r and theta
- into a results column
*/
void DataGrid::addBandIntegral(std::string name, int funcID, ldouble a)
{
  if(verbose)
    std::cout << name << std::endl;
  param_a = a;
  Res.push_back(Sample(BAND, funcID));
  Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   a      : lower bound of the frequency band
-   b      : upper bound of the frequency band
- Pushes the values for the frequency band integral for all r and theta
- into a results column
*/
void DataGrid::addBandIntegral(std::string name, int funcID, ldouble a, ldouble b)
{
  if(verbose)
    std::cout << name << std::endl;
  param_a = a;
  param_b = b;
  Res.push_back(Sample(BAND, funcID));
  Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
- Pushes the values for the mean photon energy for all r and theta
- into a results column
*/
void DataGrid::addMeanPhotonEnergy(std::string name, int funcID)
{
  if(verbose)
    std::cout << name << std::endl;
  Res.push_back(Sample(MEAN_PHOTON_ENERGY, funcID));
  Head.push_back(name);
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
- Pushes the values for the radiation force for all r and theta
- into a results column
*/
void DataGrid::addRadiationForce(std::string name, int funcID)
{
  if(verbose)
    std::cout << name << std::endl;
  Res.push_back(Sample(F_RAD, funcID));
  Head.push_back(name);
}

// =============================================================================
// Export functions
// =============================================================================

/*
- Exports the values in Res to a csv file with a header containing disk
- parameters, wind parameters, and column titles.
- This function will create a file but will NOT create a directory, so make
- sure that all folders in path already exist.
*/
void DataGrid::Export(std::string ExportPath)
{
  std::ofstream output;
  output.open(ExportPath.c_str());

  // header for data file with all parameters
  output << "Disk Parameters" << std::endl;
  output << "M" << "," << disk.M << std::endl;
  output << "M_dot" << "," << disk.M_dot << std::endl;
  output << "Gamma" << "," << disk.Gamma << std::endl;
  output << "epsilon" << "," << disk.epsilon << std::endl;
  output << "R_in" << "," << disk.R_in << std::endl;
  output << "R_out" << "," << disk.R_out << std::endl;
  output << "x" << "," << disk.x << std::endl;
  output << "beta" << "," << disk.beta << std::endl;
  output << "e_min" << "," << H * disk.nu_X_low << std::endl;
  output << "e_max" << "," << H * disk.nu_X_high << std::endl;
  output << std::endl;
  output << "Wind Parameters";
  output << "k" << "," << wind.k << std::endl;
  output << "alpha" << "," << wind.alpha << std::endl;
  output << "v_th" << "," << wind.v_th << std::endl;
  output << "rho" << "," << wind.rho << std::endl;
  output << "N_0" << "," << wind.N_0 << std::endl;
  output << "sigma_ref" << "," << wind.sigma_ref << std::endl;
  output << "eta_max" << "," << wind.eta_max << std::endl;
  output << "sz" << "," << wind.sz << std::endl;
  output << std::endl;
  output << "Data" << std::endl;

  // data head has column titles
  output << Head[0];
  for(int i = 1; i < Head.size(); ++i)
  {
    output << "," << Head[i];
  }
  output << std::endl;

  // put data into csv
  for(int i = 0; i < Res[0].size(); ++i)
  {
    output << Res[0][i];
    for(int j = 1; j < Res.size(); ++j)
    {
      output << "," << Res[j][i];
    }
    output << std::endl;
  }

  output.close();
}

void DataGrid::clear()
{
  std::string head1 = Head[0];
  std::string head2 = Head[1];
  Head.clear();
  Head.push_back(head1);
  Head.push_back(head2);

  Array1d col1 = Res[0];
  Array1d col2 = Res[1];
  Res.clear();
  Res.push_back(col1);
  Res.push_back(col2);
}

void DataGrid::clearparams()
{
  param_nu = 0;
  param_p = 0;
  param_a = 0;
  param_b = 0;
}



















//
