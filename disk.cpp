#include "disk.hpp"

// =============================================================================
// Constructors and Destructor
// =============================================================================

Disk::Disk(ldouble M, ldouble Gamma, ldouble epsilon, ldouble R_star, ldouble R_in, ldouble R_out, ldouble x, ldouble beta, ldouble nu_X_low, ldouble nu_X_high)
{
    // set all parameters
    this->M = M;
    this->Gamma = Gamma;
    this->epsilon = epsilon;
    this->R_star = R_star;
    this->R_in = R_in;
    this->R_out = R_out;
    this->beta = beta;
    this->nu_X_low = nu_X_low;
    this->nu_X_high = nu_X_high;

    // initialize the integrator
    quad = Integrator();
    nu_pts = 100;
    mu_pts = 12;
    phi_pts = 12;
    band_pts = 50;
    dist_pts = 50;
    star_pts = 50;

    // Compute accretion rate and total disk luminosity
    M_dot = 4 * PI * MP * G * M / (epsilon * C * SIGMA_T) * Gamma;
    L_D = G * M * M_dot / (2*R_in);

    // Set scale for Gauss-Laguerre quadrature
    ldouble max_T = sqrt(sqrt(3 * G * M * M_dot / (8 * PI * PI * R_in * R_in * R_in * SIGMA) * 0.0566528));
    scale = K * max_T / (2*H);

    // Set the power law normalization to 1
    A_0 = 1;
    this->x = 1;

    // Integrate the powerlaw over the range of all nu
    ldouble total = 0;
    quad.gauss(nu_pts, 0, nu_X_low, nu_X_high, star_x, star_w);
    for(int i = 0; i < star_pts; ++i)
      total += L_star(star_x[i]) * fabs(star_w[i]);

    // Adjust the normalization so that when L_star is integrated over nu
    // it returns x * L_D
    this->x = x;
    A_0 = L_D / total;
}

// =============================================================================
// Integral Bounds Handlers
// =============================================================================

/*
- Parameters
-   (r0, theta0) : polar coordinates of observer
- Sets all of the positional variables for the point being observed
- takes the position in spherical coords
*/
void Disk::setPolarObserver(ldouble r0, ldouble theta0)
{
  this->r0 = r0;
  this->theta0 = theta0;

  // compute trig values
  ct = cos(theta0);
  st = sin(theta0);
  tt = tan(theta0);

  // compute cartesian coords
  x0 = r0 * st;
  z0 = r0 * ct;

  // set the region so we know what integrals to use
  setRegion();
}

/*
- Parameters
-   (x0, z0) : polar coordinates of observer
- Sets all of the positional variables for the point being observed
- takes the position in cartesian coords
*/
void Disk::setCartesianObserver(ldouble x0, ldouble z0)
{
  this->x0 = x0;
  this->z0 = z0;

  // set polar coordinates
  r0 = sqrt(x0*x0 + z0*z0);
  if(z0 != 0)
  {
    theta0 = atan(x0 / z0);
  }
  else if(x >= 0)
  {
    theta0 = PI / 2.0;
  }
  else
  {
    theta0 = -PI / 2.0;
  }

  // set trig values
  ct = cos(theta0);
  st = sin(theta0);
  tt = tan(theta0);

  // set the region so we know what integrals to use
  setRegion();
}

/*
- The number of integrals and the equations for computing their bounds
- vary depending on the observer's position. We give each equation/# of
- integrals a region ID. This function computes that ID and sets the current
- region for the disk to it.
*/
void Disk::setRegion()
{
  // All lines of sight for the special mu values are on the phi_1 = 0 plane
  ldouble R_i = sqrt(r0*r0 + R_in*R_in + 2*r0*R_in*st);                         // distance to close inner edge of disk
  ldouble R_o1 = sqrt(r0*r0 + R_out*R_out + 2*r0*R_out*st);                     // distance to close outer edge of disk
  ldouble R_o2 = sqrt(r0*r0 + R_out*R_out - 2*r0*R_out*st);                     // distance to far outer edge of disk

  mu_i = (r0*r0 + R_i*R_i - R_in*R_in)/(2*r0*R_i);                              // mu at which line of sight is hits close inner edge of disk
  mu_o1 = (r0*r0 + R_o1*R_o1 - R_out*R_out)/(2*r0*R_o1);                        // mu at which line of sight is hits close outer edge of disk
  mu_o2 = (r0*r0 + R_o2*R_o2 - R_out*R_out)/(2*r0*R_o2);                        // mu at which line of sight is hits far outer edge of disk

  // compute mu at which the light of sight is tangent to the corona
  if(r0 >= R_in)
  {
    mu_star = sqrt(r0*r0 - R_in*R_in) / r0;
  }
  else
  {
    mu_star = 0; // if mu_star is imaginary it won't be used but set it to zero to avoid having an "undefined variable"
  }
  // compute mu at which the light of sight is tangent to a sphere
  // with radius R_out and center at the center of the disk
  if (r0 >= R_out)
  {
    mu_disk = sqrt(r0*r0 - R_out*R_out) / r0;
  }
  else
  {
    mu_disk = 0; // if mu_disk is imaginary it won't be used but set it to zero to avoid having an "undefined variable"
  }

  if (x0 < R_in && mu_star < mu_o2)                          { region = -1; int_n = 0;}   // Disk is completely obscured by corona
  else if (x0 == 0)                                          { region = 0; int_n = 1;}    // Point is directly above corona
  else if (x0 <= R_out && x0 <= R_in && mu_star > mu_o1)     { region = 1; int_n = 2;}    // Point is above star, shadow doesn't reach disk edge
  else if (x0 <= R_out && x0 <= R_in)                        { region = 2; int_n = 1;}    // Point is above star, shadow reaches disk edge
  else if (x0 <= R_out && mu_star > mu_o1)                   { region = 3; int_n = 3;}    // Point is above disk, shadow doesn't reach disk edge
  else if (x0 <= R_out)                                      { region = 4; int_n = 2;}    // Point is above disk, shadow reaches disk edge
  else if (mu_star > mu_o1 && mu_star > mu_o2)               { region = 5; int_n = 4;}    // Point is not above disk, shadow reaches disk edge
  else if (mu_star > mu_o2)                                  { region = 6; int_n = 3;}    // Point is not above disk, shadow reaches disk edge
  else                                                       { region = 7; int_n = 3;}    // Point is not above disk, virtual shadow reaches disk edge

}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- Returns the mu bounds for the mth integral in the disk's
- calculations.
*/
ldouble* Disk::getMuBounds(int m)
{

    ldouble* muBounds = new ldouble[2];

    if (m >= int_n)
    {
      std::cout << "INVALID ARGUMENT: Region " << region << " has " << int_n << " integrals; integral requested must have id from 0 to " << int_n << ", but you have requested id " << m << std::endl;
      muBounds[0] = 0; muBounds[1] = 0;
      return muBounds;
    }

    // Each region ID has a list of mu bounds indexed by m
    switch(region)
    {
      case -1:
        muBounds[0] = 0; muBounds[1] = 0;
        break;
      case 0:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_o1; }
        break;
      case 1:
        if (m == 0) {muBounds[0] = mu_star; muBounds[1] = mu_o1; }
        if (m == 1) {muBounds[0] = mu_o1; muBounds[1] = mu_o2; }
        break;
      case 2:
        if (m == 0) {muBounds[0] = mu_star; muBounds[1] = mu_o2; }
        break;
      case 3:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_star; }
        if (m == 1) {muBounds[0] = mu_star; muBounds[1] = mu_o1; }
        if (m == 2) {muBounds[0] = mu_o1; muBounds[1] = mu_o2; }
        break;
      case 4:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_star; }
        if (m == 1) {muBounds[0] = mu_star; muBounds[1] = mu_o2; }
        break;
      case 5:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_star; }
        if (m == 1) {muBounds[0] = mu_star; muBounds[1] = mu_o2; }
        if (m == 2) {muBounds[0] = mu_o1; muBounds[1] = mu_o2; }
        if (m == 3) {muBounds[0] = mu_o2; muBounds[1] = mu_disk; }
        break;
      case 6:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_star; }
        if (m == 1) {muBounds[0] = mu_star; muBounds[1] = mu_o2; }
        if (m == 2) {muBounds[0] = mu_o2; muBounds[1] = mu_disk; }
        break;
      case 7:
        if (m == 0) {muBounds[0] = mu_i; muBounds[1] = mu_o2; }
        if (m == 1) {muBounds[0] = mu_o2; muBounds[1] = mu_star; }
        if (m == 2) {muBounds[0] = mu_star; muBounds[1] = mu_disk; }
        break;
      default:
        muBounds[0] = 0; muBounds[1] = 0;
    }

    // bounds must always be lower bound < upper bound
    // this code should not be executed, but it is here as
    // a failsafe anyways.
    if (muBounds[0] > muBounds[1])
    {
      ldouble temp = muBounds[1];
      muBounds[1] = muBounds[0];
      muBounds[0] = temp;
    }

    return muBounds;
}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- Returns the ID's for the phi_1 bounds for the mth integral
- in the disk's calculations. These ID's are used to compute actual phi_1
- boundaries in the OmegaGridHelper function
*/
int* Disk::getPhiID(int m)
{

  // 0 - 0
  // 1 - phi_s
  // 2 - phi_d1
  // 3 - phi_d2
  // 4 - pi

  int* phiID = new int[2];

  switch(region)
  {
    case -1:
      phiID[0] = 0; phiID[1] = 0;
      break;
    case 0:
      if (m == 0) {phiID[0] = 0; phiID[1] = 4; }
      break;
    case 1:
      if (m == 0) {phiID[0] = 0; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 2; phiID[1] = 4; }
      break;
    case 2:
      if (m == 0) {phiID[0] = 2; phiID[1] = 4; }
      break;
    case 3:
      if (m == 0) {phiID[0] = 1; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 0; phiID[1] = 4; }
      if (m == 2) {phiID[0] = 2; phiID[1] = 4; }
      break;
    case 4:
      if (m == 0) {phiID[0] = 1; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 2; phiID[1] = 4; }
      break;
    case 5:
      if (m == 0) {phiID[0] = 1; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 0; phiID[1] = 4; }
      if (m == 2) {phiID[0] = 2; phiID[1] = 4; }
      if (m == 3) {phiID[0] = 2; phiID[1] = 3; }
      break;
    case 6:
      if (m == 0) {phiID[0] = 1; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 2; phiID[1] = 4; }
      if (m == 2) {phiID[0] = 2; phiID[1] = 3; }
      break;
    case 7:
      if (m == 0) {phiID[0] = 1; phiID[1] = 4; }
      if (m == 1) {phiID[0] = 1; phiID[1] = 3; }
      if (m == 2) {phiID[0] = 2; phiID[1] = 3; }
      break;
    default:
      phiID[0] = 0; phiID[1] = 0;
  }

  return phiID;

}

// =============================================================================
// Geometric and Radiative Functions
// =============================================================================

/*
- Parameters:
-   (mu, phi_1) : describes the line of sight.
- Returns the radius, r_D, at which a line of sight from the observer
- intersects with the disk.
*/
ldouble Disk::R_D(ldouble mu, ldouble phi_1)
{
  ldouble st_1 = sqrt(1 - mu * mu);                                             // sin(theta_1)
  ldouble R_1 = r0 * ct / (mu * ct - st_1 * cos(phi_1) * st);                   // distance from the point in the wind to the disk along the line of sight
  ldouble r_D = sqrt(R_1 * R_1 + r0 * r0 - 2 * R_1 * r0 * mu);                  // use the law of cosines to find r_D from R_1

  return r_D;
}

/*
- Parameters:
-   r_D : radius along the plane of the Disk
- Returns the temperature of a Shakura Sunyaev disk with irradiation from a
- central object as a function of disk radius.
*/
ldouble Disk::T_D(ldouble r_D)
{

  // If r_D is not on the disk then the temperature is 0
  if(r_D < R_in || r_D > R_out)
  {
    return 0;
  }

  ldouble RR = R_in / r_D;                                                      // 1/r in units of R_in
  ldouble coeff = 3 * G * M * M_dot / (8 * PI * PI * R_in * R_in * R_in);       // Shakura Sunyaev temperature coefficient

  ldouble ss = (RR*RR*RR) * (1 - sqrt(RR));                                     // Shakura Sunyaev temperature function
  ldouble ir = x / (3 * PI) * (asin(RR) - RR * sqrt(1 - RR*RR));                // Contribution from irradiation

  ldouble I_D = coeff * (ss + ir);                                              // Compute the integrated disk intensity

  return sqrt(sqrt(I_D / SIGMA));                                               // I = sigma T^4 --> T = (sigma/I)^1/4
}

/*
- Parameters:
-   nu  : frequency
-   T_r : temperature
- Returns the blackbody radiation at a frequency from a source with
- temperature T_r
*/
ldouble Disk::_B_(ldouble nu, ldouble T_r)
{

  // 0 temperature blackbody does not radiate
  if(T_r == 0)
  {
    return 0;
  }

  ldouble coeff = 2 * H / (C * C);
  ldouble ex = exp(H * nu / (K * T_r));

  if (ex - 1 == 0 || isinf(ex))                                                 // If there is an overflow or divByZero then set this to 0 because the
  {                                                                             // only cases in which the code would do that is when B is tending
    return 0;                                                                   // toward 0 anyways
  }

  ldouble B = coeff * nu*nu*nu / (ex - 1);

  return B;
}

/*
- Parameters:
-   nu  : frequency
-   T_r : temperature
- Returns a modified blackbody that can be use Gauss-Laguerre quadrature
- for integration over frequency
*/
ldouble Disk::I_D(ldouble nu, ldouble T_r)
{
  // 0 temperature blackbody does not radiate
  if(T_r == 0)
  {
    return 0;
  }

  // Modified version of Planck function is scaled so that it extends into the
  // range where Gauss-Laguerre quadrature is numerically stable. The four
  // factors of nu come from the nu^3 term and the dnu term in the integration
  ldouble coeff = 2 * H / (C * C) * scale*scale*scale*scale;
  ldouble ex1 = exp(-nu);                                                       // Gauss-Laguerre has an exponential factor, but for stability
                                                                                // I use a negative exponential and put it in the denominator
  ldouble ex2 = exp(H*nu*scale/(K*T_r)) - 1;                                    // exponent with scaled nu minus 1 (denominator of modified Planck function)
  ldouble ex = ex1 * ex2;

  if (ex == 0 || isinf(ex))                                                     // If there is a divByZero or overflow, then set this to 0 because the
  {                                                                             // only case in which the code would do that is when B is tending
    return 0;                                                                   // toward 0 anyways
  }

  return coeff / ex;
}

/*
- Parameters:
-   nu  : frequency
-   T_r : temperature
- Returns the specific luminosity of the corona
*/
ldouble Disk::L_star(ldouble nu)
{
    // As long as nu is within the accepted range
    if (nu_X_low <= nu && nu <= nu_X_high)
    {
        return x * A_0 * pow(nu, beta);  // return a scaled powerlaw
    }
    else
    {
        return 0; // corona does not radiate frequencies outside of range
    }
}

/*
- Parameters:
-   nu  : frequency
- Returns the specific intensity of the corona
*/
ldouble Disk::I_star(ldouble nu)
{
    return L_star(nu) / (4*PI*PI*R_star*R_star);
}

//==============================================================================
// Fill Quadrature Grids
//==============================================================================

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- This function gets the integral bounds for mu, computes the abscissa for
- the mu quadrature, and then computes the phi_1 bounds and abscissa for
- each mu
*/
void Disk::setOmegaGrid(int m)
{
    // get mu bounds and equation ID's for phi bounds
    ldouble* mu_bounds = getMuBounds(m);
    int* phiID = getPhiID(m);
    ldouble phi_lower, phi_upper, temp;

    // m == -1 is not handled by getMuBounds, but it denotes integration
    // over the corona
    if(m == -1)
    {
      mu_bounds[0] = mu_star;
      mu_bounds[1] = 1;
      phiID[0] = 0;
    }

    quad.gauss(mu_pts, 0, mu_bounds[0], mu_bounds[1], mu_x, mu_w);              // Compute the gaussian quadrature for mu

    // clear the phi_1 arrays
    phi_x.clear();
    phi_w.clear();

    // for each mu, compute phi_1 bounds
    for(int i = 0; i < mu_pts; ++i)
    {
        // m == -1 is denotes integration over the corona which only has two
        // cases for it's upper phiID
        if(m == -1)
        {
          phiID[1] = (mu_x[i] < mu_i) ? 1 : 4;
        }

        // get phi bounds at particular mu = mu_x[i]
        phi_lower = OmegaGridHelper(phiID[0], mu_x[i]);
        phi_upper = OmegaGridHelper(phiID[1], mu_x[i]);

        if(phi_upper < phi_lower)                                               // Make sure phi bounds are in the correct order to avoid negative integrals
        {
          temp = phi_upper;
          phi_upper = phi_lower;
          phi_lower = temp;
        }

        // Push an array into the phi_1 quadrature integrals to store
        // abscissa and weights for each mu
        phi_x.push_back(Array1d());
        phi_w.push_back(Array1d());
        quad.gauss(phi_pts, 0, phi_lower, phi_upper, phi_x[i], phi_w[i]);       // Compute the gaussian for phi_1 at mu = mu_x[i]
    }

    // free memory
    delete[] mu_bounds;
    delete[] phiID;
}

/*
- Parameters:
-   phiID : code for which phi calculation equation should be used
-   mu    : an abscissus for the mu quadrature
- This computes phi_1 integral bounds as a function of mu according to the
- equation specified by phiID
*/
ldouble Disk::OmegaGridHelper(int phiID, ldouble mu)
{
    ldouble R_x, st_1, det, D, R_lim, Alpha;
    bool d_flag;                                                                // flag determines sign of D

    switch(phiID)
    {
    case 0:                                                                     // phi = 0
      return 0;
    case 1:                                                                     // phi hits inner disk edge
      R_x = R_in;
      break;
    case 2:                                                                     // phi hits the far outer disk edge
      R_x = R_out;
      d_flag = true;
      break;
    case 3:                                                                     // phi hits the near outer disk edge
      R_x = R_out;
      d_flag = false;
      break;
    case 4:                                                                     // phi goes all the way to pi
      return PI;
    }

    st_1 = sqrt(1 - mu*mu);                                                     // sin(theta_1)

    det = R_x*R_x - r0*r0*st_1*st_1;                                            // determinant in quadratic formula form of law of cosines
    D = (det < 0) ? 0 : sqrt(det);                                              // det cannot be less than 0. If it is due to some error, snap it to 0

    R_lim = (d_flag) ? r0*mu + D : r0*mu - D;                                   // distance to ? far : near edge

    Alpha = (R_lim*mu - r0) / (tt*R_lim*st_1);                                  // Alpha = cos(phi_1) at boundary

    return (Alpha > 1) ? 0 : ((Alpha < -1) ? PI : acos(Alpha));                 // return the boundary phi_1. Alpha is always in [-1, 1]; if it isn't due to some error, snap it to [-1,1]
}

//==============================================================================
// Evaluating values over quadrature grids
//==============================================================================

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- populates mu, phi_1, and R_grid arrays.
*/
void Disk::setRGrid(int m)
{
    R_grid.clear();
    setOmegaGrid(m);

    for(int i = 0; i < mu_pts; ++i)
    {
      R_grid.push_back(Array1d());
      for(int j = 0; j < phi_pts; ++j)
      {
        // R_grid stores disk radii for each (mu, phi_1) pair
        R_grid[i].push_back(R_D(mu_x[i], phi_x[i][j]));
      }
    }
}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- populates mu, phi_1, R_grid and T_grid arrays.
*/
void Disk::setTGrid(int m)
{
  T_grid.clear();
  setRGrid(m);

  for(int i = 0; i < mu_pts; ++i)
  {
    T_grid.push_back(Array1d());
    for(int j = 0; j < phi_pts; ++j)
    {
      // T_grid stores disk temperature for each radius
      T_grid[i].push_back(T_D(R_grid[i][j]));
    }
  }
}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- populates mu, phi_1, R_grid, T_grid, and I_spec_grid arrays.
- This is used for spectra where we don't need modified intensity for
- Gauss-Laguerre quadrature
*/
void Disk::setSpecIntegrandGrid(int m, Array1d nu)
{
    I_spec_grid.clear();

    // m == -1 denotes integration over the corona
    if(m == -1)
    {
      setOmegaGrid(m);
      for(int i = 0; i < nu.size(); ++i)
      {
        // we can compute I up here because it doesn't depend on line of sight
        ldouble I = I_star(nu[i]);
        I_spec_grid.push_back(Array2d());
        for(int j = 0; j < mu_pts; ++j)
        {
          I_spec_grid[i].push_back(Array1d());
          for(int k = 0; k < phi_pts; ++k)
          {
            // I_spec_grid stores intensity for each (nu, mu, phi_1) triple
            I_spec_grid[i][j].push_back(I);
          }
        }
      }
    }
    else
    {
      // If we're looking at the disk we need to populate temperatures
      // and that will populate mu and phi_1 grids for us
      setTGrid(m);
      for(int i = 0; i < nu.size(); ++i)
      {
        I_spec_grid.push_back(Array2d());
        for(int j = 0; j < mu_pts; ++j)
        {
          I_spec_grid[i].push_back(Array1d());
          for(int k = 0; k < phi_pts; ++k)
          {
            // I_spec_grid stores intensity for each (nu, mu, phi_1) triple
            I_spec_grid[i][j].push_back(_B_(nu[i], T_grid[j][k]));
          }
        }
      }
    }
}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- populates mu, phi_1, R_grid, T_grid, and I_grid arrays. It
- uses I_D instead of _B_ because this grid will be used in Gauss-Laguerre
- quadrature.
*/
void Disk::setIntegrandGrid(int m, ldouble p)
{
    I_grid.clear();

    // m == -1 denotes integration over the corona
    if(m == -1)
    {

      // Compute gaussian quadrature for the corona frequencies
      quad.gauss(star_pts, 0, nu_X_low, nu_X_high, star_x, star_w);
      setOmegaGrid(m);
      for(int i = 0; i < star_pts; ++i)
      {
        // we can compute I up here because it doesn't depend on line of sight
        ldouble I = I_star(star_x[i]);
        I_grid.push_back(Array2d());
        for(int j = 0; j < mu_pts; ++j)
        {
          I_grid[i].push_back(Array1d());
          for(int k = 0; k < phi_pts; ++k)
          {
            // I_grid stores modified star intensity for each (nu, mu, phi_1) triple
            I_grid[i][j].push_back(I);
          }
        }
      }
    }
    else
    {
      // compute Gauss-Laguerre quadrature for frequencies
      quad.gaulag(nu_pts, p, nu_x, nu_w);
      // If we're looking at the disk we need to populate temperatures
      // and that will populate mu and phi_1 grids for us
      setTGrid(m);

      for(int i = 0; i < nu_pts; ++i)
      {
        I_grid.push_back(Array2d());
        for(int j = 0; j < mu_pts; ++j)
        {
          I_grid[i].push_back(Array1d());
          for(int k = 0; k < phi_pts; ++k)
          {
            // I_grid stores modified disk intensity for each (nu, mu, phi_1) triple
            I_grid[i][j].push_back(I_D(nu_x[i], T_grid[j][k]));
          }
        }
      }
    }
}

/*
- Parameters:
-   m : Most regions require multiple integrals. m says which integral we want.
- evaluates optical depth parameter t and force multiplier M(t) as a function
- of mu and phi_1 for all mu and phi_1 in the quadrature arrays, then stores
- those values in t_grid.
*/
void Disk::settGrid(int m, WindParams wind)
{
    // Only line of sight quadrature grid is needed for optical depth parameter
    setOmegaGrid(m);
    for(int i = 0; i < mu_pts; ++i)
    {
      for(int j = 0; j < phi_pts; ++j)
      {
        // t_grid stores optical depth parameter t for each (mu, phi_1) pair
        // Mt_grid stores the sight line force multiplier
        ldouble t = wind.get_t(mu_x[i], phi_x[i][j]);
        ldouble tau_max = wind.eta_max * t;
        t_grid[i][j] = t;
        Mt_grid[i][j] = wind.k * pow(t, -wind.alpha) * (pow(1 + tau_max, 1 - wind.alpha) - 1)/pow(tau_max, 1 - wind.alpha);
      }
    }
}

//==============================================================================
// Integration
//==============================================================================

/*
- Parameters:
-   nu : frequency
- Returns the flux as seen by a distant observer, scaled to the distance of
- the observer that was set by setPolarObserver or setCartesianObserver
*/
ldouble Disk::F_distant(ldouble nu)
{
    ldouble total = 0;
    for(int i = 0; i < dist_pts; ++i)
    {
        total += _B_(nu, T_D(dist_x[i])) * dist_x[i] * fabs(dist_w[i]);              // integrate R*I dR from R_in to R_out
    }

    std::cout << nu << " " << 2*PI * ct/(r0*r0) * total << std::endl;
    return 2*PI * ct/(r0*r0) * total;                                           // scale by 2pi * cos(theta)/r^2
}

/*
- Parameters:
-   nu : frequency
- Returns the flux as seen by a distant observer, scaled to the distance of
- the observer that was set by setPolarObserver or setCartesianObserver
- using I_D instead of _B_ so that it can be integrated with Gauss-Laguerre
- quadrature
*/
ldouble Disk::F_distant_gl(ldouble nu)
{
    ldouble total = 0;
    for(int i = 0; i < dist_pts; ++i)
    {
        total += I_D(nu, T_D(dist_x[i])) * dist_x[i] * fabs(dist_w[i]);                      // integrate R*I dR from R_in to R_out
    }

    return 2*PI*ct/(r0*r0) * total;                                             // scale by 2pi * cos(theta)/r^2
}

/*
- Parameters:
-    funcID : id of the function that you would like to evaluate
-    nu     : frequency
- Returns the specific frequency value of the specified function. Computed
- by integrating over the solid angle occupied by the disk.
- see disk.hpp for funcID's
*/
ldouble Disk::dOmega_integral(int funcID, ldouble nu)
{
  // This function just creates an array of length 1 and uses the array
  // version of dOmega_integral. Potential soruce of furture optimization.
  Array1d nu_vec;
  nu_vec.push_back(nu);

  return dOmega_integral(funcID, nu_vec)[0];
}

/*
- Parameters:
-    funcID : id of the function that you would like to evaluate
-    nu     : frequency array
- Returns an array of specific frequency values of the specified function.
- Computed by integrating over the solid angle occupied by the disk.
- see disk.hpp for funcID's
*/
Array1d Disk::dOmega_integral(int funcID, Array1d nu)
{
    Array1d integ_0;                                                            // phi_1 integral sum
    Array1d total;                                                              // total integral sum
    Array1d st1;                                                                // sin(theta_1)
    Array2d cph, sph;                                                           // cos(phi_1), sin(phi_1)
    int int_n_temp = int_n;

    // initialize total and integ_0 to all 0's
    for(int i = 0; i < nu.size(); ++i)
    {
      integ_0.push_back(0);
      total.push_back(0);
    }
    // if distant observer just compute F_distant
    if(funcID == F_DIST)
    {
      // compute the gaussian quadrature for distant F
      quad.gauss(dist_pts, 0, R_in, R_out, dist_x, dist_w);

      // compute distant F for each nu
      for(int i = 0; i < nu.size(); ++i)
        total[i] = F_distant(nu[i]);

      // corona distant F is identical to computation for F_R so switch that
      funcID = F_R;
      int_n = 0;
    }

    // start m at -1 only if there is radiation from the corona
    for(int m = (x != 0) ? -1 : 0; m < int_n; ++m)                              // For each integral id in the region
    {
      // set the integrand grid
      setSpecIntegrandGrid(m, nu);

      // Precompute trig function values
      st1.clear();
      cph.clear();
      sph.clear();
      for(int i = 0; i < mu_pts; ++i)
      {
        st1.push_back(sqrt(1 - mu_x[i]*mu_x[i]));
        cph.push_back(Array1d());
        sph.push_back(Array1d());
        for(int j = 0; j < phi_pts; ++j)
        {
          cph[i].push_back(cos(phi_x[i][j]));
          sph[i].push_back(sin(phi_x[i][j]));
        }
      }

      // loop over nu, mu, and phi_1
      for(int i = 0; i < nu.size(); ++i)
      {
        for(int j = 0; j < mu_pts; ++j)                                         // mu quadrature loop
        {
          std::fill(integ_0.begin(), integ_0.end(), 0);                         // zero out the counters
          for(int k = 0; k < phi_pts ; ++k)                                     // phi_1 quadrature loop
          {
            switch(funcID)                                                      // funcID determines whether we compute J_nu, F_r, F_theta, or F_phi
            {
              case J_NU:
                integ_0[i] += I_spec_grid[i][j][k]/(4*PI) * fabs(phi_w[j][k]);
                break;
              case F_R:
                integ_0[i] += mu_x[j] * I_spec_grid[i][j][k] * fabs(phi_w[j][k]);
                break;
              case F_THETA:
                integ_0[i] += -st1[j] * cph[j][k] * I_spec_grid[i][j][k] * fabs(phi_w[j][k]);
                break;
              // For axisymmetric disk F_PHI will go to zero so we just ignore it
              // case F_PHI:
              //   integ_0[i] += st1[i] * sin(phi_x[i][j]) * I_spec_grid[i][j][k] * fabs(phi_w[i][j]);
              //   break;
            }
          }
          total[i] += 2 * integ_0[i] * fabs(mu_w[j]);                           // double because integration is reflected over disk (axisymmetry)
        }
      }
    }

    // if int_n got set to zero for corona calculations, set it back to what it's supposed to be
    int_n = int_n_temp;
    return total;
}

/*
- Parameters:
-    funcID : id of the function that you would like to evaluate
- Returns the frequency integrated value of the specified functions by
- integrating first over solid angle using gaussian quadrature, then over
- frequency using Gauss-Laguerre quadrature. Frequency is taken over [0, inf)
- see disk.hpp for funcID's
*/
ldouble Disk::dnu_integral(int funcID)
{
    // this is just the special case of power_integral with no power of nu
    return power_integral(funcID, 0);
}

/*
- Parameters:
-    funcID : id of the function that you would like to evaluate
-    p      : power of the exponential nu factor
- Returns the frequency integrated value of (nu^p)*f(nu) where f is the
- specified functions. Computed by integrating first over solid angle using
- gaussian quadrature, then over frequency using Gauss-Laguerre quadrature.
- Frequency is taken over [0, inf)
- see disk.hpp for funcID's
*/
ldouble Disk::power_integral(int funcID, ldouble p)
{
  int int_n_temp = int_n;
  ldouble integ_0, integ_1, total;                                              // internal integral sums
  total = 0;
  Array1d st1;
  Array2d cph, sph;

  // The Planck function has a factor of nu^3 so we add this into the power.
  p += 3;

  if(funcID == F_DIST)
  {
    quad.gaulag(nu_pts, p, nu_x, nu_w);                                         // Computes Gauss-Laguerre quadrature for frequency
    quad.gauss(dist_pts, 0, R_in, R_out, dist_x, dist_w);                       // Computes Gaussian quadrature for radii

    // integral loop for distant F
    for(int i = 0; i < nu_pts; ++i)
    {
      total += F_distant_gl(nu_x[i]) * fabs(nu_w[i]);
    }

    // corona distant F is identical to computation for F_R so switch that
    funcID = F_R;
    int_n = 0;
  }

  // start m at -1 only if there is radiation from the corona
  for(int m = ((x != 0) ? -1 : 0); m < int_n; ++m)                              // For each integral id in the region
  {
    // set the Gauss-Laguerre integrand grid
    setIntegrandGrid(m, p);

    // Precompute trig function values
    st1.clear();
    cph.clear();
    sph.clear();
    for(int i = 0; i < mu_pts; ++i)
    {
      st1.push_back(sqrt(1 - mu_x[i]*mu_x[i]));
      cph.push_back(Array1d());
      sph.push_back(Array1d());
      for(int j = 0; j < phi_pts; ++j)
      {
        cph[i].push_back(cos(phi_x[i][j]));
        sph[i].push_back(sin(phi_x[i][j]));
      }
    }

    // loop over all frequencies for disk and corona frequencies for corona
    for(int i = 0; i < ((m != -1) ? nu_pts : star_pts); ++i)
    {
      integ_1 = 0;                                                              // zero mu counter
      for(int j = 0; j < mu_pts; ++j)                                           // mu quadrature loop
      {
        integ_0 = 0;                                                            // zero phi_1 counter
        for(int k = 0; k < phi_pts; ++k)                                        // phi_1 quadrature loop
        {
          switch(funcID)                                                        // funcID determines whether we compute J_nu, F_r, F_theta, or F_phi
          {
            case J_NU:
              integ_0 += I_grid[i][j][k]/(4*PI) * fabs(phi_w[j][k]);
              break;
            case F_R:
              integ_0 += mu_x[j] * I_grid[i][j][k] * fabs(phi_w[j][k]);
              break;
            case F_THETA:
              integ_0 += -st1[j] * cph[j][k] * I_grid[i][j][k] * fabs(phi_w[j][k]);
              break;
            // For axisymmetric disk F_PHI will go to zero so we just ignore it
            // case F_PHI:
            //   integ_0 += st1[j] * sph[j][k] * I_grid[i][j][k] * fabs(phi_w[j][k]);
            //   break;
          }
        }
        integ_1 += integ_0 * fabs(mu_w[j]);
      }
      total += 2 * integ_1 * fabs(((m != -1) ? nu_w : star_w)[i]);              // double because integration is reflected over disk (axisymmetry)
    }
  }

  // if int_n got set to zero for corona calculations, set it back to what it's supposed to be
  int_n = int_n_temp;

  // Need extra scaling for each factor of nu beyond the normal 3 from the Planck function
  return pow(scale, p - 3) * total;
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   a      : lower bound of the integral
- Returns the frequency integrated value of the specified functions by
- integrating first over solid angle using gaussian quadrature, then over
- frequency using Gauss-Laguerre quadrature. Frequency is taken over [a, inf)
- see disk.hpp for funcID's
*/
ldouble Disk::band_integral(int funcID, ldouble a)
{
  ldouble total_0 = power_integral(funcID, 0);                                  // Integrate from 0 to infinity using Gauss-Laguerre
  Array1d total_1;
  Array1d nu_vec;

  // compute Gaussian quadrature from 0 to a
  quad.gauss(band_pts, 0, 0, a, band_x, band_w);
  for(int i = 0; i < band_pts; ++i)
  {
    nu_vec.push_back(band_x[i]);
  }
  // Compute geometric integral over the band spectrum
  total_1 = dOmega_integral(funcID, nu_vec);
  // Integrate from 0 to a using Gaussian quadrature but instead of adding
  // to a total on each iteration, subtract from the [0,inf) total.
  // [0, inf) - [0,a] --> [a,inf)
  for(int i = 0; i < band_pts; ++i)
  {
    total_0 -= total_1[i] * fabs(band_w[i]);
  }

  return total_0;
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-   a      : lower bound of the band
-   b      : upper bound of the band
- Returns the frequency integrated value of the specified functions by
- integrating first over solid angle using gaussian quadrature, then over
- frequency using Gauss-Laguerre quadrature. Frequency is taken over [a, b]
- see disk.hpp for funcID's
*/
ldouble Disk::band_integral(int funcID, ldouble a, ldouble b)
{
  ldouble total;
  Array1d integ_0;
  Array1d nu_vec;

  // Compute Gaussian quadrature from a to b
  quad.gauss(band_pts, 0, a, b, band_x, band_w);
  for(int i = 0; i < band_pts; ++i)
  {
    nu_vec.push_back(band_x[i]);
  }
  // Compute geometric integral over the band spectrum
  integ_0 = dOmega_integral(funcID, nu_vec);

  // Intagration loop
  for(int i = 0; i < band_pts; ++i)
  {
    total += integ_0[i] * fabs(band_w[i]);
  }

  return total;
}

/*
- Parameters:
-   funcID : id of the weighting function f(nu)
- Computes the mean photon energy <h nu> = int[h*nu*f(nu) dnu] / int[f(nu) dnu]
*/
ldouble Disk::mean_photon_energy(int funcID)
{
    // integral [h*nu*f(nu) dnu] is just h times power_integral with p = 1
    ldouble num = H * power_integral(funcID, 1);

    // integral [f(nu) dnu] is just power_integral with p = 1
    ldouble den = power_integral(funcID, 0);

    return num / den;
}

/*
- Parameters:
-   funcID : id of the function that you would like to evaluate
-     Note this is restricted to F_R and F_THETA
- Returns the specified component of the radiation force due to line absorption
- given a wind near the disk.
*/
ldouble Disk::F_rad(int funcID, WindParams wind)
{
    // set observation point in the wind at our observation point
    wind.setObserver(r0, theta0, x0, z0, ct, st, tt);

    ldouble integ_0, integ_1, integ_2, integ_3, total, J_tot;
    ldouble nu_proper;
    Array1d st1;
    Array2d cph, sph;

    Array1d nu_vec_x, nu_vec_w;
    Array1d J_nu;
    total = 0;

    // Compute total mean intensity
    J_tot = 4 * PI * power_integral(J_NU, 0);

    // start m at -1 only if there is contribution from the corona
    // loop over all integrals in the region
    for(int m = (x != 0 ? 1 : 0); m < int_n; ++m)
    {
      // set the quadrature grids for the Gauss-Laguerre integrands and optical depth parameter
      setIntegrandGrid(m, 3);
      settGrid(m, wind);

      // Precompute trig function values
      st1.clear();
      cph.clear();
      sph.clear();
      for(int i = 0; i < mu_pts; ++i)
      {
        st1.push_back(sqrt(1 - mu_x[i]*mu_x[i]));
        cph.push_back(Array1d());
        sph.push_back(Array1d());
        for(int j = 0; j < phi_pts; ++j)
        {
          cph[i].push_back(cos(phi_x[i][j]));
          sph[i].push_back(sin(phi_x[i][j]));
        }
      }

      // loop over all frequencies for disk and corona frequencies for corona
      for(int i = 0; i < (m != -1 ? nu_pts : star_pts); ++i)
      {
        integ_3 = 0;
        integ_2 = 0;
        for(int j = 0; j < mu_pts; ++j)
        {
          integ_1 = 0;
          integ_0 = 0;
          for(int k = 0; k < phi_pts; ++k)
          {
            // Use un-scaled nu to compute J_nu
            nu_proper = nu_x[i] * scale;

            // Integration of J_nu using _B_ intead of I_D because we want
            // Gaussian quadrature, not Gauss-Laguerre
            integ_0 += nu_proper * _B_(nu_proper, T_grid[j][k]) * fabs(phi_w[j][k]);
            // if there is contribution from the corona add it.
            if (x != 0)
              integ_0 += nu_proper * I_star(nu_proper) * fabs(phi_w[j][k]);

            switch(funcID)
            {
            // Radial component
            case F_R:
              integ_1 += Mt_grid[j][k] * mu_x[j] * I_grid[i][j][k] * fabs(phi_w[j][k]);
              break;
            // Theta component
            case F_THETA:
              integ_1 += Mt_grid[j][k] * -st1[j] * cph[j][k] * I_grid[i][j][k] * fabs(phi_w[j][k]);
              break;
            // For axisymmetric disk F_PHI will go to zero so we just ignore it
            // case F_PHI:
            //   integ_1 += Mt_grid[j][k] * st1[j] * sph[j][k] * I_grid[i][j][k] * fabs(phi_w[j][k]);
            //   break;
            }
          }
          integ_2 += integ_0 * fabs(mu_w[j]);
          integ_3 += integ_1 * fabs(mu_w[j]);
        }
        total += integ_2 * integ_3 * fabs(nu_w[i]);
      }
    }

    return 2 * wind.sigma_ref / (C * J_tot) * total;                            // double because integration is reflected over disk (axisymmetry)
}










//
