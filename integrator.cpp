#include "integrator.hpp"

#define ldouble long double

/*
   From: "COMPUTATIONAL PHYSICS, 2nd Ed"
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Wiley-VCH, 2007.
   Electronic Materials copyright: R Landau, Oregon State Univ, 2007;
   MJ Paez, Univ Antioquia, 2007; & CC Bordeianu, Univ Bucharest, 2007
   Support by National Science Foundation
*/
void Integrator::gauss(int n, int job, ldouble a, ldouble b, ldouble* x, ldouble* w)
{

//     n     number of points
//     job = 0  rescaling uniformly between (a,b)
//           1  for integral (0,b) with 50% pts inside (0, ab/(a+b)
//           2  for integral (a,inf) with 50% inside (a,b+2a)
//     x, w     output grid points and weights.

    int     m, i, j;
    ldouble  t, t1, pp, p1, p2, p3;
    ldouble  eps = 3.e-16;
  	                                         // eps = accuracy to adjust
    m = (n + 1) / 2;
    for (i = 1; i <= m; ++i) {
        t = cos(PI * (i - 0.25) / (n + 0.5));
        t1 = 1;
        while((fabs(t - t1)) >= eps) {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; ++j) {
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) * t * p2 - (j-1) * p3) / j;
            }
            pp = n * (t * p1 - p2) / (t*t - 1);
            t1 = t;
            t = t1 - p1 / pp;
        }
        x[i-1] = -t;
        x[n-i] = t;
        w[i-1]   = 2.0 / ((1 - t*t) * pp*pp);
        w[n-i] = w[i-1];
    }
    if (job == 0)   {
        for (i = 0; i < n ; ++i)  {
            x[i] = x[i] * (b - a) / 2.0 + (b + a) / 2.0;
            w[i] = w[i] * (b - a) / 2.0;
        }
    }
    if (job == 1)  {
        for (i = 0; i < n; ++i)  {
          x[i] = a * b * (1 + x[i]) / (b + a - (b - a) * x[i]);
          w[i] = w[i] * 2 * a * b * b /((b + a - (b - a) * x[i]) * (b + a - (b - a) * x[i]));
       }
    }
    if (job == 2)   {
        for (i = 0; i < n; ++i)  {
            x[i] = (b * x[i] + b + a + a) / (1 - x[i]);
            w[i] =  w[i] * 2 * (a + b)  / ((1 - x[i]) * (1 - x[i]));
        }
    }
}

void Integrator::gauss(int n, int job, ldouble a, ldouble b, std::vector<ldouble>& x, std::vector<ldouble>& w)
{
  ldouble* x_ = new ldouble[n];
  ldouble* w_ = new ldouble[n];
  gauss(n, job, a, b, x_, w_);

  x.clear();
  w.clear();
  for(int i = 0; i < n; ++i)
  {
    x.push_back(x_[i]);
    w.push_back(w_[i]);
  }

  delete[] x_;
  delete[] w_;
}

void Integrator::gaulag(int n, ldouble alpha, ldouble* x, ldouble* w)
{
	// The following code is full of magic that I cannot fully explain
  // it is from Numerical Recipes on pages 184 and 185

	int MAXIT = 10;
	ldouble eps = 1.0e-14;

  int ai;
  ldouble z, z1, p1, p2, p3, pp;

	for(int i = 0; i < n; ++i)
  {
		// These create close and stable guesses for the nth root of the general laguerre polynomial L_n^alpha(x)
		if(i == 0)
    {
			z = (1 + alpha) * (3.0 + 0.92 * alpha) / (1 + 2.4 * n + 1.8 * alpha);
    }
    else if(i == 1)
    {
			z += (15.0 + 6.25*alpha) / (1.0 + 0.9*alpha + 2.5*n);
		}
    else
    {
			ai = i - 1;
			z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alpha / (1.0 + 3.5 * ai)) * (z - x[i-2]) / (1.0 + 0.3 * alpha);
    }

		// And this is the root finder
		for(int its = 0; its < MAXIT; ++its)
    {
			p1 = 1;
			p2 = 0;
			for(int j = 0; j < n; ++j)
      {
				p3 = p2;
				p2 = p1;
				p1 = ((2 * j + 1 + alpha - z) * p2 - (j + alpha) * p3) / (j + 1);
      }

			pp = (n * p1 - (n + alpha) * p2) / z;
			z1 = z;
			z = z1 - p1/pp;
			if(fabs(z - z1) <= eps) { break; }
    }

		x[i] = z; 														                        // Abscissa x_i is the ith root of L_n^alpha(x)
		w[i] = -exp(lgamma(alpha + n) - lgamma(n)) / (pp * n * p2);	// Weight at abscissa x_i
  }
}

void Integrator::gaulag(int n, ldouble alpha, std::vector<ldouble>& x, std::vector<ldouble>& w)
{
  ldouble* x_ = new ldouble[n];
  ldouble* w_ = new ldouble[n];
  gaulag(n, alpha, x_, w_);

  x.clear();
  w.clear();
  for(int i = 0; i < n; ++i)
  {
    x.push_back(x_[i]);
    w.push_back(w_[i]);
  }

  delete[] x_;
  delete[] w_;
}
