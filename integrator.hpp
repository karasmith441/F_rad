#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <math.h>
#include <vector>

#include "consts.hpp"

#define ldouble long double

class Integrator
{

public:

  Integrator() {}

  void gauss(int, int, ldouble, ldouble, ldouble*, ldouble*);
  void gauss(int, int, ldouble, ldouble, std::vector<ldouble>&, std::vector<ldouble>&);

  void gaulag(int, ldouble, ldouble*, ldouble*);
  void gaulag(int, ldouble, std::vector<ldouble>&, std::vector<ldouble>&);

private:

};

#endif
