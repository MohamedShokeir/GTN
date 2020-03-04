#ifndef _NEWTON_HPP_
#define _NEWTON_HPP_

#include <functional>
#include <vector>
#include <math.h>
#include "userFunction.hpp"
using namespace std; 

inline double scalarNewtonSolver(std::function<double (double&, const std::vector<double>&)> f, std::function<double (double&, const std::vector<double>&)> df, const std::vector<double>& param, const double& x1){

  double epsilon = 1.e-6;
  int niter_max  = 100;
  int niter      = 0;
  double x       = x1;
  double resinit = f(x, param);
  double res     = resinit;
  //cout << "resinit= " << resinit << endl;
  do{
    x  -= f(x, param) / df(x, param);
    res = f(x, param);
    niter++;
    //cout << "res= " << res << endl;
    //cout << "niter= " << niter << endl;
    //cout << "---------------------------" << endl;
  }while((fabs(res/resinit) > epsilon) && (niter < niter_max));
  
  return x;
}

#endif
