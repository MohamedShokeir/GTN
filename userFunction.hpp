#ifndef _USER_FUNCTION_HPP_
#define _USER_FUNCTION_HPP_

#include <vector>
#include <math.h>

inline double myfunction(double& variable, const std::vector<double> parameters){
    
  double seq   = parameters[0];
  double tr    = parameters[1];
  double q1    = parameters[2];
  double q2    = parameters[3];
  double fstar = parameters[4];

  //return ((seq * seq) / (variable * variable)) + (2. * q1 * fstar * cosh(q2 * tr * 0.5 / variable)) - (1.) - ((q1 * q1) * (fstar * fstar));
  return (seq * seq) + (variable * variable) * ((2. * q1 * fstar * cosh(q2 * tr * 0.5 / variable)) - (1.) - ((q1 * q1) * (fstar * fstar)));
  //return ((seq * seq) * (variable * variable)) + (2. * q1 * fstar * cosh(q2 * tr * 0.5 * variable)) - (1.) - ((q1 * q1) * (fstar * fstar));
}

inline double derivMyFunction(double& variable, const std::vector<double> parameters){

//  double seq   = parameters[0];
  double tr    = parameters[1];
  double q1    = parameters[2];
  double q2    = parameters[3];
  double fstar = parameters[4];
  
  //return (-2. * (seq * seq) / (variable * variable * variable)) - (q1 * q2 * fstar / (variable * variable) * tr * sinh(q2 * tr * 0.5 / variable));
  return (2. * variable) * ((2. * q1 * fstar * cosh(q2 * tr * 0.5 / variable)) - (1.) - ((q1 * q1) * (fstar * fstar))) - (q1 * q2 * fstar * tr * sinh(q2 * tr * 0.5 / variable));
  //return (-2. * (seq * seq) * (variable * variable * variable)) - (q1 * q2 * fstar * (variable * variable) * tr * sinh(q2 * tr * 0.5 * variable));
}

#endif
