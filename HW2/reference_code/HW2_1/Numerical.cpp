#include <stdio.h>
#include <math.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include "Numerical.h"

using namespace std;

double DE_product_transform(Integrand& f, Integrand& g, double t){
  double c_constant =M_PI /2.0;
  double sinht = sinh(t);
  double x = sinh(c_constant *sinht);
  double dx_dt = c_constant * cosh(t) * cosh(c_constant*sinht);
  double fxgx = f.eval(x) * g.eval(x) ;
  return fxgx*dx_dt;
}


double DEInfinity(Integrand& f1, Integrand& f2, const int numpoints_oneside, const double stepsize){
  double Sum = DE_product_transform(f1, f2, 0.0 );
  for(int k = 1; k <= numpoints_oneside; k++){
    Sum += DE_product_transform(f1, f2, k*stepsize);
    Sum += DE_product_transform(f1, f2,  -k*stepsize );
  }
  return Sum * stepsize;
}

double Trapzd(Integrand& f1, Integrand& f2, const int numpoints_oneside, const double stepsize){
  double midpoint = (f1.get_alpha()*f1.get_x0()+f2.get_alpha()*f2.get_x0())/(f1.get_alpha()+f2.get_alpha()); //Start from the center of new gaussian function
  double x;
  double Sum = f1.eval(midpoint) * f2.eval(midpoint);
  for(int k = 1; k < numpoints_oneside; k++){
    x = midpoint + k*stepsize;
    Sum += f1.eval(x) * f2.eval(x);
    x = midpoint - k*stepsize;
    Sum += f1.eval(x) * f2.eval(x);
  }
  Sum += 0.5 * (f1.eval(midpoint + numpoints_oneside*stepsize) * f2.eval(midpoint + numpoints_oneside*stepsize) 
  + f1.eval(midpoint - numpoints_oneside*stepsize) * f2.eval(midpoint - numpoints_oneside*stepsize));
  return Sum * stepsize;
}

double Gaussian::eval(double x)
{
  double xd = x - x0;  
  double y = -xd*xd * alpha;
  double l_part = pow(xd , l);

  return l_part * exp(y);
}

void Gaussian::Reset(double x0_input, double alpha_input, int l_input)
{
  x0 = x0_input;
  alpha = alpha_input;
  l = l_input;
}

void ReadGaussianparameter(Gaussian & s1, Gaussian & s2, string &fname)
{
  fstream in(fname, ios::in);
  string line1, line2;
  double x0, alpha;
  int l;
  getline(in, line1);
  getline(in, line2);
  istringstream iss1(line1);
  if (!(iss1 >> x0 >> alpha >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  s1.Reset(x0, alpha, l);
  istringstream iss2(line2);
  if (!(iss2 >> x0 >> alpha >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  s2.Reset(x0, alpha, l);
  in.close();  
}
