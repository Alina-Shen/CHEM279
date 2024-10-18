#ifndef NUMERICAL_H
#define NUMERICAL_H

#ifndef M_PI
#define M_PI 3.14159265
#endif

#include <iostream>

//Base class for integrand
class Integrand
{
public:
virtual double eval(double x) = 0;
virtual double get_x0() const = 0;
virtual double get_alpha() const = 0;
virtual int get_l() const = 0;
};

//DE rule for numerical integration, it should converge faster than trapezoidal rule
double DEInfinity(Integrand& f1, Integrand& f2,const int numpoints_oneside, const double stepsize);

//Trapezoidal rule
double Trapzd(Integrand& f1, Integrand& f2,const int numpoints_oneside, const double stepsize);

//Derived class for 1-D gaussian type orbitals
class Gaussian: public Integrand
{
    private:
    double x0;//center
    double alpha; //exponent
    int l; //1-D angular momentum
    public:
    Gaussian(double x0_input, double alpha_input, int l_input):
    x0(x0_input), alpha(alpha_input), l(l_input) {}
    Gaussian():x0(0.0), alpha(0.5), l(0){}
    ~Gaussian(){}
    void Reset(double x0_input, double alpha_input, int l_input);
    double eval(double x);
    double get_x0() const {return x0;}
    double get_alpha() const {return alpha;}
    int get_l() const {return l;}
};

void ReadGaussianparameter(Gaussian & s1, Gaussian & s2, std::string &fname);

#endif // NUMERICAL_H
