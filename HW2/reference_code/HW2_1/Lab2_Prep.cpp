// Lab2 Prep.

/* Structure of HW submittion on github
    5 HW in the same repo
    HW# directory should includes:
        code folder
        sample input folder
        sample output folder
        output folder
        README.md
        Discussion file
        makefile
*/

/* Report submitted onto github
    link of your github repo
    instruction of how to use your code (README.md)
    discussion about your code and results (README.md)
*/

/*
    won't be harsh on grading! no worries
    main point is to see if you understand the concepts and algorithms
    and be able to implement it into code
    other things like formats/error reports won't take lots of points
*/

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <armadillo>

// Q1

// Base class for integrand
// An abstract base class cannot be directly instantiated
class Integrand
{
    public:
    virtual double eval(double x) = 0;
    virtual double get_x0() const = 0;
    virtual double get_alpha() const = 0;
    virtual int get_l() const = 0;
};


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

void ReadGaussianpara(Gaussian & s1, Gaussian & s2, string &fname)
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
}


double Trapzd(Integrand& f1, Integrand& f2,const int numpoints_oneside, const double stepsize)
{
    ...
}

// how to choose integral range? 95%?


// Q2

//The class for shell of primitive GTOs
class Shell
{
    private:
    arma::vec R0; //center of the shell
    double alpha; //exponents
    int l; // total angular momentum
    
    public:
    Shell(double x0_input, double y0_input, double z0_input, double alpha_input, int l_input):
    alpha(alpha_input), l(l_input) {R0={x0_input, y0_input, z0_input};}//{R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input;}
    Shell():alpha(0.5), l(0) {R0.zeros(3);}
    ~Shell(){}
    void Reset(double x0_input, double y0_input, double z0_input, double alpha_input, int l_input);
    void printinfo();
    // give the numbers of functions in this shell
    int dim_func();
    int get_l(){ return l;}
    double get_alpha(){ return alpha;}
    arma::vec get_R0(){ return R0;}
    double eval(double x);
};

void ReadShellparameter(Shell & sh1, Shell & sh2, std::string &fname);
double OverlapOneDim(double xa, double xb, double alphaa, double alphab, int la, int lb) // 1-D analytical overlap integral
{
    // prefator
    // combination
    // double factorial
    // numerator
    // denominator
    // how many nested loops?
}

void EvalOverlap(arma::mat &Overlap, Shell& sh1, Shell& sh2) // Evaluate the whole overlap matrix between two shells
{
    // use OverlabOneDim
    // how many AOs?
    // how many nested loops?
    // Overlap of (la_x, la_y, la_z) & (lb_x, lb_y, lb_z) -> S(a,b)
}
