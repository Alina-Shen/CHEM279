#ifndef SHELL_H
#define SHELL_H

#ifndef M_PI
#define M_PI 3.14159265
#endif

#include <iostream>
#include <armadillo>

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
double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb); // 1-D analytical overlap integral

void Eval_Ov(arma::mat &Overlap, Shell& sh1, Shell& sh2); // Evaluate the whole overlap matrix between two shells

#endif // SHELL_H
