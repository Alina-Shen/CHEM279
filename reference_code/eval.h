#ifndef EVAL_H
#define EVAL_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <armadillo>

#include "Atom.h"

double E_LJ(const std::vector<Atom> &Atoms);

void F_LJ_analytical(arma::mat & F, const std::vector<Atom> &Atoms);

void F_LJ_forward_difference(arma::mat & F, const std::vector<Atom> &Atoms, double stepsize);

void F_LJ_central_difference(arma::mat & F, const std::vector<Atom> &Atoms, double stepsize);

void Steepest_descent(std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms, double fdstepsize, double searchstepsize, double fthresh);

void Steepest_descent_line_search(std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms, double fdstepsize, double searchstepsize, double fthresh, double ethresh);

#endif // EVAL_H
