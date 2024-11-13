#include "CNDO.h"
#include <stdio.h>
#include <math.h>
#include <cassert>
#include <unordered_map>
#include <armadillo>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

double hartree_to_ev = 27.211396641308;
struct CNDO_para
{
  std::map<std::string, double> IA;
  double beta;
} C_para = {{{"s", 14.051}, {"p", 5.572}}, 21.}, H_para{{{"s", 7.176}}, 9.};
CNDO_para N_para = {{{"s", 19.316}, {"p", 7.275}}, 25.};
CNDO_para O_para = {{{"s", 25.390}, {"p", 9.111}}, 31.};
CNDO_para F_para = {{{"s", 32.272}, {"p", 11.080}}, 39.};

std::map<std::string, CNDO_para> CNDO_para_map{
    {"H", H_para},
    {"C", C_para},
    {"O",O_para},
    {"N",N_para},
    {"F", F_para},
};

CNDO::CNDO(Molecule_basis &mol_i, int max_it, double tolerence) : mol(mol_i), max_iter(max_it), tol(tolerence)
{
  dim = mol.getnumAOs();
  Pa = arma::zeros(dim, dim);
  Pb = arma::zeros(dim, dim);
  Ga = arma::zeros(dim, dim);
  Gb = arma::zeros(dim, dim);
  Ca.set_size(dim, dim);
  Cb.set_size(dim, dim);
  Ea.set_size(dim);
  Eb.set_size(dim);

  H_core.set_size(dim, dim);
  S.set_size(dim, dim);
  mol.eval_OVmat(S);
  int num_atoms = mol.getnum_atoms();
  gamma.set_size(num_atoms, num_atoms);
  mol.eval_gammamat(gamma);
  gamma *= hartree_to_ev;
  // std::cout << std::setprecision(3);
  gamma.print("gamma");
  S.print("Overlap");
  // for (size_t k = 0; k < num_atoms; k++)
  //   gamma(k,k) *= 0.9365244;
  q = mol.num_ele / 2;
  p = mol.num_ele - q;
  cout << " p =  " << p << " q =  " << q  << endl;
  Ee = 0.;
  Etotal = 0.;
  Ec = 0.;
  for (size_t k = 0; k < num_atoms; k++)
    for (size_t j = 0; j < k; j++)
    {
      arma::vec Ra = mol.mAtoms[k].mAOs[0].get_R0(), Rb = mol.mAtoms[j].mAOs[0].get_R0();
      double Rd = arma::norm(Ra - Rb, 2);
      Ec += mol.mAtoms[k].VAN * mol.mAtoms[j].VAN / Rd;
    }
  Ec *= hartree_to_ev;
}

int CNDO::init()
{
  int num_atoms = mol.getnum_atoms();
  //Do H_core part
  size_t k_AO = 0;
  for (size_t k = 0; k < num_atoms; k++)
  {
    Atom &A_AO = mol.mAtoms[k];
    double ZA = double(A_AO.VAN);
    double gammaAA = gamma(k, k);
    CNDO_para A_para = CNDO_para_map[A_AO.name];
    for (auto ao_A : A_AO.mAOs)
    {
      H_core(k_AO, k_AO) = -A_para.IA[ao_A.get_lable()] - (ZA - 0.5) * gammaAA;
      // H_core.print("H_core");
      size_t j_AO = 0;
      for (size_t j = 0; j < num_atoms; j++)
      {
        Atom &B_AO = mol.mAtoms[j];
        CNDO_para B_para = CNDO_para_map[B_AO.name];
        double averagebeta = -(A_para.beta + B_para.beta) / 2.;
        if (k != j)
          H_core(k_AO, k_AO) -= double(B_AO.VAN) * gamma(k, j);
        // cout<< k_AO << " " << H_core(k_AO, k_AO)<< endl;
        for (auto ao : B_AO.mAOs)
        {
          if (k_AO != j_AO)
            H_core(k_AO, j_AO) = averagebeta * S(k_AO, j_AO);
          j_AO++;
        }
      }
      k_AO++;
    }
  }
  if (k_AO != dim)
  {
    cout << "warn! the number of AOs is wrong." << endl;
    return 1;
  }
  H_core.print("H_core");

  return 0;
}

int CNDO::updateG()
{
  // std::cout << std::setprecision(3);
  int num_atoms = mol.getnum_atoms();
  arma::vec P_t = arma::zeros(num_atoms);
  size_t k_AO = 0;
  for (size_t k = 0; k < num_atoms; k++)
    for (auto ao : mol.mAtoms[k].mAOs)
    {
      P_t(k) += Pa(k_AO, k_AO) + Pb(k_AO, k_AO);
      k_AO++;
    }
  P_t.print("P_t");
  k_AO = 0;
  for (size_t k = 0; k < num_atoms; k++)
  {
    double gammaAA = gamma(k, k);
    for (auto ao_A : mol.mAtoms[k].mAOs)
    {
      Ga(k_AO, k_AO) = (P_t(k) - Pa(k_AO, k_AO)) * gammaAA;
      Gb(k_AO, k_AO) = (P_t(k) - Pb(k_AO, k_AO)) * gammaAA;
      size_t j_AO = 0;
      for (size_t j = 0; j < num_atoms; j++)
      {
        double gammaAB = gamma(k, j);
        if (k != j)
        {
          Ga(k_AO, k_AO) += P_t(j) * gammaAB;
          Gb(k_AO, k_AO) += P_t(j) * gammaAB;
        }
        for (auto ao : mol.mAtoms[j].mAOs)
        {
          if (k_AO != j_AO)
          {
            Ga(k_AO, j_AO) = -gammaAB * Pa(k_AO, j_AO);
            Gb(k_AO, j_AO) = -gammaAB * Pb(k_AO, j_AO);
            // cout << k_AO << "  " << j_AO<< "  " << gammaAB<<"  " << Pa(k_AO, j_AO)<< "  " << Ga(k_AO, j_AO)<< endl;
          }
          j_AO++;
        }
      }
      k_AO++;
    }
  }
  Ga.print("Ga");
  Gb.print("Gb");
  return 0;
}

int CNDO::run()
{
  arma::mat Fa = H_core + Ga;
  arma::mat Fb = H_core + Gb;
  arma::mat Pa_old, Pb_old;
  size_t k = 0;
  for (; k < max_iter; k++)
  {
    cout << "Iteration: " << k << endl;
    // Pa.print("Pa");
    // Ga.print("Ga");
    Fa.print("Fa");
    Fb.print("Fb");
    Pa_old = Pa;
    Pb_old = Pb;
    arma::eig_sym(Ea, Ca, Fa);
    arma::eig_sym(Eb, Cb, Fb);
    cout << "after solving eigen equation: " << k << endl;
    Ca.print("Ca");
    Cb.print("Cb");
    cout << " p =  " << p << " q =  " << q  << endl;
    Pa = Ca.cols(0, p - 1) * Ca.cols(0, p - 1).t();
    if (q > 0)
      Pb = Cb.cols(0, q - 1) * Cb.cols(0, q - 1).t();
    else
      Pb.zeros();
    if (arma::approx_equal(Pa, Pa_old, "absdiff", tol) && arma::approx_equal(Pb, Pb_old, "absdiff", tol))
      break;
    Pa.print("Pa_new");
    Pb.print("Pb_new");
    updateG();
    Fa = H_core + Ga;
    Fb = H_core + Gb;
  }
  if (k == max_iter)
  {
    cout << "Error: the job could not be finished in " << max_iter << "iterations.\n";
    return 1;
  }
  Ea.print("Ea");
  Eb.print("Eb");
  Ca.print("Ca");
  Cb.print("Cb");
  return 0;
}

double CNDO::getEnergy()
{
  arma::mat Ptotal = Pa + Pb;
  Ee = arma::dot(Pa, Ga) / 2. + arma::dot(Pb, Gb) / 2.;
  Ee += arma::dot(Ptotal, H_core);
  Etotal = Ee + Ec;
  cout << "Nuclear Repulsion Energy is " << Ec << " eV." << endl;
  cout << "Electron Energy is " << Ee << " eV." << endl;
  return Etotal;
}

