#include "AO.h"
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <cassert>
#include <string>
#include "util.h"
#include <vector>

using namespace std;

//AO functions
void AO::printinfo(){
  printf("This AO info: %s, R( %1.2f, %1.2f, %1.2f), with angular momentum: %lld %lld %lld\n", lable.c_str(),
    R0(0), R0(1), R0(2), lmn(0), lmn(1), lmn(2));
  d_coe.print("d_coe");
  alpha.print("alpha");
}

AO::AO(arma::vec &R0_input, arma::vec &alpha_input,  arma::vec &d_input, arma::uvec &lmn_input, string lable_input):
R0(R0_input), alpha(alpha_input), d_coe(d_input), lmn(lmn_input), lable(lable_input){
    assert(R0.n_elem == 3);
    assert(lmn.n_elem == 3);
    len = alpha.n_elem;
    assert(d_coe.n_elem == len);
    for (size_t k = 0; k <len; k++){
      double Overlap_Self = Overlap_3d(R0, R0, alpha(k), alpha(k), lmn, lmn);
      d_coe(k) /= sqrt(Overlap_Self);
    }
}


std::map<std::string, int> VAN_map { {"H", 1}, {"C", 4}, {"N", 5}, {"O", 6}, {"F", 7}, };
Atom GenerateAtom(std::string atomname, arma::vec R0){
  std::string basisname = atomname + "_STO3G.txt";
  arma::mat basis;
  basis.load(basisname);
  arma::uvec lmn = {0, 0, 0};
  arma::vec alpha = basis.col(0);
  arma::vec d_coe = basis.col(1);
  string lable("s");
  AO AO_s(R0, alpha, d_coe, lmn, lable);
  if (atomname == string("H")){
    Atom aatom(atomname, 1);
    aatom.addAO(AO_s);
    return aatom;
  }
  if(VAN_map.find(atomname) == VAN_map.end()){
    throw invalid_argument("Do not support this kind of atom.");
  }
  int atomicnumber = VAN_map[atomname];
  Atom aatom(atomname, atomicnumber);
  aatom.addAO(AO_s);
  for(size_t j = 0; j < 3; j++){
    d_coe = basis.col(2);
    lmn.zeros();
    lmn(j) = 1;
    string lable("p");
    AO readedAOp(R0, alpha, d_coe, lmn, lable);
    aatom.addAO(readedAOp);
  }
  return aatom;
}

double Eval_Ov_AOs(AO &ao1, AO &ao2)
{
  int len = ao1.get_len();
  assert(ao2.get_len() == len);
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec Ra = ao1.get_R0(), Rb = ao2.get_R0();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();

  double sum = 0.;
  for (size_t k = 0; k < len; k++)
  {
    double alpha_k = alphaa(k);
    for (size_t j = 0; j < len; j++)
    {
      double alpha_j = alphab(j);
      double Overlap = Overlap_3d(Ra, Rb, alpha_k, alpha_j, la, lb);
      // printf("%ld %ld = %1.10f\n", k, j, Overlap);
      sum += da(k) * db(j) * Overlap;
    }
  }
  return sum;
}

double Eval_2eI_sAO(AO& ao1, AO& ao2){
  // ao1.printinfo();
  // ao2.printinfo();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  if (!(arma::accu(la) == 0 && arma::accu(lb) == 0 ))
    throw invalid_argument("Eval_2eI_sAO is only used for s Orbitals.");
  int len = ao1.get_len();
  assert(ao2.get_len() == len);
  arma::vec Ra = ao1.get_R0(), Rb = ao2.get_R0();
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();
  
  double sum = 0.;
  for (size_t k1 = 0; k1 < len; k1++)
  for (size_t k2 = 0; k2 < len; k2++)
  {
    double sigmaA = 1.0/ (alphaa(k1)+ alphaa(k2));
    for (size_t j1 = 0; j1 < len; j1++)
    for (size_t j2 = 0; j2 < len; j2++)
    {
      double sigmaB = 1.0/ (alphab(j1)+ alphab(j2));
      double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB);
      // printf("%ld %ld %ld %ld = %1.10f\n", k1, k2, j1, j2, I2e);
      // return I2e;
      sum += da(k1) * da(k2) * db(j1) *db(j2) * I2e;
    }
  }
  return sum;    
}

std::map<std::string, std::string> AN_map { {"1", "H"}, {"6", "C"}, {"7", "N"}, {"8", "O"}, {"9", "F"} };
Molecule_basis::Molecule_basis(string &fname, arma::mat &H_basis, arma::mat &C_basis)
{
  int basislen = H_basis.n_rows;
  assert(C_basis.n_rows == basislen);
  int num_charge, num_Atoms;

  ifstream in(fname, ios::in);
  // cout << fname;

  string line;
  getline(in, line);
  istringstream iss(line);
  if (!(iss >> num_Atoms >> num_charge))
    throw invalid_argument("There is some problem with molecule format.");
  int count_atoms = 0;

  while (getline(in, line))
  {
    istringstream iss(line);
    arma::vec R0(3);
    // int AtomicN = 0;
    string atomnumber;
    if (!(iss >> atomnumber >> R0[0] >> R0[1] >> R0[2]))
      throw invalid_argument("There is some problem with AO format.");
    
    if(AN_map.find(atomnumber) == AN_map.end()){
        throw invalid_argument("Do not support this kind of atom.");
    }

    string atomname = AN_map[atomnumber];

    arma::uvec lmn = {0, 0, 0};
    arma::vec alpha(basislen);
    arma::vec d_coe(basislen);
    Atom readAtom =GenerateAtom(atomname, R0);
    mAtoms.push_back(readAtom);
    // cout << readAO << std::endl;
    count_atoms ++;
  }
  if(count_atoms != num_Atoms){
    throw invalid_argument("Number of AOs are not consistent ");
  }
  in.close();
  num_ele = 0;
  for(auto atom : mAtoms)
    num_ele += atom.VAN;
  num_ele -= num_charge;
}

void Atom::PrintAOs(){
  printf("THis is a %s atom\n", name.c_str());
  for(auto ao : mAOs)
    ao.printinfo();
  printf("\n");
}
void Molecule_basis::PrintAtoms(){
  for(auto atom : mAtoms)
    atom.PrintAOs();
}
int Molecule_basis::getnumAOs(){
  int numAOs = 0;
  for(auto atom : mAtoms)
    numAOs += atom.mAOs.size();
  return numAOs;
}


void Molecule_basis::eval_OVmat(arma::mat &OV_mat){
  int dim = getnumAOs();
  assert(OV_mat.n_rows == dim && OV_mat.n_cols == dim);
  vector<AO> allAOs;
  for(auto atom : mAtoms)
    for(auto ao: atom.mAOs)
      allAOs.push_back(ao);
  for (size_t k = 0; k <dim; k++){
    for (size_t j = 0; j <= k; j++){
      double OV_1AO = Eval_Ov_AOs(allAOs[k], allAOs[j]);
      OV_mat(k,j) = OV_1AO;
      OV_mat(j,k) = OV_1AO;
    }
  }
}

void Molecule_basis::eval_gammamat(arma::mat &gamma_mat){
  int dim = mAtoms.size();
  assert(gamma_mat.n_rows == dim && gamma_mat.n_cols == dim);
  for (size_t k = 0; k <dim; k++){
    for (size_t j = 0; j <= k; j++){
      double OV_1AO = Eval_2eI_sAO(mAtoms[k].mAOs[0], mAtoms[j].mAOs[0]);
      gamma_mat(k,j) = OV_1AO;
      gamma_mat(j,k) = OV_1AO;
    }
  }
}











