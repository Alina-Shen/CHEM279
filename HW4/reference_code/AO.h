#ifndef AO_H
#define AO_H

#include <iostream>
#include <armadillo>
#include <vector>

class AO
{
    private:
        arma::vec R0;
        arma::uvec lmn;
        arma::vec alpha;
        arma::vec d_coe;
        int len;
        std::string lable;
    public:
        AO(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input, arma::uvec &lmn_input, std::string lable_input);
        ~AO(){}
        void printinfo();
        arma::uvec get_lmn(){ return lmn;}
        arma::vec get_alpha(){ return alpha;}
        arma::vec get_d_coe(){ return d_coe;}
        arma::vec get_R0(){ return R0;}
        int get_len(){ return len;}
        std::string get_lable(){ return lable;}
};

class Atom{
    public:
        std::vector<AO> mAOs;
        std::string name;
        int VAN; // Valence atomic number
        Atom():name("0"){}
        Atom(std::vector<AO> AOs): mAOs(AOs), name("0"){}
        Atom(std::vector<AO> AOs, std::string atomname, int VAN_i): mAOs(AOs), name(atomname), VAN(VAN_i){}
        Atom(std::string atomname, int VAN_i): name(atomname), VAN(VAN_i){}
        void addAO(AO aAO){
            mAOs.push_back(aAO);
        }
        void PrintAOs();
};

class Molecule_basis{
    public:
        std::vector<Atom> mAtoms;
        int num_ele;
        Molecule_basis(): num_ele(0){}
        Molecule_basis(std::vector<Atom> Atoms): mAtoms(Atoms), num_ele(0){}
        Molecule_basis(std::vector<Atom> Atoms, int cha): mAtoms(Atoms), num_ele(cha){}
        Molecule_basis(std::string &fname, arma::mat &H_basis, arma::mat &C_basis);
        void addAtom(Atom aAtom){
            mAtoms.push_back(aAtom);
        }
        void setnum_ele(int cha){num_ele = cha;}
        void PrintAtoms();
        int getnum_ele(){return num_ele;}
        int getnum_atoms(){return mAtoms.size();}
        int getnumAOs();
        void eval_OVmat(arma::mat &OV_mat);
        // void eval_Hmat(arma::mat &OV_mat, arma::mat &H_mat);
        void eval_gammamat(arma::mat &gamma_mat);
};

double Eval_Ov_AOs(AO& sh1, AO& sh2);

double Eval_2eI_sAO(AO& sh1, AO& sh2);


#endif // AO_H
