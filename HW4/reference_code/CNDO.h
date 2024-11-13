#if !defined CNDO_H
#define CNDO_H
#include <armadillo>
#include <cassert>
#include "AO.h"


class SCF{
    public:
        virtual int init()=0;
        virtual int run()=0;
        virtual arma::mat getPa()=0;
        virtual arma::mat getPb()=0;
        virtual double getEnergy()=0;
};

class CNDO: public SCF{
    private:
        Molecule_basis &mol;
        int dim, max_iter;
        int p, q; //num of alpha and beta electrons
        double tol;
        arma::mat S;
        arma::mat gamma;
        arma::mat H_core;
        arma::mat Pa, Pb;
        arma::mat Ga, Gb;
        arma::mat Ca, Cb;
        arma::vec Ea, Eb;
        double Ee, Ec, Etotal;
    public:
        CNDO(Molecule_basis &mol_i, int max_it, double tolerence);
        virtual int init();
        int updateG();
        virtual int run();
        virtual arma::mat getPa() {return Pa;}
        virtual arma::mat getPb() {return Pb;}
        virtual double getEnergy();
};



double Solve_CNDO(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat, arma::vec &energy_vec, int num_ele);




#endif // CNDO_H
