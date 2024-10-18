#include "eval.h"
#include <vector>
#include <math.h>
#include <armadillo>

using namespace std;

const double r_Au = 2.951;
const double e_Au = 5.29;

// Initializes and populates a matrix with 3D coordinates from a vector of Atom objects.
void Convert_atoms_to_corr(arma::mat & coor, const vector<Atom> & Atoms) {
    coor.zeros(3,Atoms.size());
    for (int i = 0; i < Atoms.size(); i++) {
        for (int k = 0; k < 3; k++) {
            coor(k,i) = Atoms[i].r[k];
        }
    }
}

// Updates the 3D coordinates of Atom objects in a vector using values from a matrix.
void Convert_corr_to_atoms(vector<Atom> & Atoms, const arma::mat & coor) {
    for (int i = 0; i < Atoms.size(); i++) {
        for (int k = 0; k < 3; k++) {
            Atoms[i].r[k] = coor(k,i);
        }
    }
}

//Calculate the energy of the system with LJ potential
double E_LJ(const std::vector<Atom> &Atoms) {
    double E = 0;
    for (int i = 0; i < Atoms.size(); i++)
    {
        for (int j = i + 1; j < Atoms.size(); j++)
        {
            double R2_ij = (Atoms[i].r[0] - Atoms[j].r[0])*(Atoms[i].r[0] - Atoms[j].r[0])
            + (Atoms[i].r[1] - Atoms[j].r[1])*(Atoms[i].r[1] - Atoms[j].r[1])
            + (Atoms[i].r[2] - Atoms[j].r[2])*(Atoms[i].r[2] - Atoms[j].r[2]);
            double r2overR2 = (r_Au * r_Au) / R2_ij;
            E += e_Au * (pow(r2overR2, 6) - 2*pow(r2overR2, 3));
        }
    }
    return E;
}

//Calculate LJ pontential for the matrix of 3D coordinates
double E_LJ(arma::mat coor) {
    int n_atoms = coor.n_cols;
    std::vector<Atom> Atoms(n_atoms);
    Convert_corr_to_atoms(Atoms, coor);
    return E_LJ(Atoms);
}

//Analytical force for LJ potential
void F_LJ_analytical(arma::mat & F, const vector<Atom> & Atoms) {
    F.zeros(3,Atoms.size()); 
    for (int i = 0; i < Atoms.size(); i++){
        for (int j = 0; j < Atoms.size(); j++){
            if (i != j) {
                double R2_ij = (Atoms[i].r[0] - Atoms[j].r[0])*(Atoms[i].r[0] - Atoms[j].r[0])
                + (Atoms[i].r[1] - Atoms[j].r[1])*(Atoms[i].r[1] - Atoms[j].r[1])
                + (Atoms[i].r[2] - Atoms[j].r[2])*(Atoms[i].r[2] - Atoms[j].r[2]);
                for (int k = 0; k < 3; k++) {
                    F(k,i) += 12 * e_Au * (pow(r_Au, 12) / pow(R2_ij, 7) - pow(r_Au, 6) / pow(R2_ij, 4)) * (Atoms[i].r[k] - Atoms[j].r[k]);
                }
            }
        }
    }
}

//Calculate force using forward difference
void F_LJ_forward_difference(arma::mat & F, const vector<Atom> & Atoms, double stepsize) {
    F.zeros(3,Atoms.size());
    double E_0 = E_LJ(Atoms);
    for (int i = 0; i < Atoms.size(); i++){
        for (int k = 0; k < 3; k++) {
            vector<Atom> Atoms_forward = Atoms;
            Atoms_forward[i].r[k] += stepsize;
            double E_forward = E_LJ(Atoms_forward);
            F(k,i) = - (E_forward - E_0)/stepsize;
        }
    }
}

//Calculate force using central difference
void F_LJ_central_difference(arma::mat & F, const vector<Atom> & Atoms, double stepsize) {
    F.zeros(3,Atoms.size());
    for (int i = 0; i < Atoms.size(); i++){
        for (int k = 0; k < 3; k++) {
            vector<Atom> Atoms_forward = Atoms;
            Atoms_forward[i].r[k] += stepsize;
            double E_forward = E_LJ(Atoms_forward);
            vector<Atom> Atoms_backward = Atoms;
            Atoms_backward[i].r[k] -= stepsize;
            double E_backward = E_LJ(Atoms_backward);
            F(k,i) = - (E_forward - E_backward)/(2*stepsize);
        }
    }
}

//Standard steepest descent
void Steepest_descent(vector<Atom> &opt_Atoms, const vector<Atom> &Atoms, double fdstepsize, double searchstepsize, double fthresh) {
    arma::mat old_coor, new_coor;
    double E_LJ_before = E_LJ(Atoms);
    Convert_atoms_to_corr(old_coor, Atoms);
   
    cout << "Initial energy: "<< E_LJ_before << endl;
    cout << "Stepsize for central difference is:" << fdstepsize << ";Initial stepsize for steepest descent is:" << searchstepsize << ";Threshold for convergence in force is:" << fthresh << endl; 
    arma::mat aForce, fForce, cForce;
    F_LJ_analytical(aForce, Atoms);
    aForce.print("Analytical Force");
    F_LJ_forward_difference(fForce, Atoms, fdstepsize);
    fForce.print("Forward Difference Force");
    F_LJ_central_difference(cForce, Atoms, fdstepsize);
    cForce.print("Central Difference Force");
    int count = 0;
	
    cout << "Start steepest descent with central difference force." << endl;
    
    //The main loop, I used the Frobenius norm of force as the criterion for convergence
    while(arma::norm(cForce,"fro") > fthresh) {
        new_coor = old_coor + searchstepsize * cForce/arma::norm(cForce,"fro"); //find a new point, arma::norm calculates the norm of a matrix, which is used to get a unit vector along the search direction
        double E_LJ_after = E_LJ(new_coor);
        if (E_LJ_after < E_LJ_before) {
            old_coor = new_coor;
            E_LJ_before = E_LJ_after;
            Convert_corr_to_atoms(opt_Atoms, new_coor);
            F_LJ_analytical(cForce, opt_Atoms); //calculate the force at the new point
            searchstepsize *= 1.2; //increase the stepsize if a new point with lower energy is found
        }
        else {
            searchstepsize /= 2; //Halve the step size if a point with lower energy is not found
        }
        count++;
    }
    
    cout << "Total iterations: " << count << endl;
    printf("Final energy: %.3e\n", E_LJ_before);
}

//Steepest descent with 1d line search
void Steepest_descent_line_search(vector<Atom> &opt_Atoms, const vector<Atom> &Atoms, double fdstepsize, double searchstepsize, double fthresh, double ethresh) {
    double golden_ratio = 0.38197;
    // first bracket A<B<D (A is the inital point), such that E(B)<E(A),E(B)<E(D)
    // during line search, we keep the order A<B<C<D
    arma::mat A_point(3, Atoms.size());
    Convert_atoms_to_corr(A_point, Atoms);
    arma::mat B_point(3, Atoms.size());
    arma::mat C_point(3, Atoms.size());
    arma::mat D_point(3, Atoms.size());
    // length between points
    double AB, BD, AD;
    double golden_tol = 1e-7; // tolerance to stop golden search
    double init_searchstepsize = searchstepsize;
    
    double E_LJ_A = E_LJ(Atoms);
    cout << "Initial energy: "<< E_LJ_A << endl;
    arma::mat cForce, unitForce;
    cout << "Stepsize for central difference is:" << fdstepsize << ";Initial stepsize for line search is:" << searchstepsize << ";Threshold for convergence in force is:" << fthresh << endl;
    F_LJ_central_difference(cForce, Atoms, fdstepsize);
    cForce.print("Central Difference Force");

    double old_e = 1e308;
    double cur_e = E_LJ_A;

    int count = 0;
    cout << "Start steepest descent with golden section line search using central difference force" << endl;
    //main loop
    while(arma::norm(cForce,"fro") > fthresh || abs(cur_e - old_e) > ethresh) { //Use both norm of Force and different of energy as the criterion to better confirm convergence
        unitForce = cForce/arma::norm(cForce,"fro"); //Use unit force (unit vector) in case the force is too large or small
        count++;
	//search for point B
        int bracket_count = 0;
        bool FindB = true;
        B_point = A_point + searchstepsize * unitForce;
        double E_LJ_B = E_LJ(B_point);
        while (E_LJ_B > E_LJ_A) {
            searchstepsize /= 2; //Halve the step size if B is not found in this iteration
            B_point = A_point + searchstepsize * unitForce;
            E_LJ_B = E_LJ(B_point);
            bracket_count++;
	    //break the loop if we cannot find B in finite iterations
            if (bracket_count > 100) {
                cout << "Cannot find point B" << endl;
                FindB = false;
                break;
            }
        }
        //search for point D if point B is found
        bracket_count = 0;
        bool FindD = true;
        if (FindB) {
            D_point = B_point + searchstepsize * unitForce;
            double E_LJ_D = E_LJ(D_point);
            while (E_LJ_D < E_LJ_B) {
                searchstepsize *= 1.2; //Increase the stepsize if D is not found in this iteration
                D_point = B_point + searchstepsize * unitForce;
                E_LJ_D = E_LJ(D_point);
                bracket_count++;
		//break the loop if we cannot find D in finite iterations
                if (bracket_count > 100) {
                    cout << "Cannot find point D" << endl;
                    FindD = false;
                    break;
                }
            }
        }

        // if we fail to find B or D, it is possible there's no local minima along the
        // negative gradient direction, in this case we move one step towards the negative
        // gradient to decrease energy, then use line search in the future steps
        if (!FindB || !FindD) {
            cout << "Cannot find B or D, using normal steepest descent algorithm" << endl;
	    searchstepsize = init_searchstepsize;
            arma::mat new_point = A_point + searchstepsize * unitForce;
            while (E_LJ(new_point) > E_LJ_A){
                searchstepsize /= 2;
                new_point = A_point - searchstepsize * unitForce;
            }
            A_point = new_point;
            A_point.print("new_point");
            Convert_corr_to_atoms(opt_Atoms, A_point);
            printf("current energy: %.3e\n", E_LJ(opt_Atoms));
            F_LJ_central_difference(cForce, opt_Atoms, fdstepsize);
            cForce.print("Central Difference Force");
	    old_e = cur_e;
	    cur_e = E_LJ(opt_Atoms);
            continue;
        }
        else { //If we find both B and D, start the Golden section search
            cout << "Start golden section search" << endl;
            AB = arma::norm(B_point - A_point, "fro");
            BD = arma::norm(D_point - B_point, "fro");
            AD = arma::norm(D_point - A_point, "fro");
            if (AB < BD) {
                C_point = D_point + golden_ratio * (A_point - D_point);
            }
            else {
                C_point = B_point;
                B_point = A_point + golden_ratio * (D_point - A_point);
            }
            while (AD > golden_tol) {
                if (E_LJ(B_point) > E_LJ(C_point)) {
                    A_point = B_point;
                    B_point = C_point;
                } else {
                    D_point = C_point;
                }
                AB = arma::norm(B_point-A_point, "fro");
                BD = arma::norm(D_point-B_point, "fro");
                AD = arma::norm(D_point-A_point, "fro");
                if (AB < BD) {
                    C_point = D_point + golden_ratio*(A_point-D_point);
                } else {
                    C_point = B_point;
                    B_point = A_point + golden_ratio*(D_point-A_point);
                }
            }
	    //renew the point according to the relative energy between B and C
            if (E_LJ(B_point) > E_LJ(C_point)) {
                C_point.print("new_point");
                Convert_corr_to_atoms(opt_Atoms, C_point);
                printf("current energy: %.3e\n", E_LJ(opt_Atoms));
            } else {
                B_point.print("new_point");
                Convert_corr_to_atoms(opt_Atoms, B_point);
                printf("current energy: %.3e\n", E_LJ(opt_Atoms));
            }
            F_LJ_central_difference(cForce, opt_Atoms, fdstepsize);
            cForce.print("Central Difference Force");    
	    old_e = cur_e;
	    cur_e = E_LJ(opt_Atoms); 
        }
    }

    cout << "Total iterations: " << count << endl;
    printf("Final energy: %.3e\n", E_LJ(opt_Atoms));
}

