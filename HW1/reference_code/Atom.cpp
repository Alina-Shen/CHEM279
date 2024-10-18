#include <fstream>
#include <sstream>

#include "Atom.h"

using namespace std;

void read_atoms_from_file(std::vector<Atom> &Atoms, const std::string &fname)
{
    ifstream in;
    in.open(fname, ios::in);
    if (!in.is_open()) {
        throw std::runtime_error("Unable to open file: " + fname);
    }
    string line;
    getline(in, line); //get the number of atoms
    size_t num_atoms = stoull(line);
    //get the coordinates
    while (getline(in, line)){
        istringstream iss(line);
        Atom readatom;
        if (!(iss>>readatom.AtomicNumber>>readatom.r[0]>>readatom.r[1]>>readatom.r[2])){
            throw invalid_argument("There is some problem with the Atom format.");
        }
	//Check if the element is gold
        if (readatom.AtomicNumber != 79){
            throw invalid_argument("There are atoms other than Au atom.");
        }
        Atoms.push_back(readatom);
    }
    if (Atoms.size() != num_atoms){
        throw invalid_argument("Number of atoms are not consistent.");
    }
    in.close();
}

void write_atoms_to_file(const std::vector<Atom> &Atoms, const std::string &fname)
{
    ofstream out;
    out.open(fname, ios::out);
    out << Atoms.size() << endl;
    for (auto x : Atoms){
        out << x << endl;
    }
    out.close();
}
