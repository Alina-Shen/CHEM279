#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <vector>

struct Atom
{
    int AtomicNumber;
    double r[3]; //xyz coordinates
    friend std::ostream &operator<<(std::ostream &os, const Atom p)
    {
        os << p.AtomicNumber << "(" << p.r[0] << ", "
           << p.r[1] << ", " << p.r[2] << ")";
        return os;
    }
};

void read_atoms_from_file(std::vector<Atom> &Atoms, const std::string &fname);
void write_atoms_to_file(const std::vector<Atom> &Atoms, const std::string &fname);

#endif // ATOM_H
