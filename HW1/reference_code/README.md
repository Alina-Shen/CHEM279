# Description
This is the reference code for HW1. It reads atomic coordinates from a file, calculates the Lennard-Jones energy, evaluates analytical and numerical forces on atoms, and performs optimization using the steepest descent method with options for a standard approach or enhanced with a golden section line search.

#Build the code on Datahub
`g++ -o hw1 main.cpp Atom.cpp eval.cpp -larmadillo`

#Usage
`./hw1 <filename>`
`<filename>`: The path to the input file containing atomic numbers and xyz coordinates.

