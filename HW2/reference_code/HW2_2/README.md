# Description
This is the reference code for HW2.2. It reads two shells of primitive gaussian type orbitals and returned the analytical overlap integrals between the shells.

# Build the code on Datahub
```
g++ -I. -c util.cpp; g++ -I. -c Shell.cpp; g++ -o hw2_2 -I. hw2_2_main.cpp Shell.o util.o -larmadillo
```

# Usage
```
./hw2_2 <filename>
```
`<filename>`: The path to the input file containing information of 2 shells of primitive gaussian type orbitals.

