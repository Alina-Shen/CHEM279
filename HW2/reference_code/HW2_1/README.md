# Description
This is the reference code for HW2.1. It reads two 1-D primitive gaussian type orbitals and returned the 1-D numerical overlap integral.

# Build the code on Datahub
```
g++ -I. -c Numerical.cpp; g++ -o hw2_1 -I. hw2_1_main.cpp Numerical.o -larmadillo
```

# Usage
```
./hw2_1 <filename>
```
`<filename>`: The path to the input file containing information of 2 1-D gaussian type orbitals.

