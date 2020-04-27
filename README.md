# LDPCGaussian
Optimized version of Gaussian elimination meant for solving LDPC codes.

The source code can be found under the name Gaussianfields.java. You can run the Gaussianfields.exe program from command line like this:

Gaussianfields.exe input.txt output.txt

---

The input format is as follows:

exponent m\
decimal p\
equations N, unknowns M\
Matrix A

---

An example input would be:

4\
19\
4, 4\
10 15 0 0 3\
0 0 1 6 3\
0 10 0 6 8\
13 0 6 0 13

---

Note that the last column of matrix A contains b values. Input and output files have to be in the same directory

---

The output format is as follows:

vector x

---

An example output would be:

11 8 4 6

---

For any additional questions please contact the author (owner of repository).
