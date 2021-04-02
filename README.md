# Configuration-Interaction


The CI method is a method for solving the many-body Schrödinger eqution based on the Ritz method aka the method of linear variations.
This program solves numerically the first eigenstate of the 1D infinite potential well with two interacting electrons. More precisely,
we solve the time-independent Schrödinger equation

![sc1](https://github.com/plsda/Configuration-Interaction/blob/main/tise.png)

by expanding the wavefunction in a chosen basis and solving the eigenvalue problem

![sc2](https://github.com/plsda/Configuration-Interaction/blob/main/eigenEq.png)

where the Hamiltonian is now represented as a matrix in a basis formed by a finite set of functions

![sc3](https://github.com/plsda/Configuration-Interaction/blob/main/basisFunction.png)

There are many other possible choices of basis set, but this choice is easily motivated starting from the independent-electron model 
and by the fact that we target a symmetric spatial state.

One of the components needed for the CI method is an efficient iterative diagonalization procedure. One such method that
has emerged from the needs of quantum chemistry is Davidson's method and its variations. 

Davidson's algorithm is a numerical algorithm for finding the lowest(highest) eigenvalue(s) of diagonally dominant (real symmetric) matrices. Davidson's algorithm
and its derivatives are widely used in configuration interaction calculations. This program includes an implementation of Davidson's method
for the lowest eigenpair of a matrix satisfying the conditions mentioned above.

Davidson's method example output:

Let A be a N by N matrix defined by 

![sc4](https://github.com/plsda/Configuration-Interaction/blob/main/code/eqn.PNG)

(the example used in [3])

For N = 1000

![sc5](https://github.com/plsda/Configuration-Interaction/blob/main/code/sampleOutput1000.PNG)

For N = 10000

![sc6](https://github.com/plsda/Configuration-Interaction/blob/main/code/sampleOutput10000.PNG)


It turns out that the ground state energy of the system is approximately 2.2971 a.u.

![sc7](https://github.com/plsda/Configuration-Interaction/blob/main/CIExampleOutput.PNG)



Dependencies:

  - Intel MKL

  - GSL (not strictly necessary for the core functionality)
  

Building requires C++ 17 or newer. The program has only been tested on Windows 10 and build.bat builds the program using MSVC, but the code should be platform independent as long
as the dependencies listed are supported on the platform.


References:

[1] The iterative calculation of a few of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, R. Davidson, 1975, J. Chem. Phys.

[2] Numerical Methods For Large Eigenvalue Problems, Yousef Saad, 2011, p.203-204

[3] The simultaneous expansion method for the iterative solution of several of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, B. 
Liu, Technical Report LBL-8158, Lawrence Berkeley Laboratory, University of California, Berkeley, 1978. 

[4] The Configuration Interaction Method: Advances in Highly Correlated Approaches, C. David Sherrillt and Henry F. Schaefer I11, 1999, p.182-186



