# Configuration-Interaction

Davidson's algorithm is a numerical algorithm for finding the lowest eigenvalue(s) of diagonally dominant (real symmetric) matrices. Davidson's algorithm
and its derivatives are widely used in configuration interaction calculations.

Example output:


Let A be a N by N matrix defined by 

<img src="https://latex.codecogs.com/png.latex?A_{ij} = A_{ji} = 1, \\A_{ii} = 1 + 0.1(i-1),\quad\text{i = 1, 2..,min(5, N)}, \\A_{ii} = 2i - 1,\quad\text{otherwise}"/> 
(the example used in [3])

For N = 1000

![sc1](https://github.com/plsda/Configuration-Interaction/blob/main/code/sampleOutput1000.PNG)

For N = 10000

![sc1](https://github.com/plsda/Configuration-Interaction/blob/main/code/sampleOutput10000.PNG)





References:

[1] The iterative calculation of a few of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, R. Davidson, 1975, J. Chem. Phys.

[2] Numerical Methods For Large Eigenvalue Problems, Yousef Saad, 2011, p.203-204

[3] The simultaneous expansion method for the iterative solution of several of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, B. 
Liu, Technical Report LBL-8158, Lawrence Berkeley Laboratory, University of California, Berkeley, 1978. 

[4] The Configuration Interaction Method: Advances in Highly Correlated Approaches, C. David Sherrillt and Henry F. Schaefer I11, 1999, p.182-186



