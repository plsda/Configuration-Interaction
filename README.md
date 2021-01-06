# Configuration-Interaction

Davidson's algorithm for finding the lowest eigenvalue(s) of diagonally dominant (real symmetric) matrices. For use in the CI method.

Examples for small matrices(still fails for larger matrices):

Output for matrix <img src="https://render.githubusercontent.com/render/math?math=\begin{pmatrix}10&0.01&0&00&.0001\\0.01&20.1&0&0.02&0\\0&0&40.3&0&0\\0&0.02&0&39.8&0.1\\0.0001&0&0&0.1&30\end{pmatrix}">
![sc1](sampleOutput.png)

Two iterations yield lowest eigenvalue 7.7827176025 for <img src="https://render.githubusercontent.com/render/math?math=\begin{pmatrix}11.711783196508682&0.0& 0.71016571550645469&4.4279305398724329&0.0&0.0\\13.160405285805840&0.0&0.0&0.0&0.71016571550645469&0.0\\18.901150547807244&0.0&5.4240546891689814&4.4279305398724329&0.0&0.0\\13.332956938383129&1.0312204351939451&0.0&0.0&5.4240546891689814&1.0312204351939451&14.0\end{pmatrix}">, while full diagonalization gives 7.7548029175 as the lowest eigenvalue.



References:
[1] The iterative calculation of a few of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, R. Davidson, 1975, J. Chem. Phys.
[2] Numerical Methods For Large Eigenvalue Problems, Yousef Saad, 2011, p.203-204
[3] The simultaneous expansion method for the iterative solution of several of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices, B. Liu, Technical Report LBL-8158, Lawrence Berkeley Laboratory, University of California, Berkeley, 1978. 
[4] The Configuration Interaction Method: Advances in Highly Correlated Approaches, C. David Sherrillt and Henry F. Schaefer I11, 1999, p.182-186



