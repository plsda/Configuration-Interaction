#include <stdio.h>
#include <stdint.h>
#include <cfloat>
#include <ctime>
#include <cmath>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>

#include "mathutils.h"

int main()
{

	srand((int)time(0));
	r64 scale = 1.0/(r64)RAND_MAX;

	const int matrixDim = 5;
	const int maxIterations = 5;
	const int guessVectorCount = 2;

	//NOTE: Hermitian(real symmetric) matrix
	//Row-major matrices
	Matrix testMatrix(matrixDim, matrixDim);
	r64 diagonal[matrixDim];
	r64 diagonalForSort[matrixDim];

	//Test matrix generation
	for(int column = 0; column < matrixDim; column++)
	{
		for(int row = column; row < matrixDim; row++)
		{
			int randomNumber = rand();
			r64 matrixElement = (randomNumber % 2 == 0 ? 0 : scale*randomNumber);
			testMatrix[row][column] = matrixElement;
			testMatrix[column][row] = matrixElement;
		}

		r64& diagonalElement = testMatrix[column][column];
		diagonalElement += column + 10;
		diagonal[column] = diagonalElement;
		diagonalForSort[column] = diagonalElement;
	}

	Matrix subspaceBasis(maxIterations + guessVectorCount, matrixDim, guessVectorCount, matrixDim);	//Rows of the matrix form the vectors and components of a vector are stored contiguously. 

	//Generate unit vectors according to the largest diagonal elements as the initial guess vectors
	for(int vectorCount = 0; vectorCount < guessVectorCount; vectorCount++)
	{
		r64 largestDiagonal = DBL_MIN+1;
		int largestDiagonalIndex = 0;

		for(int i = 0; i < matrixDim; i++)
		{
			if(diagonalForSort[i] >= largestDiagonal)
			{
				largestDiagonal = diagonalForSort[i];
				largestDiagonalIndex = i;
			}
		}

		diagonalForSort[largestDiagonalIndex] = DBL_MIN; 
		subspaceBasis[vectorCount][largestDiagonalIndex] = 1.0;
	}

	//Allocate workspace
	size_t subspaceMaxDim = maxIterations + guessVectorCount;
	Matrix sigmaMatrix(subspaceMaxDim , subspaceBasis.columns); //Rows of sigmaMatrix are the sigma vectors
	Matrix projectedMatrix(subspaceMaxDim, subspaceMaxDim);
	Matrix eigenLHS(1, subspaceBasis.columns);					   //LHS of Av = b
	Matrix eigenRHS(1, subspaceBasis.columns);					   //RHS of Av = b, the Ritz vector
	r64* eigenvalues = (r64*)calloc(subspaceMaxDim, sizeof(r64));


	//### Begin Davidson ###
	
	r64 tolerance = 0.00001;
	r64 eigenvalue = 0;
	
	for(int iteration = 0; iteration < maxIterations; iteration++)
	{
		r64 absErrorUpperBound = 0;

		mmul(sigmaMatrix, subspaceBasis, testMatrix, false, true); 				
		mmul(projectedMatrix, sigmaMatrix, subspaceBasis, false, true); //projectedMatrix is column-major for convenience

		assert(projectedMatrix.rows == projectedMatrix.columns);
		Matrix& eigenvectors = projectedMatrix;

		//NOTE: Out-of-place calculation of eigenvectors and eigenvalues is not possible with MKL(or LAPACK) -- either copy projectedMatrix or reconstruct it each iteration
		int info = (int)LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', projectedMatrix.rows, projectedMatrix.elements, projectedMatrix.columns, eigenvalues);

		eigenvalue = eigenvalues[0]; 		  		 		 //Seek the lowest eigenvalue
		eigenvectors.dim = {1, projectedMatrix.rows}; //Pick only the 1st eigenvector since we seek only the lowest eigenvalue

		mmul(eigenLHS, eigenvectors, sigmaMatrix);	
		mmul(eigenRHS, eigenvectors, subspaceBasis); 
		Vector correctionVector(matrixDim, subspaceBasis[subspaceBasis.rows]);

		for(int i = 0; i < correctionVector.dim; i++)
		{
			r64 residual = eigenLHS[0][i] - eigenvalue*eigenRHS[0][i];
			r64 preconditioner = 1.0/(eigenvalue - diagonal[i]);
			correctionVector[i] = preconditioner*residual;
			absErrorUpperBound += residual*residual;
		}

		if(sqrt(absErrorUpperBound) < tolerance && iteration >= 2)
		{
			break;
		}

		subspaceBasis.rows++;
		orthonormalizeVectors(subspaceBasis, eigenvalues);

	}

	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', testMatrix.rows, testMatrix.elements, testMatrix.columns, eigenvalues);

	printf("Calculated eigenvalue:    %.*f\nActual lowest eigenvalue: %.*f", 10, eigenvalue, 10, eigenvalues[0]);

	getchar();
	return 0;
}
