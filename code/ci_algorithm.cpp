#include <stdio.h>
#include <stdint.h>
#include <memory>
#include <cfloat>
#include <ctime>
#include <cmath>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>

//#define CI_DEBUG //Define to enable debug functions and macros such as printMatrix and assert; Undefine to disable

#include "mathutils.h"

enum GUESS_TYPE
{
	GUESS_SMALLEST,	//Guess vectors are chosen according to the smallest diagonal elements
	GUESS_LARGEST,		//Guess vectors are chosen according to the largest diagonal elements
	GUESS_BLOCK,		//Guess vectors are the eigenvalues of a leading block of the matrix
};

void generateTestMatrix(Matrix& result, r64* diagonal, r64 sparseness, r64 diagonalWeight);
template <typename M> void generateTestMatrix(Matrix& result, r64* diagonal, M matrix);
void generateGuessVectors(Matrix& matrix, r64* diagonal, Matrix& guessVectors, int guessVectorCount, GUESS_TYPE guessType);

int main()
{

	srand((int)time(0));
	r64 scale = 1.0/(r64)RAND_MAX;

	const int matrixDim = 10000;
	const int maxIterations = 10;
	const int guessVectorCount = 4;

	//NOTE: Hermitian(real symmetric) matrix
	//Row-major matrices
	Matrix testMatrix(matrixDim, matrixDim);
	Matrix subspaceBasis(maxIterations + guessVectorCount, matrixDim, guessVectorCount, matrixDim);	//Rows of the matrix form the vectors and components of a vector are stored contiguously. 
	r64 diagonal[matrixDim];
	r64 sortedDiagonal[matrixDim];

	//generateTestMatrix(testMatrix, diagonal, scale, 10);
	generateTestMatrix(testMatrix, diagonal, 
							 [](int row, int column) { 
							 	return (row == column) ? (row < 5 ? (1 + 0.1*row) : (2*(row + 1) - 1)) : 1.0; 
							 });
	generateGuessVectors(testMatrix, diagonal, subspaceBasis, guessVectorCount, GUESS_BLOCK);



	printMatrix(testMatrix, true);
	printMatrix(subspaceBasis, true);

	//Allocate workspace
	size_t subspaceMaxDim = maxIterations + guessVectorCount;
	Matrix sigmaMatrix(subspaceMaxDim , subspaceBasis.columns); //Rows of sigmaMatrix are the sigma vectors
	Matrix projectedMatrix(subspaceMaxDim, subspaceMaxDim);		
	Matrix projectedMatrixCopy(subspaceMaxDim, subspaceMaxDim);		
	Matrix eigenvectors(subspaceMaxDim, subspaceMaxDim);		
	Matrix eigenLHS(1, subspaceBasis.columns); //LHS of Av = b
	Matrix eigenRHS(1, subspaceBasis.columns); //RHS of Av = b, the Ritz vector
	r64* eigenvalues = (r64*)calloc(matrixDim, sizeof(r64));

	//Rows of sigmaMatrix are the sigma vectors
	mmul(sigmaMatrix, subspaceBasis, testMatrix, false, true); 				
	mmul(projectedMatrix, subspaceBasis, sigmaMatrix, false, true);
	assert(projectedMatrix.rows == projectedMatrix.columns);

	printMatrix(sigmaMatrix, true);
	printMatrix(projectedMatrix, true);

	//Pack projectedMatrix(store lower triangular part)
	size_t flatIndex = 0;
	for(int row = 0; row < projectedMatrix.rows; row++)
	{
		for(int column = 0; column < row+1; column++)
		{
			projectedMatrix[flatIndex++] = projectedMatrix(row,column);
		}
	}

	printMatrix(projectedMatrix, true);
	printPacked(projectedMatrix, true);

	r64* scratchpad = (r64*)malloc(100*sizeof(r64));

	//### Begin Davidson ###

	r64 tolerance = 0.00001; //Squared
	r64 eigenvalue = 0;
	int iteration = 0;

	clock_t clockBegin = clock();

	for(iteration; iteration < maxIterations; iteration++)
	{
		r64 absErrorUpperBound = 0;

		eigenvectors.dim = projectedMatrix.dim;
		size_t elementsInLowerTriangle = projectedMatrix.rows*(projectedMatrix.rows + 1)/2;
		cblas_dcopy(elementsInLowerTriangle, projectedMatrix.elements, 1, projectedMatrixCopy.elements, 1);

		int info = (int)LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'L', eigenvectors.rows, projectedMatrixCopy.elements, 
												eigenvalues, eigenvectors.elements, eigenvectors.columns);

		printMatrix(eigenvectors, true);

		eigenvalue = eigenvalues[0]; 		  		 		 //Seek the lowest eigenvalue
		cblas_dcopy(projectedMatrix.rows, eigenvectors.elements, eigenvectors.rows, eigenvectors.elements, 1);
		eigenvectors.dim = {1, projectedMatrix.rows}; //Pick only the 1st eigenvector since we seek only the lowest eigenvalue

		printMatrix(eigenvectors, true);

		mmul(eigenLHS, eigenvectors, sigmaMatrix);	
		mmul(eigenRHS, eigenvectors, subspaceBasis); 
		Vector correctionVector(subspaceBasis(subspaceBasis.rows));

		printMatrix(eigenLHS, true);
		printMatrix(eigenRHS, true);
		printMatrix(correctionVector, true);

		for(int i = 0; i < correctionVector.dim; i++)
		{
			r64 residual = eigenLHS(0,i) - eigenvalue*eigenRHS(0,i);
			r64 preconditioner = 1.0/(eigenvalue - diagonal[i]);
			correctionVector[i] = preconditioner*residual;
			absErrorUpperBound += residual*residual;
		}

		if(absErrorUpperBound < tolerance)
		{
			break;
		}

		orthonormalizeVectorAgainst(correctionVector, subspaceBasis);
		subspaceBasis.rows++;

		Vector newBasisVector(subspaceBasis(subspaceBasis.rows-1));

		Matrix basisNorm = mmul(subspaceBasis, subspaceBasis, false, true);

		printMatrix(basisNorm, true);

		Vector newSigmaVector(sigmaMatrix(sigmaMatrix.rows++));
		newSigmaVector = testMatrix*newBasisVector; 

		//Projected matrix is symmetric, so we only need to store the new row(using lower triangular part)
		Vector projectedMatrixNewRow(projectedMatrix.columns, &projectedMatrix[elementsInLowerTriangle]); 
		projectedMatrixNewRow = subspaceBasis*newSigmaVector;
		projectedMatrix.rows++;
		projectedMatrix.columns++;

		printMatrix(subspaceBasis, true);
		printMatrix(projectedMatrixNewRow, true);
		printMatrix(projectedMatrix, true);
		printPacked(projectedMatrix, true);
	}

	clock_t clockDavidsonEnd = clock();

	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'L', testMatrix.rows, testMatrix.elements, testMatrix.columns, eigenvalues); //NOTE: Should use dsyer and only calculate the lowest eigenvalue to save time

	clock_t clockEnd = clock();
	r64 davidsonTime = ((r64)(clockDavidsonEnd - clockBegin)) / CLOCKS_PER_SEC;
	r64 fullDiagonalizationTime = ((r64)(clockEnd - clockDavidsonEnd)) / CLOCKS_PER_SEC;

	printf("Calculated eigenvalue:    %.*f\nActual lowest eigenvalue: %.*f", 10, eigenvalue, 10, eigenvalues[0]);
	printf("\n\nCompleted iterations: %i", iteration);
	printf("\n\nDavidson time: %fs", davidsonTime);
	printf("\n\ndsyev time: %fs", fullDiagonalizationTime);

	getchar();
	return 0;
}

void findLowestEigenvalue(Matrix& matrix, r64* diagonal)
{
}

void generateTestMatrix(Matrix& result, r64* diagonal, r64 sparseness, r64 diagonalWeight)
{
	size_t matrixDim = result.rows;
	for(int column = 0; column < matrixDim; column++)
	{
		for(int row = column; row < matrixDim; row++)
		{
			int randomNumber = rand();
			r64 matrixElement = sparseness*randomNumber;
			result(row, column) = matrixElement;
			result(column,row) = matrixElement;
		}

		r64& diagonalElement = result(column,column);
		diagonalElement += column + diagonalWeight;
		diagonal[column] = diagonalElement;
	}
}

template <typename M>
void generateTestMatrix(Matrix& result, r64* diagonal, M matrix)
{
	size_t matrixDim = result.rows;
	for(int column = 0; column < matrixDim; column++)
	{
		for(int row = column; row < matrixDim; row++)
		{
			result(row, column) = matrix(row, column);
			result(column,row) = result(row, column);
		}
		diagonal[column] = result(column, column);
	}
}

void generateGuessVectors(Matrix& matrix, r64* diagonal, Matrix& guessVectors, int guessVectorCount, GUESS_TYPE guessType)
{
	if(guessType == GUESS_SMALLEST || guessType == GUESS_LARGEST)
	{
		int* largestDiagonalIndices= (int*)malloc(guessVectorCount*sizeof(int));
		r64 compareLimitValue;
		r64 defaultCompare;
		if(guessType == GUESS_SMALLEST)
		{
			compareLimitValue = DBL_MAX;
			defaultCompare = compareLimitValue - 1.0;
		}
		else
		{
			compareLimitValue = DBL_MIN;
			defaultCompare = compareLimitValue + 1.0;
		}
		
		r64* sortedDiagonal = (r64*)malloc(matrix.rows*sizeof(r64));
		memcpy(sortedDiagonal, diagonal, matrix.rows*sizeof(r64));

		for(int vectorCount = 0; vectorCount < guessVectorCount; vectorCount++)
		{
			r64 compare = defaultCompare;
			int compareIndex = 0;

			for(int i = 0; i < matrix.rows; i++)
			{
				if((guessType == GUESS_SMALLEST && sortedDiagonal[i] <= compare)||
					(guessType == GUESS_LARGEST && sortedDiagonal[i] >= compare))
				{
					compare = sortedDiagonal[i];
					compareIndex = i;
				}
			}

			sortedDiagonal[compareIndex] = compareLimitValue;
			largestDiagonalIndices[vectorCount] = compareIndex;
		}


		sort(largestDiagonalIndices, guessVectorCount);
		for(int i = 0; i < guessVectorCount; i++)
		{
			guessVectors(i,largestDiagonalIndices[i]) = 1.0;
		}

		free(largestDiagonalIndices);
		free(sortedDiagonal);
	}
	else
	{
		for(int row = 0; row < guessVectorCount; row++)
		{
			for(int column = 0; column < guessVectorCount; column++)
			{
				guessVectors(row, column) = matrix(row,column);
				guessVectors(column, row) = guessVectors(row, column);
			}
		}

		printMatrix(guessVectors, true);
		r64* eigenvalues = (r64*)malloc(guessVectorCount*sizeof(r64));
		int info = (int) LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', guessVectorCount, guessVectors.elements, guessVectors.columns, eigenvalues);

		for(int column = 0; column < guessVectorCount; column++)
		{
			for(int row = column; row < guessVectorCount; row++)
			{
				printMatrix(guessVectors, true);
				r64 element = guessVectors(row,column);
				guessVectors(row,column) = guessVectors(column,row);
				guessVectors(column,row) = element;
			}
		}
	}
}


