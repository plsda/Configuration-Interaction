#include <stdio.h>
#include <stdint.h>
#include <memory>
#include <cfloat>
#include <ctime>
#include <cmath>

//#define CI_DEBUG //Define to enable debug functions and macros such as debugPrintMatrix and assert

#ifndef CI_DEBUG
   #define NDEBUG
#endif

#include "mathutils.h"
#include "ci_integrals.h"

int main()
{

   // In atomic units
   r64 wellWidth = 4.0;


   // Simply take the first *dim* basis functions, no fancy stuff here
   size_t maxWavenumber = 2;
   size_t dim = choose2(maxWavenumber);
   Matrix hamiltonian(dim, dim);

   size_t maxSteps = 10;
   r64* rombergArray = (r64*)malloc(4*maxSteps*sizeof(*rombergArray));

   size_t i = 1;
   size_t j = 2;
   size_t k = 1;
   size_t l = 2;

   for(int row = 0; row < dim; row++)
   {
      for(int column = row; column < dim; column++)
      {

         // Integrals are recalculated each time
         r64 oneElectronContribution = ((i == k)*(j == l) - (j == k)*(i == l));

         r64 twoElectronContribution = 0;

         int ijParity = ((i ^ j) & 1);
         int klParity = ((k ^ l) & 1);
         if((ijParity && klParity) || (!ijParity && !klParity))
         {
            clock_t clockBegin = clock();
            twoElectronContribution = getTwoElectronIntegral(i, j, k, l, wellWidth, rombergArray, maxSteps);
            clock_t clockEnd = clock();

            r64 intTime = ((r64)(clockEnd - clockBegin)) / CLOCKS_PER_SEC;

            //if(intTime > 0.03)
            //{
            //   printf("%fs \t\t i = %zd, j = %zd, k = %zd, l = %zd, \t %f\n", intTime, i, j, k, l, twoElectronContribution);
            //}
         }

         // Normalization from one-electron wavefunctions & physical constants
         oneElectronContribution *= 1.0/2.0*(PI64*PI64/(wellWidth*wellWidth)) * (k*k + l*l);
         twoElectronContribution *= 2.0/(wellWidth*wellWidth);

         hamiltonian(row, column) = oneElectronContribution + twoElectronContribution;
         hamiltonian(column, row) = hamiltonian(row, column);

         debugPrintMatrix(hamiltonian);

         if(j++ == maxWavenumber)
         {
            i++;
            j = i + 1;
         }

      }

      if(l++ == maxWavenumber)
      {
         k++;
         l = k + 1;
      }

      i = k;
      j = l;

      printf("\n\n %d/%zd \n\n", row + 1, dim);
   }


   free(rombergArray);

#if 1
   const size_t matrixDim = dim;
   const size_t maxIterations = 20;
   const size_t guessVectorCount = 1;
   const r64 tolerance = 10e-12; //Squared

   //NOTE: Hermitian(real symmetric) matrix
   Matrix subspaceBasis(maxIterations + guessVectorCount, matrixDim, guessVectorCount, matrixDim); //Rows of the matrix form the vectors and components of a vector are stored contiguously. 

   generateGuessVectors(hamiltonian, subspaceBasis, guessVectorCount, GUESS_BLOCK);

   debugPrintMatrix(hamiltonian);
   debugPrintMatrix(subspaceBasis);

   DavidsonInfo info = {maxIterations, guessVectorCount, 0, 0, tolerance};

   r64 eigenvalue = findLowestEigenpair(hamiltonian, subspaceBasis, info);
   Vector eigenvector(subspaceBasis(0));
   printMatrix(eigenvector);

   FILE* outputFile = fopen("E:\\Thesis\\code_visualization\\code\\CI_vector.txt", "wb");
   fwrite(eigenvector.elements, eigenvector.dim*sizeof(*eigenvector.elements), 1, outputFile);
   fclose(outputFile);

   // for testing
   
   r64* eigenvalues = (r64*)calloc(matrixDim, sizeof(r64));
   LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', hamiltonian.rows, hamiltonian.elements, hamiltonian.columns, eigenvalues); 
   for(int n = 0; n < dim; n++)
   {
      hamiltonian(0, n) = hamiltonian(n, 0);
   }
   Vector actualEigenvector(hamiltonian(0));
   printMatrix(actualEigenvector);
   //FILE* outputFile = fopen("CI_vector.txt", "wb");
   //fwrite(actualEigenvector.elements, actualEigenvector.dim*sizeof(*actualEigenvector.elements), 1, outputFile);
   //fclose(outputFile);
   //

   maxSteps = 15;
   rombergArray = (r64*)malloc(4*maxSteps*sizeof(*rombergArray));
   size_t rombergEvals = 0;
   auto varIntegrand = [eigenvalue, eigenvector, wellWidth, maxSteps, tolerance, rombergArray, &rombergEvals](r64 x1)
   {
      auto inner = [eigenvalue, eigenvector, wellWidth, x1, &rombergEvals](r64 x2)
      {

         r64 value = 0;
         if(!floatCompare(x1, x2))
         {
            r64 wavefunction = 0;
            r64 sum1 = 0;
            r64 sum2 = 0;

            size_t coeffIndex = 0;
            size_t maxWavenumber = (size_t)(1 + (int)sqrt(1 + 8*eigenvector.dim)) / 2;

            for(size_t i = 1; i <= maxWavenumber; i++)
            {
               for(size_t j = i + 1; j <= maxWavenumber; j++)
               {
                  r64 n = i*PI64/wellWidth;
                  r64 m = j*PI64/wellWidth;

                  r64 cn = eigenvector.elements[coeffIndex++];
                  r64 waveTerm = cn*SIGN(x1 - x2)*(sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1));
                  wavefunction += waveTerm;
                  sum1 += i*i*waveTerm;
                  sum2 += j*j*waveTerm;
               }
            }

            wavefunction *= 2.0/(wellWidth*sqrt(2.0));
            value = 1.0/(wellWidth*sqrt(2.0))*PI64*PI64/(wellWidth*wellWidth)*(sum1 + sum2) + wavefunction/(x1 - x2);

            value -= eigenvalue*wavefunction;
            value *= value;

            return value;
         }

         return value;
      };

      //r64 value = 2.0*rombergIntegrate(inner, 0.0, x1, rombergArray, maxSteps, tolerance);
      //r64 value = 2.0*adaptiveSimpson(inner, 0.0, x1, tolerance, maxSteps);
      r64 value = 2.0*tanhSinhQuadrature(inner, 0.0, x1, 1e-12);

      return value;
   };

   //r64 localEnergyVariance = rombergIntegrate(varIntegrand, 0, wellWidth, rombergArray + 2*maxSteps, maxSteps, tolerance);
   //r64 localEnergyVariance = adaptiveSimpson(varIntegrand, 0.0, wellWidth, tolerance, maxSteps);
   r64 localEnergyVariance = tanhSinhQuadrature(varIntegrand, 0.0, wellWidth, 1e-12);

   printf("Basis size: %zd\n\n", dim);
   printf("Calculated eigenvalue: %.10f\n", eigenvalue);
   printf("Actual eigenvalue: %.10f\n", eigenvalues[0]);
   printf("Variance: %.10f\n",localEnergyVariance);
   printf("\n\nCompleted iterations: %i", info.completedIterations);
   printf("\n\nDavidson time: %fs", info.elapsedTime);
#endif

   getchar();
   return 0;
}


