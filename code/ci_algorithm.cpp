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
   // Simply pick the first *dim* basis functions for the CI expansion, no fancy stuff here

   // In atomic units
   r64 wellWidth = 4.0;

   gsl_integration_glfixed_table* glTable = gsl_integration_glfixed_table_alloc(100);

   size_t maxWavenumber = 10;
   size_t dim = choose2(maxWavenumber);
   r64 intTolerance = 1e-10;
   r64* integralLookup = 0;

   //printf("Read integrals from file(y/n)? ");
   printf("\n\t1 - Compute on-the-fly\n\t2 - Read integrals from integrals.txt\n\t3 - Calculate and write integrals to file\n\n>");
   char choice;
   scanf("%c", &choice);
   switch(choice)
   {
      case '1':
      {
         printf("Max wavenumber: ");
         scanf("%zd", &maxWavenumber);
         printf("\nTolerance: ");
         scanf("%lf", &intTolerance);
      } break;

      case '2':
      {
         FILE* fileID = fopen("integrals.txt", "rb");

         IntegralsHeader integralsInfo;
         fread(&integralsInfo, sizeof(integralsInfo), 1, fileID);

         if(wellWidth != integralsInfo.wellWidth)
         {
            printf("The specified well width doesn't match the well width in the integrals file!\n");
            return 0;
         }

         maxWavenumber = integralsInfo.maxWavenumber;
         intTolerance = integralsInfo.tolerance;
         size_t lookupSize = integralsInfo.integralsWritten*sizeof(r64);
         integralLookup = (r64*)malloc(lookupSize);

         size_t bytesRead = fread(integralLookup, lookupSize, 1, fileID);
         if(bytesRead < 1)
         {
            printf("Could not read integrals from file!\n");
         }
         fclose(fileID);
      } break;

      case '3':
      {
         r64 tolerance = 1e-10;
         printf("Max wavenumber: ");
         scanf("%zd", &maxWavenumber);
         printf("\nTolerance: ");
         scanf("%lf", &tolerance);
         writeTwoElectronIntegrals(maxWavenumber, wellWidth, tolerance);

         return 0;
      } break;

      default:
      {
         return 0;
      } break;
   }

   dim = choose2(maxWavenumber);
   Matrix hamiltonian(dim, dim);

   size_t maxSteps = 10;
   r64* rombergArray = (r64*)malloc(4*maxSteps*sizeof(*rombergArray));

   size_t i = 1;
   size_t j = 2;
   size_t k = 1;
   size_t l = 2;
   size_t integralIndex = 0;

   printf("\n");
   for(int row = 0; row < dim; row++)
   {
      for(int column = row; column < dim; column++)
      {

         r64 oneElectronContribution = ((i == k)*(j == l) - (j == k)*(i == l));

         r64 twoElectronContribution = 0;

         int ijParity = ((i ^ j) & 1);
         int klParity = ((k ^ l) & 1);
         if((ijParity && klParity) || (!ijParity && !klParity))
         {
            if(choice == '2')
            {
               twoElectronContribution = integralLookup[integralIndex++];
            }
            else
            {
               //clock_t clockBegin = clock();
               twoElectronContribution = getTwoElectronIntegral(i, j, k, l, wellWidth, rombergArray, maxSteps, glTable, intTolerance);
               //clock_t clockEnd = clock();

               //r64 intTime = ((r64)(clockEnd - clockBegin)) / CLOCKS_PER_SEC;
               //printf("%f\n", intTime);

               //////if(intTime > 0.03)
               //////{
               //   printf("\n\t%fs \t\t i = %zd, j = %zd, k = %zd, l = %zd, \t %f\n", intTime, i, j, k, l, twoElectronContribution);
               //////}
            }
         }

         // Normalization from one-electron wavefunctions and other constants
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
   
      printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      printf("%d/%zd", row + 1, dim);
   }


   free(rombergArray);
   free(integralLookup);

   const size_t matrixDim = dim;
   const size_t maxIterations = 20;
   const size_t guessVectorCount = 1;
   const r64 diagTolerance = 1e-12; //Squared

   //NOTE: Hermitian(real symmetric) matrix
   Matrix subspaceBasis(maxIterations + guessVectorCount, matrixDim, guessVectorCount, matrixDim); //Rows of the matrix form the vectors and components of a vector are stored contiguously. 

   generateGuessVectors(hamiltonian, subspaceBasis, guessVectorCount, GUESS_BLOCK);

   debugPrintMatrix(hamiltonian);
   debugPrintMatrix(subspaceBasis);

   DavidsonInfo info = {maxIterations, guessVectorCount, 0, 0, diagTolerance};

   r64 eigenvalue = findLowestEigenpair(hamiltonian, subspaceBasis, info);
   Vector eigenvector(subspaceBasis(0));

   FILE* outputFile = fopen("CI_vector.txt", "wb");
   fwrite(eigenvector.elements, eigenvector.dim*sizeof(*eigenvector.elements), 1, outputFile);
   fclose(outputFile);

   maxSteps = 15;
   rombergArray = (r64*)malloc(4*maxSteps*sizeof(*rombergArray));
   size_t rombergEvals = 0;
   auto varIntegrand = [eigenvalue, eigenvector, wellWidth, maxSteps, intTolerance, rombergArray, &rombergEvals, glTable](r64 x1)
   {
      auto inner = [eigenvalue, eigenvector, wellWidth, x1, &rombergEvals](r64 x2)
      {

         r64 value = 0;
         if(!floatCompare(x1, x2))
         {
            r64 wavefunction = 0;
            r64 kineticTerm = 0;

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
                  kineticTerm += (i*i + j*j)*waveTerm; 
               }
            }

            wavefunction *= 2.0/(wellWidth*sqrt(2.0));
            value = 1.0/(wellWidth*sqrt(2.0))*PI64*PI64/(wellWidth*wellWidth)*kineticTerm + wavefunction/(x1 - x2);

            value *= value;

         }

         return value;
      };

      //r64 value = rombergIntegrate(inner, 0.0, x1, rombergArray, maxSteps, intTolerance);
      r64 value = tanhSinhQuadrature(inner, 0.0, x1, intTolerance);
      
      return value;
   };

   //r64 localEnergyVariance = 2.0*rombergIntegrate(varIntegrand, 0, wellWidth, rombergArray + 2*maxSteps, maxSteps, intTolerance) - eigenvalue*eigenvalue;
   r64 localEnergyVariance = 2.0*tanhSinhQuadrature(varIntegrand, 0.0, wellWidth, intTolerance) - eigenvalue*eigenvalue;


   printf("\n\nWell width: %f\n", wellWidth);
   printf("Basis size: %zd\n", dim);
   printf("Integral tolerances: %.2e\n\n", intTolerance);
   printf("Calculated eigenvalue: %.10f\n", eigenvalue);
   printf("Variance: %.10f\n",localEnergyVariance);
   printf("\n\nCompleted iterations: %i", info.completedIterations);
   printf("\n\nDavidson time: %fs", info.elapsedTime);
#endif

   getchar();
   getchar();
   return 0;
}


