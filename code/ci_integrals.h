#include <stdio.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


// Returns 2^(stage-1)-point refinement to the result of 2^(stage-1)-subinterval trapezoidal rule.
// Together with previous stages, the result can be used to compute trapezoidal rule with 2^stage subintervals.
template <typename f>
r64 trapezoidalRuleStage(f& integrand, r64 lowerBound, r64 upperBound, size_t stage)
{
   r64 result = 0;

   if(stage == 0)
   {
      r64 intervalLength = upperBound - lowerBound;
      result = 0.5*intervalLength*(integrand(lowerBound) + integrand(upperBound));
   }
   else
   {
      size_t newPoints = 1ULL << (stage - 1);
      r64 intervalLength = (upperBound - lowerBound) / ((r64)(2*newPoints));
      for(int interval = 0; interval < newPoints; interval++)
      {
         result += integrand(lowerBound + (2*interval + 1)*intervalLength);
      }

      result *= intervalLength;
   }

   return result;
}

template <typename f>
r64 trapezoidalRule(f integrand, r64 lowerBound, r64 upperBound, size_t subdivisions)
{
   assert(lowerBound < upperBound);

   r64 result = 0.5*(integrand(lowerBound) + integrand(upperBound));
   r64 intervalLength = (upperBound - lowerBound) / (r64)subdivisions;

   for(int intervalCount = 1; intervalCount < subdivisions; intervalCount ++)
   {
      result += integrand(lowerBound + intervalCount*intervalLength);
   }

   result *= intervalLength;

   return result;
}

//NOTE: rombergArray must have length of, at least, 2*maxSteps
template <typename f>
r64 rombergIntegrate(f& integrand, r64 lowerBound, r64 upperBound, r64* rombergArray, size_t maxSteps, r64 absTolerance, size_t minSteps = 2) 
{
   r64 result = 0;
   r64* __restrict prevRow = rombergArray;
   r64* __restrict newRow = rombergArray+maxSteps;

   prevRow[0] = trapezoidalRuleStage(integrand, lowerBound, upperBound, 0);

   for(int row = 1; row < maxSteps; row++)
   {
      newRow[0] = 0.5*prevRow[0] + trapezoidalRuleStage(integrand, lowerBound, upperBound, row);

      r64 weight = 4.0;
      for(int column = 1; column <= row; column++)
      {
         newRow[column] = newRow[column-1] + (newRow[column-1] - prevRow[column-1]) / (weight - 1.0);
         weight *= 4.0;
      }

      result = newRow[row];
      if(fabs(newRow[row] - prevRow[row-1]) < absTolerance && row > minSteps) //NOTE: This check is CRUCIAL for avoiding early termination and getting garbage results!
      {
         break;
      }

      r64* temp = prevRow;
      prevRow = newRow;
      newRow = temp;
   }


   return result;
}

// A rudimentary implementation of tanh-sinh quadrature for one-parameter integrand
template <typename f>
r64 tanhSinhQuadrature(f& integrand, r64 lowerBound, r64 upperBound, r64 absTolerance) 
{
   static r64* abscissae;
   static r64* weights;
   static size_t storedNonnegativePointCount;
   size_t maxStages = 8; 
   const r64 windowSize = 1.0;

   if(!storedNonnegativePointCount)
   {
      size_t maxNonnegativePoints = 2000;
      r64 minSubintervalLength = windowSize/(1ULL << maxStages);
      abscissae = (r64*)malloc(maxNonnegativePoints*sizeof(*abscissae));
      weights = (r64*)malloc(maxNonnegativePoints*sizeof(*weights));

      for(storedNonnegativePointCount = 0; 
          storedNonnegativePointCount < maxNonnegativePoints; 
          storedNonnegativePointCount++)
      {
         r64 pointCosh = cosh(HPI64*sinh(storedNonnegativePointCount*minSubintervalLength));
         abscissae[storedNonnegativePointCount] = tanh(HPI64*sinh(storedNonnegativePointCount*minSubintervalLength));
         weights[storedNonnegativePointCount] = HPI64*cosh(storedNonnegativePointCount*minSubintervalLength)/(pointCosh*pointCosh);

         if(weights[storedNonnegativePointCount] < 1.0e-13)
         {
            storedNonnegativePointCount--;
            break;
         }
      }
   }


   r64 result = 0;
   size_t stride = 1ULL << maxStages;

   for(; stride >= storedNonnegativePointCount; stride >>= 1, maxStages--){}

   auto transformedIntegrand = [integrand, lowerBound, upperBound](r64 x)
   {
      r64 normalizedX = 0.5*(upperBound - lowerBound)*x + 0.5*(upperBound + lowerBound); // "normalized" to [-1, 1]

      r64 value = integrand(normalizedX);

      return value;
   };


   size_t stage;
   r64 subintervalLength = windowSize;
   result = subintervalLength*weights[0]*transformedIntegrand(0);

   //Evaluate the integrand at every first-stage abscissa
   size_t nodeIndex = stride;
   for(int i = 1; 
       nodeIndex < storedNonnegativePointCount; 
       i++, nodeIndex = i*stride)
   {
      r64 abscissa = abscissae[nodeIndex]; 
      result += weights[nodeIndex]*(transformedIntegrand(abscissa) + transformedIntegrand(-abscissa));
   }

   stride >>= 1;
   for(stage = 2; stage <= maxStages; stage++)
   {
      subintervalLength *= 0.5;
      //subintervalLength = windowSize / (1 << (stage-1));

      r64 stageSum = 0;
      nodeIndex = stride;
      //Now we only need the abscissas and weights with odd indices
      for(int i = 0; 
          nodeIndex < storedNonnegativePointCount;
          i++, nodeIndex = (2*i + 1)*stride)
      {
         r64 abscissa = abscissae[nodeIndex]; 
         stageSum += weights[nodeIndex]*(transformedIntegrand(abscissa) + transformedIntegrand(-abscissa));
      }

      if(subintervalLength*fabs(stageSum) < absTolerance && stage > 1)
      {
         stage--;
         break;
      }

      result += stageSum;
      stride >>= 1;
   }

   result = 0.5*(upperBound - lowerBound)*subintervalLength*result;

   return result;
}


struct IntegrandParametersInner
{
   r64 x1, n, m, p, q; 
};

struct IntegrandParametersOuter
{
   r64 n, m, p, q; 
   r64 lowerBound, upperBound;
   gsl_integration_glfixed_table* glTable;
};

r64 twoElectronIntegrandInner(r64 x2, void* params)
{
   IntegrandParametersInner* paramsIn = (IntegrandParametersInner*)params;
   r64 x1 = paramsIn->x1;
   r64 n = paramsIn->n;
   r64 m = paramsIn->m;
   r64 p = paramsIn->p;
   r64 q = paramsIn->q;

	r64 value = floatCompare(x1, x2) ? 0.0 : (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/fabs(x1 - x2);

	return value;
}

r64 twoElectronIntegrandOuter(r64 x1, void* params)
{
   IntegrandParametersOuter* paramsIn = (IntegrandParametersOuter*)params;

   IntegrandParametersInner paramsOut = {x1, paramsIn->n, paramsIn->m, paramsIn->p, paramsIn->q};

   gsl_function F;
   F.function = &twoElectronIntegrandInner;
   F.params = &paramsOut;

   // Using symmetry
   r64 result = gsl_integration_glfixed(&F, paramsIn->lowerBound, x1, paramsIn->glTable);

	return result;
}

// NOTE: Only for integrals that don't vanish by group theory.
r64 twoElectronGaussLegendre(size_t i, size_t j, size_t k, size_t l, r64 iWellWidth, r64 lowerBound, 
                             r64 upperBound, size_t maxNodes = 32, gsl_integration_glfixed_table* glTable = 0)
{
   if(!glTable)
   {
      glTable = gsl_integration_glfixed_table_alloc(maxNodes);
   }

   IntegrandParametersOuter params = 
      {PI64*iWellWidth*i, PI64*iWellWidth*j, PI64*iWellWidth*k, PI64*iWellWidth*l, 
       lowerBound, upperBound, glTable};

   gsl_function F;
   F.function = &twoElectronIntegrandOuter;
   F.params = &params;
   r64 result = gsl_integration_glfixed(&F, lowerBound, upperBound, glTable);

   return result;
}



#pragma pack(push, 1)
struct IntegralsHeader
{
   size_t maxWavenumber;
   size_t integralsWritten;
   size_t maxSteps;
   r64 tolerance;
   r64 wellWidth;
};
#pragma pack(pop)

// NOTE: Only for integrals that don't vanish by group theory.
r64 getTwoElectronIntegral(size_t i, size_t j, size_t k, size_t l, r64 wellWidth, r64* scratch, size_t maxSteps, 
                           gsl_integration_glfixed_table* glTable, r64 tolerance = 1e-8)
{
   r64 iWellWidth = 1.0/wellWidth;
   size_t rombergEvals = 0;
   r64 n = PI64*iWellWidth*i; 
   r64 m = PI64*iWellWidth*j; 
   r64 p = PI64*iWellWidth*k; 
   r64 q = PI64*iWellWidth*l;

   r64* rombergArray1 = scratch;
   r64* rombergArray2 = scratch + 2*maxSteps;

   r64 value = 0;


   auto integrand = [n, m, p, q, rombergArray1, maxSteps, tolerance, wellWidth, &rombergEvals](r64 x1)
   {

      auto inner = [n, m, p, q, x1](r64 x2)
      {
         r64 value = 0.0;

         if(!floatCompare(x1, x2))
         {
            value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/(x1 - x2);
         }

         return value;
      };

      //r64 value = rombergIntegrate(inner, 0.0, x1, rombergArray1, maxSteps, tolerance, rombergEvals, 2);
      r64 value = tanhSinhQuadrature(inner, 0.0, x1, tolerance);

      return value;
   };

   //value = 2.0*twoElectronGaussLegendre(i, j, k, l, iWellWidth, 0, wellWidth, 1000, glTable);
   //value = 2.0*rombergIntegrate(integrand, 0, wellWidth, rombergArray2, maxSteps, tolerance, rombergEvals);
   value = 2.0*tanhSinhQuadrature(integrand, 0.0, wellWidth, tolerance);


   return value;
}


int writeTwoElectronIntegrals(size_t maxWavenumber, r64 wellWidth, r64 tolerance)
{
   r64 iWellWidth = 1.0/wellWidth;
   size_t maxSteps = 20;

   r64* rombergArray = (r64*)malloc(4*maxSteps*sizeof(*rombergArray));
   gsl_integration_glfixed_table* glTable = gsl_integration_glfixed_table_alloc(100);

   FILE* outputFile = fopen("integrals.txt", "wb");

   if(outputFile)
   {

      size_t dim = choose2(maxWavenumber);
      size_t dimCounter = 1;
      r64* rowBuffer = (r64*)malloc(dim*sizeof(*rowBuffer));
      size_t rowCounter = 0;
      size_t integralCount = 0;

      size_t i = 1;
      size_t j = 2;
      size_t k = 1;
      size_t l = 2;

      fseek(outputFile, sizeof(IntegralsHeader), SEEK_SET);
      for(int row = 0; row < dim; row++)
      {
         for(int column = row; column < dim; column++)
         {
            // Skip integrals vanishing by group theory
            int ijParity = ((i ^ j) & 1);
            int klParity = ((k ^ l) & 1);
            if((ijParity && klParity) || (!ijParity && !klParity))
            {

               r64 n = PI64*iWellWidth*i; 
               r64 m = PI64*iWellWidth*j; 
               r64 p = PI64*iWellWidth*k; 
               r64 q = PI64*iWellWidth*l;

               r64 integralValue = getTwoElectronIntegral(i, j, k, l, wellWidth, rombergArray, maxSteps, glTable, tolerance);
               rowBuffer[rowCounter++] = integralValue;
            }

            if(j++ == maxWavenumber)
            {
               i++;
               j = i + 1;
            }

         }

         fwrite(rowBuffer, rowCounter*sizeof(*rowBuffer), 1, outputFile);
         integralCount += rowCounter;
         rowCounter = 0;

         if(l++ == maxWavenumber)
         {
            k++;
            l = k + 1;
         }

         i = k;
         j = l;

         printf("%zd/%zd\n", dimCounter++, dim);
      }

      IntegralsHeader header = {};
      header.maxWavenumber = maxWavenumber;
      header.integralsWritten = integralCount;
      header.maxSteps = maxSteps;
      header.tolerance = tolerance;
      header.wellWidth = wellWidth;
      fseek(outputFile, 0, SEEK_SET);
      fwrite(&header, sizeof(header), 1, outputFile);

      free(rombergArray);
      free(rowBuffer);
      fclose(outputFile);

      return 0;
   }
   else
   {
      return -1;
   }
}

template<typename F>
double adaptiveSimpson(F f, double lower, double upper, double tolerance, size_t maxSteps)
{
   double midPoint = 0.5*(lower + upper);
   double result = (upper - lower)/6.0 * (f(lower) + 4.0*f(midPoint) + f(upper));

   result = simpsonStep(f, 1, result, lower, upper, tolerance, maxSteps);

   return result;
}

template<typename F>
double simpsonStep(F f, size_t step, double previousEstimate, double lower, double upper, double tolerance, size_t maxSteps)
{
   double h = 0.5*(upper - lower);
   double midPoint = 0.5*(lower + upper);
   double fMid = f(midPoint);

   double leftSimpson = h/6.0 * (f(lower) + 4.0*f(0.5*(lower + midPoint)) + fMid);
   double rightSimpson = h/6.0 * (fMid + 4.0*f(0.5*(midPoint + upper)) + f(upper));

   double result = leftSimpson + rightSimpson;
   double diff = result - previousEstimate;

   if(fabs(diff) >= 15.0*tolerance && step < maxSteps)
   {
      result = simpsonStep(f, step + 1, leftSimpson, lower, midPoint, 0.5*tolerance, maxSteps) + 
               simpsonStep(f, step + 1, rightSimpson, midPoint, upper, 0.5*tolerance, maxSteps);
   }
   else
   {
      result += diff / 15.0;
   }

   return result;
}

void testIntegration(int i, int j, int k, int l, r64 wellWidth, r64 tolerance)
{
   r64 iWellWidth = 1.0/wellWidth;
   r64 n = PI64*iWellWidth*i; 
   r64 m = PI64*iWellWidth*j; 
   r64 p = PI64*iWellWidth*k; 
   r64 q = PI64*iWellWidth*l;

   size_t maxSteps = 20;
   r64* rombergArray1 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray1));
   r64* rombergArray2 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray2));

   auto integrandRomberg = [n, m, p, q, rombergArray1, maxSteps, tolerance, wellWidth](r64 x1)
   {

      auto inner = [n, m, p, q, x1](r64 x2)
      {
         r64 value = 0.0;

         if(!floatCompare(x1, x2))
         {
            value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/(x1 - x2);
         }

         return value;
      };

      r64 value = rombergIntegrate(inner, 0.0, x1, rombergArray1, maxSteps, tolerance);

      return value;
   };

   auto integrandTanhSinh = [n, m, p, q, maxSteps, tolerance, wellWidth](r64 x1)
   {
      auto inner = [n, m, p, q, x1](r64 x2)
      {
         r64 value = 0.0;

         if(!floatCompare(x1, x2))
         {
            value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/(x1 - x2);
         }

         return value;
      };

      r64 value = tanhSinhQuadrature(inner, 0.0, x1, tolerance);

      return value;
   };

   r64 rombergValue = 2.0*rombergIntegrate(integrandRomberg, 0, wellWidth, rombergArray2, maxSteps, tolerance);
   r64 tanhSinhValue = 2.0*tanhSinhQuadrature(integrandTanhSinh, 0.0, wellWidth, tolerance);
   r64 gaussLegendreValue = 2.0*twoElectronGaussLegendre(i, j, k, l, iWellWidth, 0, wellWidth, 100);
   printf("i = %d, j = %d, k = %d, l = %d\n\n", i, j, k, l);
   printf("Romberg: %.10e\n", rombergValue);
   printf("tanh-sinh: %.10e\n", tanhSinhValue);
   printf("GL: %.10e\n", gaussLegendreValue);

   if(fabs(rombergValue - tanhSinhValue) > 0.001*rombergValue && rombergValue > 1e-10)
   {
      printf("Warning: integral values differ significantly(i=%zd, j=%zd, k=%zd, l=%zd): \n", i, j, k, l);
      printf("\t\tRomberg: %e, tanh-sinh: %e, GL: %e\n", rombergValue, tanhSinhValue, gaussLegendreValue);
   }

}
