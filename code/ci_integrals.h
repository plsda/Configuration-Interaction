#include <stdio.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


// Returns 2^(stage-1)-point refinement to the result of 2^(stage-1)-subinterval trapezoidal rule.
// Together with previous stages, the result can be used to compute trapezoidal rule with 2^stage subintervals.
template <typename f>
r64 trapezoidalRuleStage(f& integrand, r64 lowerBound, r64 upperBound, size_t stage)
{
   //assert(lowerBound < upperBound);

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

//A rudimentary implementation of tanh-sinh quadrature for one-parameter integrand
template <typename f>
r64 tanhSinhQuadrature(f& integrand, r64 lowerBound, r64 upperBound, r64 absTolerance) 
{
   static r64* myAbscissae;
   static r64* myWeights;
   static size_t storedNonnegativePointCount;
   size_t maxStages = 8; 
   const r64 windowSize = 1.0;

   if(!storedNonnegativePointCount)
   {
      size_t maxNonnegativePoints = 2000;
      r64 minSubintervalLength = windowSize/(1ULL << maxStages);
      myAbscissae = (r64*)malloc(maxNonnegativePoints*sizeof(*myAbscissae));
      myWeights = (r64*)malloc(maxNonnegativePoints*sizeof(*myWeights));

      for(storedNonnegativePointCount = 0; 
          storedNonnegativePointCount < maxNonnegativePoints; 
          storedNonnegativePointCount++)
      {
         r64 pointCosh = cosh(HPI64*sinh(storedNonnegativePointCount*minSubintervalLength));
         myAbscissae[storedNonnegativePointCount] = tanh(HPI64*sinh(storedNonnegativePointCount*minSubintervalLength));
         myWeights[storedNonnegativePointCount] = HPI64*cosh(storedNonnegativePointCount*minSubintervalLength)/(pointCosh*pointCosh);

         if(myWeights[storedNonnegativePointCount] < 1.0e-13)
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
   result = subintervalLength*myWeights[0]*transformedIntegrand(0);

   //printf("stage 1\n");
   //Evaluate the integrand at every first-stage abscissa
   size_t nodeIndex = stride;
   for(int i = 1; 
       nodeIndex < storedNonnegativePointCount; 
       i++, nodeIndex = i*stride)
   {
      r64 abscissa = myAbscissae[nodeIndex]; 
      result += myWeights[nodeIndex]*(transformedIntegrand(abscissa) + transformedIntegrand(-abscissa));

      //printf("\t%d: abscissa = %f, weight = %f\n", i, abscissa, myWeights[nodeIndex]);
   }

   stride >>= 1;
   for(stage = 2; stage <= maxStages; stage++)
   {
      //printf("stage %d\n", (int)stage);

      subintervalLength *= 0.5;
      //subintervalLength = windowSize / (1 << (stage-1));

      r64 stageSum = 0;
      nodeIndex = stride;
      //Now we only need the abscissas and weights with odd indices
      for(int i = 0; 
          nodeIndex < storedNonnegativePointCount;
          i++, nodeIndex = (2*i + 1)*stride)
      {
         r64 abscissa = myAbscissae[nodeIndex]; 
         stageSum += myWeights[nodeIndex]*(transformedIntegrand(abscissa) + transformedIntegrand(-abscissa));

         //printf("\t%d: abscissa = %f, weight = %f\n", i, abscissa, myWeights[nodeIndex]);
      }

      if(subintervalLength*fabs(stageSum) < absTolerance && stage > 1)
      {
         stage--;
         break;
      }

      result += stageSum;
      stride >>= 1;

      //printf("\ttanh-sinh 1: %.15f\n\n", result);
   }

   //printf("\n ----- \n");

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

   r64 result = gsl_integration_glfixed(&F, paramsIn->lowerBound, paramsIn->upperBound, paramsIn->glTable);

	return result;
}

//Same bounds for both integrals
r64 twoElectronGaussLegendre(int i, int j, int k, int l, r64 iWellWidth, r64 lowerBound, 
                             r64 upperBound, size_t maxSteps = 32)
{
   gsl_integration_glfixed_table* glTable = gsl_integration_glfixed_table_alloc(maxSteps);

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
   size_t integralCount;
   size_t maxSpatialIndex;
   size_t maxSteps;
   u32 method;
   r64 wellWidth;
   r64 tolerance;
};
struct TwoElectronIntegral
{
   //Could later infer these indices, but good to have them for now
   size_t i, j, k, l;
   r64 value;
};
#pragma pack(pop)

int writeTwoElectronIntegrals(size_t maxSpatialIndex, r64 wellWidth, r64 tolerance, TwoElectronIntegral* storage, size_t maxIntegrals)
{
   r64 iWellWidth = 1.0/wellWidth;
   size_t maxSteps = 20;//50;

   r64* rombergArray1 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray1));
   r64* rombergArray2 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray2));

   //FILE* outputFile = fopen("integrals.dat", "wb");

   //if(outputFile)
   //{
   if(storage)
   {

   //for(int i = 1; i < 4; i++)
   //{
   //   for(int j = 1; j < i+2; j++)
   //   {
   //      for(int k = 1; k < i+2; k++)
   //      {
   //         for(int l = 1; l < i+2; l++)
   //         {
   //Only calculate unique integrals(i > j etc.)
      int integralCount = 0;
      for(int i = 1; i < maxSpatialIndex && integralCount < maxIntegrals; i++)
      {
         for(int j = i; j < maxSpatialIndex && integralCount < maxIntegrals; j++)
         {
            for(int k = i; k < maxSpatialIndex && integralCount < maxIntegrals; k++)
            {
               for(int l = k; l < maxSpatialIndex && integralCount < maxIntegrals; l++)
               {

                  r64 n = PI64*iWellWidth*i; 
                  r64 m = PI64*iWellWidth*j; 
                  r64 p = PI64*iWellWidth*k; 
                  r64 q = PI64*iWellWidth*l;

                  size_t rombergEvals = 0;
                  size_t tanhSinhEvals = 0;

                  auto integrandRomberg = [n, m, p, q, rombergArray1, maxSteps, tolerance, wellWidth, &rombergEvals](r64 x1)
                  {

                     auto inner = [n, m, p, q, x1](r64 x2)
                     {
                        r64 value = 0.0;

                        if(!floatCompare(x1, x2))
                        {
                           value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/(x1 - x2);
                           //value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/fabs(x1 - x2);
                        }

                        return value;
                     };

                     r64 value = 2.0*rombergIntegrate(inner, 0.0, x1, rombergArray1, maxSteps, tolerance);

                     return value;
                  };

                  //auto integrandTanhSinh = [n, m, p, q, maxSteps, tolerance, wellWidth](r64 x1)
                  //{

                  //   auto inner = [n, m, p, q, x1](r64 x2)
                  //   {
                  //      r64 value = 0.0;

                  //      if(!floatCompare(x1, x2))
                  //      {
                  //         value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/(x1 - x2);
                  //         //value = (sin(n*x1)*sin(m*x2) - sin(n*x2)*sin(m*x1))*(sin(p*x1)*sin(q*x2) - sin(p*x2)*sin(q*x1))/fabs(x1 - x2);
                  //      }

                  //      return value;
                  //   };

                  //   //r64 value = 2.0*tanhSinhQuadrature(inner, 0.0, wellWidth, tolerance);

                  //   return value;
                  //};

                  r64 rombergValue = rombergIntegrate(integrandRomberg, 0, wellWidth, rombergArray2, maxSteps, tolerance, rombergEvals);
                  //r64 tanhSinhValue = tanhSinhQuadrature(integrandTanhSinh, 0.0, wellWidth, 10e-14);
                  //r64 integralValue = 
                  //   (floatCompare(rombergValue, tanhSinhValue) || rombergValue < 10.0e-12)? rombergValue : (tanhSinhValue, printf("\n\tRomberg discarded\n"));
                  //printf("Romberg: %e, tanh-sinh: %e\n", rombergValue, tanhSinhValue);
                  //
                  TwoElectronIntegral& integral = storage[integralCount++];
                  integral.i = i;
                  integral.j = j;
                  integral.k = k;
                  integral.l = l;
                  integral.value = rombergValue;

                  //TwoElectronIntegral integral;
                  //integral.i = i;
                  //integral.j = j;
                  //integral.k = k;
                  //integral.l = l;
                  //integral.value = rombergValue;
                  //fwrite(&integral, sizeof(integral), 1, outputFile);

                  //printf("i = %d, j = %d, k = %d, l = %d\n\t%.10e, %.10e\n", i, j, k, l, rombergValue, tanhSinhValue);
                  //r64 gaussLegendreValue = twoElectronGaussLegendre(i, j, k, l, iWellWidth, 0, wellWidth, 5000);
                  //printf("i = %d, j = %d, k = %d, l = %d\n\t%.10e, %.10e, %.10e\n", i, j, k, l, rombergValue, tanhSinhValue, gaussLegendreValue);
                  
                  //printf("%d - i = %d, j = %d, k = %d, l = %d\n", runningIndex++, i, j, k, l);
               }
            }
         }
      }

      free(rombergArray1);
      free(rombergArray2);
      //fclose(outputFile);

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

r64 getTwoElectronIntegral(size_t i, size_t j, size_t k, size_t l, r64 wellWidth, r64* scratch, size_t maxSteps)
{
   r64 iWellWidth = 1.0/wellWidth;
   r64 tolerance = 1e-8;
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

      //r64 value = 2.0*rombergIntegrate(inner, 0.0, x1, rombergArray1, maxSteps, tolerance, rombergEvals, 2);
      //r64 value = 2.0*adaptiveSimpson(inner, 0.0, x1, tolerance, maxSteps);
      r64 value = 2.0*tanhSinhQuadrature(inner, 0.0, x1, tolerance);

      return value;
   };

   //value = rombergIntegrate(integrand, 0, wellWidth, rombergArray2, maxSteps, tolerance, rombergEvals);
   //value = adaptiveSimpson(integrand, 0.0, wellWidth, tolerance, maxSteps);
   value = tanhSinhQuadrature(integrand, 0.0, wellWidth, tolerance);

   return value;
}


void testIntegration(int i, int j, int k, int l, r64 wellWidth)
{
   r64 iWellWidth = 1.0/wellWidth;
   r64 n = PI64*iWellWidth*i; 
   r64 m = PI64*iWellWidth*j; 
   r64 p = PI64*iWellWidth*k; 
   r64 q = PI64*iWellWidth*l;
   size_t rombergEvals = 0;
   size_t tanhSinhEvals = 0;

   r64 tolerance = 1e-8;
   size_t maxSteps = 50;
   r64* rombergArray1 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray1));
   r64* rombergArray2 = (r64*)malloc(2*maxSteps*sizeof(*rombergArray2));

   auto integrandRomberg = [n, m, p, q, rombergArray1, maxSteps, tolerance, wellWidth, &rombergEvals](r64 x1)
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

      r64 value = 2.0*rombergIntegrate(inner, 0.0, x1, rombergArray1, maxSteps, tolerance);

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

      r64 value = 2.0*tanhSinhQuadrature(inner, 0.0, wellWidth, tolerance);

      return value;
   };

   r64 rombergValue = rombergIntegrate(integrandRomberg, 0, wellWidth, rombergArray2, maxSteps, tolerance, rombergEvals);
   r64 tanhSinhValue = tanhSinhQuadrature(integrandTanhSinh, 0.0, wellWidth, tolerance);//10e-14);
   r64 gaussLegendreValue = twoElectronGaussLegendre(i, j, k, l, iWellWidth, 0, wellWidth, 5000);
   printf("i = %d, j = %d, k = %d, l = %d\n\n", i, j, k, l);
   printf("Romberg: %e\n", rombergValue);
   printf("tanh-sinh: %e\n", tanhSinhValue);
   printf("GL: %e\n", gaussLegendreValue);
}
