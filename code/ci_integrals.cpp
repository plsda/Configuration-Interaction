

void integralTest()
{
	r64 i = 7;
	r64 j = 1;
	r64 k = 3;
	r64 l = 1;
	r64 integralParams[] = 
	{
		0,
		0.5*i*PI64/wellWidth, 0.5*j*PI64/wellWidth,
		0.5*k*PI64/wellWidth, 0.5*l*PI64/wellWidth
	}; //y(blank/reserved),n,m,p,q

	gsl_function integrand;
	integrand.params = integralParams;

	r64 epsabs = 0;
	r64 epsrel = 1e-6;
	r64 result1 = 0;
	r64 result2 = 0;
	r64 resultErrorEstimate1 = 0;
	r64 resultErrorEstimate2 = 0;

	r64 lowerBound1 = 0;
	r64 upperBound1 = wellWidth; 
	r64 lowerBound2 = -wellWidth;
	r64 upperBound2 = 0; 

	integrand.function = &twoElectronIntegrandOuter1;
	gsl_integration_qags(&integrand, lowerBound1, upperBound1, epsabs, epsrel, 1000, 
								globalWorkspace, &result1, &resultErrorEstimate1);
	integrand.function = &twoElectronIntegrandOuter2;
	gsl_integration_qags(&integrand, lowerBound2, upperBound2, epsabs, epsrel, 1000, 
								globalWorkspace, &result2, &resultErrorEstimate2);

	r64 resultTotal = 0.5*result1 - 0.5*result2;
	printf("Result: %f\nError 1:%f\nError 2:%f", resultTotal, resultErrorEstimate1, resultErrorEstimate2);
}

static gsl_integration_workspace* globalWorkspace = gsl_integration_workspace_alloc(4000);
static r64 wellWidth = 1;
static r64 inverseWellWidth = 1;

r64 twoElectronIntegrandInner(r64 x, void* params)
{
	r64 y = *(r64*)params;
	r64 n = *(((r64*)params)+1);
	r64 m = *(((r64*)params)+2);
	r64 p = *(((r64*)params)+3);
	r64 q = *(((r64*)params)+4);
	r64 value = (1.0/y)*sin(n*(x+y))*sin(m*(x-y))*sin(p*(x+y))*sin(q*(x-y));

	return value;
}

r64 twoElectronIntegrandOuter1(r64 y, void* params)
{
	gsl_function integrand;
	integrand.function = &twoElectronIntegrandInner;
	((r64*)params)[0] = y;
	integrand.params = params;

	r64 lowerBound = y;
	r64 upperBound = 2*wellWidth - y;
	r64 epsabs = 0;
	r64 epsrel = 1e-7;

	r64 value;
	r64 estimatedError;
	gsl_integration_qags(&integrand, lowerBound, upperBound, epsabs, epsrel, 1000, globalWorkspace, &value, &estimatedError);

	return value;
}

r64 twoElectronIntegrandOuter2(r64 y, void* params)
{
	gsl_function integrand;
	integrand.function = &twoElectronIntegrandInner;
	((r64*)params)[0] = y;
	integrand.params = params;

	r64 lowerBound = -y;
	r64 upperBound = 2*wellWidth + y;
	r64 epsabs = 0;
	r64 epsrel = 1e-7;

	r64 value;
	r64 estimatedError;
	gsl_integration_qags(&integrand, lowerBound, upperBound, epsabs, epsrel, 1000, globalWorkspace, &value, &estimatedError);

	return value;
}


//TODO: Find implications between the conditions for a and b for different a,b. No need to check all the conditions.
r64 getJTrigIntegral(s32 aInt, s32 bInt, r64 iWellWidth, b32 invert) //NOTE: a and b without the pi/(2L) factor
{
	r64 value = 0;
	r64 a = 0.5*PI64*inverseWellWidth*(r64)aInt;
	r64 b = 0.5*PI64*inverseWellWidth*(r64)bInt;
	//TODO: Express I_ijkl in terms of trigonometric integrals
	//TODO: Use suitable J integral and get trigonometric integral value
	
	//TODO: Integrate these expressions
	if(aInt == 0)
	{
		if(bInt == 0)
		{
			value = 2.0*(wellWidth - c)/y; //a,b = 0
		}
		else
		{
			value = (sin(b*(2.0*wellWidth - c + y)) - sin(a*(y+c)))/(y*b); // a = 0, b != 0
		}
	}
	else if(bInt == 0)
	{
		value = (sin(a*(2.0*wellWidth - c + y)) - sin(a*(y+c)))/(y*a); // a != 0, b = 0
	}
	else if((aInt + bInt) != 0 && (aInt - bInt) != 0)
	{
		value = 1.0/(a*a - b*b)*(a*() - b*());// a + b != 0, a - b != 0
	}
	else if((aInt == -bInt) || (aInt == bInt))
	{
		value = (sin(2.0*a*(2.0*wellWidth - c)) - sin(2.0*a*c))/(4.0*a*y) + (wellWidth - c)*cos(2.0*a*y)/y; // a = -b != 0 or a = b != 0
	}

	return value;
}
//TODO: Define determinant/spin-orbital struct, include indices and index parity

//TODO: Overload 2-electron version
//TODO: Build integral look-up table
r64 getSpinOrbitalIntegral(u32 leftIndex, u32 rightIndex, r64 L)
{
	assert((leftIndex != 0) && (rightIndex != 0))
	r64 value = 0;
	r64 lIndex = (r64)leftIndex;
	r64 rIndex = (r64)leftIndex;
	
	if(((leftIndex % 2 == 0) && (rightIndex % 2 == 0)) || 
		((leftIndex % 2 != 0) && (rightIndex % 2 != 0)))
	{
		if(leftIndex != rightIndex)
		{
			value = IPI64*L*(1.0/(rIndex*rIndex- lIndex*lIndex))*(lIndex*cos(PI64*lIndex)*sin(PI64*rIndex) - rIndex*cos(PI64*rIndex)*sin(PI64*lIndex));
		}
		else
		{
			value = 0.5*L - 0.25*IPI64*L*sin(2.0*PI64*lIndex)/lIndex;
		}
	}

	return value;
}
r64 getSpinOrbitalIntegral(u32 leftIndex1, u32 leftIndex2, u32 rightIndex1, u32 rightIndex2, r64 iWellWidth)
{
	assert((leftIndex1 != 0) && (rightIndex1 != 0) && (leftIndex2 != 0) && (rightIndex2 != 0))
	r64 value = 0;
	r64 lIndex1 = (r64)leftIndex1;
	r64 lIndex2 = (r64)leftIndex2;
	r64 rIndex1 = (r64)leftIndex1;
	r64 rIndex2 = (r64)leftIndex2;

	if((isEven(leftIndex1) && isEven(leftIndex2) && isEven(rightIndex1) && isEven(rightIndex2)) ||
		(!isEven(leftIndex1) && !isEven(leftIndex2) && !isEven(rightIndex1) && !isEven(rightIndex2)) ||
		(!isEven(leftIndex1) && isEven(leftIndex2) && !isEven(rightIndex1) && isEven(rightIndex2)) ||
		(isEven(leftIndex1) && !isEven(leftIndex2) && isEven(rightIndex1) && !isEven(rightIndex2)))
	{
		value = 0.125*(  getJTrigIntegral(leftIndex1-rightIndex1, leftIndex2-rightIndex2, inverseWellWidth, false) 
							- getJTrigIntegral(leftIndex1+rightIndex1, leftIndex2-rightIndex2, inverseWellWidth, false) 
							- getJTrigIntegral(leftIndex1-rightIndex1, leftIndex2+rightIndex2, inverseWellWidth, false) 
							+ getJTrigIntegral(leftIndex1+rightIndex1, leftIndex2+rightIndex2, inverseWellWidth, false) 
							+ getJTrigIntegral(leftIndex1-rightIndex1, leftIndex2-rightIndex2, inverseWellWidth, true) 
							- getJTrigIntegral(leftIndex1+rightIndex1, leftIndex2-rightIndex2, inverseWellWidth, true) 
							- getJTrigIntegral(leftIndex1-rightIndex1, leftIndex2+rightIndex2, inverseWellWidth, true) 
							+ getJTrigIntegral(leftIndex1+rightIndex1, leftIndex2+rightIndex2, inverseWellWidth, true)
						  );
	}

	return value;
}


r64 integrateOneElectron(u32 det1i, u32 det1j, u32 det2i, u32 det2j, r64 iWellWidth)
{
	r64 value = 0;
	
	//TODO: Need to align determinants. For now only assert that the indices must be smaller first
		//Need to include cases det1i == det2j && det1j == det2i (or maybe not this) and det1i == det2j || det1j == det2i 
		//(i.e. check if the determinants have at least one index in common)
	assert((det1i <= det1j) && (det2i <= det2j));
	
	if(det1i == det2i && det1j == det2j)
	{
		value = 0.5*(getSpinOrbitalIntegral(det1i,det1i, inverseWellWidth) + getSpinOrbitalIntegral(det1j,det1j, inverseWellWidth));
	}
	else if(det1i == det2i || det1i == det2j || det1j == det2i || det1j == det2j)
	{
		value = 0.5*getSpinOrbitalIntegral(det1i, det2i, inverseWellWidth); //NOTE: Requires alignment of indices
	}

	return value;
}


r64 integrateTwoElectron(u32 det1i, u32 det1j, u32 det2i, u32 det2j)
{
	r64 value = 0;
	return value;
}

