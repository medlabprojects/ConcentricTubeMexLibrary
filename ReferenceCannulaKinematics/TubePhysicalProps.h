/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

namespace CTR
{
	inline double AnnularSecondAreaMoment(double OD, double ID)
	{
		return M_PI / 64.0*(pow(OD, 4) - pow(ID, 4));
	}

	inline double AnnularBendingStiffness(double OD, double ID, double E)
	{
		return AnnularSecondAreaMoment(OD, ID)*E;
	}

	inline double AnnularTorsionalCompliance(double OD, double ID, double G)
	{
		return 1.0/(2.0*AnnularSecondAreaMoment(OD, ID)*G);
	}
}