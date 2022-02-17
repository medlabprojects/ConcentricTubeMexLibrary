/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <Eigen/Dense>

//If this is the microsoft compiler, use __declspec(align(16)) for memory-aligned objects
#ifdef _MSC_VER
#define ALIGN_SPEC	__declspec(align(16))
#elif defined(__GLIBC__)
#define ALIGN_SPEC __attribute__((aligned(16)))
#elif defined(__clang__)
#define ALIGN_SPEC __attribute__((aligned(16)))
#endif

namespace CTR
{

	typedef Eigen::Matrix3d	RotationMatrix;

	template <int N>
	struct Vector
	{
		typedef Eigen::Matrix<double, N, 1> type;
	};

	template <int M, int N>
	struct Matrix
	{
		typedef Eigen::Matrix<double, M, N> type;
	};
}
