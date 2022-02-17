/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <utility>

#ifndef STRONG_INLINE
#define STRONG_INLINE
#endif

namespace CTR
{
	namespace internal
	{
		struct plus
		{
			template <typename T>
			struct result;

			template <typename F, typename Left, typename Right>
			struct result < F(Left, Right) >
			{
				typedef decltype(std::declval<Left>() + std::declval<Right>())	type;
			};

			template <typename Left, typename Right>
			typename result<plus(Left, Right)>::type operator()(Left&& l, Right&& r) const
			{
#ifdef COUNT_OPS
            g_nadd += std::decay<Left>::type::RowsAtCompileTime * std::decay<Left>::type::ColsAtCompileTime;
#endif
				return l + r;
			}


		};

		struct minus
		{
			template <typename T>
			struct result;

			template <typename F, typename Left, typename Right>
			struct result < F(Left, Right) >
			{
				typedef decltype(std::declval<Left>() - std::declval<Right>())	type;
			};

			template <typename Left, typename Right>
			typename result<minus(Left, Right)>::type operator()(Left&& l, Right&& r) const
			{
#ifdef COUNT_OPS
            g_nadd += std::decay<Left>::type::RowsAtCompileTime * std::decay<Left>::type::ColsAtCompileTime;
#endif
				return l - r;
			}
		};

		struct times
		{
			template <typename T>
			struct result;

			template <typename F, typename Vector>
			struct result < F(Vector, double) >
			{
            typedef decltype(std::declval<double>() * std::declval<Vector>())	type;
			};

			template <typename F, typename Vector>
			struct result < F(double, Vector) >
			{
				typedef decltype( std::declval<double>() * std::declval<Vector>() )	type;
			};

			template <typename F>
			struct result < F(double, double) >
			{
				typedef double type;
			};

			result<times(double, double)>::type	operator()(double a, double b) const
			{
#ifdef COUNT_OPS
            g_nmult += 1;
#endif
				return a*b;
			}

			template <typename Vector>
			typename result<times(Vector, double)>::type operator()(Vector&& v, double a) const
			{
#ifdef COUNT_OPS
            g_nmult += std::decay<Vector>::type::RowsAtCompileTime;
#endif
				return a*v;
			}

			template <typename Vector>
			typename result<times(double, Vector)>::type operator()(double a, Vector&& v) const
			{
#ifdef COUNT_OPS
            g_nmult += std::decay<Vector>::type::RowsAtCompileTime;
#endif
				return a*v;
			}
		};
	}
}