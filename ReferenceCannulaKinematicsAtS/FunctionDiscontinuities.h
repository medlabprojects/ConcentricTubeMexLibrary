#pragma once

#include <boost/container/static_vector.hpp>
#include "BasicFunctions.h"
#include "BumpFunction.h"
#include "IndicatorFunction.h"
#include "SmoothStepFunction.h"

#define MAX_DISCONTINUITIES	50

namespace CTR
{
	namespace Functions
	{
		template <typename F>
		struct discontinuity_lister;

		template <>
		struct discontinuity_lister < indicator_function >
		{
			template <typename Cont>
			static void call(indicator_function const& f, Cont &pts)
			{
				pts.push_back(f.left);
				pts.push_back(f.right);
			}
		};

		template <typename T>
		struct discontinuity_lister < step_function<T> >
		{
			template <typename Cont>
			static void call(step_function<T> const& f, Cont &pts)
			{
				pts.push_back(f.sw);
			}
		};

		template <typename T>
		struct discontinuity_lister < constant_fun<T> >
		{
			template <typename Cont>
			static void call(constant_fun<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < zero_fun<T> >
		{
			template <typename Cont>
			static void call(zero_fun<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < smooth_step_function<T> >
		{
			template <typename Cont>
			static void call(smooth_step_function<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < mollified_indicator_fun<T> >
		{
			template <typename Cont>
			static void call(mollified_indicator_fun<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < quintic_poly_step<T> >
		{
			template <typename Cont>
			static void call(quintic_poly_step<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < quintic_poly_step_der<T> >
		{
			template <typename Cont>
			static void call(quintic_poly_step_der<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < bump_function<T> >
		{
			template <typename Cont>
			static void call(bump_function<T> const &f, Cont &pts)
			{
			}
		};

		template <typename T>
		struct discontinuity_lister < bump_function_der<T> >
		{
			template <typename Cont>
			static void call(bump_function_der<T> const &f, Cont &pts)
			{
			}
		};


		template <typename F, typename Cont>
		void add_discontinuities(F const &f, Cont &pts)
		{
			discontinuity_lister<F>::template call<Cont>(f, pts);
		}
	}




}