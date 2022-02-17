/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "TaggedInterval.h"
#include "BasicFunctions.h"
#include "SmoothStepFunction.h"
#include "IndicatorFunction.h"
#include "BumpFunction.h"

//
// TODO: make interval computers for the following function types
//		quintic_poly_step
//		constant_fun
//		zero_fun
//		mollified_indicator_fun

namespace CTR
{
	namespace Functions
	{
		template <typename T>
		struct interval_computer;

		template <typename T>
		struct interval_computer < zero_fun<T> >
		{
			static void call(const zero_fun<T> &f, TInterval::IntervalList& il)
			{
			}
		};

		template <typename T>
		struct interval_computer < constant_fun<T> >
		{
			static void call(const constant_fun<T>& f, TInterval::IntervalList& il)
			{
			}
		};

		template <typename T>
		struct interval_computer < smooth_step_function<T> >
		{
			static void call(const smooth_step_function<T> &f, TInterval::IntervalList &il)
			{
				il.push_back(TInterval::TaggedInterval(f.GetSw() - f.GetWidth() / 2.0, f.GetSw() + f.GetWidth() / 2.0, TInterval::DENSE));
			}
		};

		template <typename T>
		struct interval_computer < smooth_step_fun_der<T> >
		{
			static void call(const smooth_step_fun_der<T> &f, TInterval::IntervalList &il)
			{
				il.push_back(TInterval::TaggedInterval(f.GetSw() - f.GetWidth() / 2.0, f.GetSw() + f.GetWidth() / 2.0, TInterval::DENSE));
			}
		};

		template <typename T>
		struct interval_computer < mollified_indicator_fun<T> >
		{
			static void call(const mollified_indicator_fun<T> &f, TInterval::IntervalList &il)
			{
				const smooth_step_function<T> &fLeft = f.GetLeftStep();
				const smooth_step_function<T> &fRight = f.GetRightStep();
	
				il.push_back(TInterval::TaggedInterval(fLeft.GetSw() - fLeft.GetWidth() / 2.0,
					fLeft.GetSw() + fLeft.GetWidth() / 2.0,
					TInterval::DENSE));
				il.push_back(TInterval::TaggedInterval(fRight.GetSw() - fRight.GetWidth() / 2.0,
					fRight.GetSw() + fRight.GetWidth() / 2.0,
					TInterval::DENSE));
			}
		};

		template <typename T>
		struct interval_computer < bump_function<T> >
		{

			static void call(const bump_function<T> &f, TInterval::IntervalList &il)
			{
				const smooth_step_function<T> &fLeft = f.GetLeftStep();
				const smooth_step_function<T> &fRight = f.GetRightStep();

				il.push_back(TInterval::TaggedInterval(fLeft.GetSw() - fLeft.GetWidth() / 2.0,
					fLeft.GetSw() + fLeft.GetWidth() / 2.0,
					TInterval::DENSE));
				il.push_back(TInterval::TaggedInterval(fRight.GetSw() - fRight.GetWidth() / 2.0,
					fRight.GetSw() + fRight.GetWidth() / 2.0,
					TInterval::DENSE));
			}
		};

		template <>
		struct interval_computer < indicator_function >
		{
			static void call(const indicator_function &f, TInterval::IntervalList &il)
			{
			}
		};


		template <typename F>
		void add_intervals(const F& f, TInterval::IntervalList &il)
		{
			interval_computer<F>::call(f, il);
		}
	}

}
