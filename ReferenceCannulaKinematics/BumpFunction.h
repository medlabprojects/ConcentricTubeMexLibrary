/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "SmoothStepFunction.h"

namespace CTR
{
	namespace Functions
	{
		template <typename T> struct bump_function;

		template <typename T> struct bump_function_der;
	}

	namespace internal
	{
		struct bump_function_der_maker;
	}

	namespace internal
	{
		template <typename T> struct traits;

		template <typename T>
		struct traits < Functions::bump_function<T> >
		{
			typedef bump_function_der_maker	derivative_maker;
		};
	}

	namespace Functions
	{
		template <typename T>
		struct bump_function
		{
			bump_function(T location, T width) :
				m_leftStep(0.0, 1.0, location - width / 4.0, width / 2.0),
				m_rightStep(0.0, -1.0, location + width / 4.0, width / 2.0),
				m_A(1.0 / 0.5 / width)
			{}

			T operator()(double s) const
			{
				return m_A*(m_leftStep(s) + m_rightStep(s));
			}

			const smooth_step_function<T>& GetLeftStep() const { return m_leftStep; }
			const smooth_step_function<T>& GetRightStep() const { return m_rightStep; }

		private:
			smooth_step_function<T>	m_leftStep;
			smooth_step_function<T> m_rightStep;
			double m_A;

			friend struct bump_function_der < T > ;
		};

		template <typename T>
		struct bump_function_der
		{
			bump_function_der(const bump_function<T>& f) :
				m_leftStepDer(f.m_leftStep),
				m_rightStepDer(f.m_rightStep),
				m_A(f.m_A)
			{
			}

			T operator()(double s) const
			{
				return m_A*(m_leftStepDer(s) + m_rightStepDer(s));
			}

		private:
			smooth_step_fun_der<T> m_leftStepDer;
			smooth_step_fun_der<T> m_rightStepDer;
			double m_A;
		};
	}

	namespace internal
	{
		struct bump_function_der_maker
		{
			template <typename T>
			struct result;

			template <typename T>
			struct result < bump_function_der_maker(T const&) >
			{
				typedef Functions::bump_function_der<T>	type;
			};

			template <typename T>
			Functions::bump_function_der<T>	operator()(Functions::bump_function<T> const& bf)
			{
				return Functions::bump_function_der<T>(bf);
			}

		};
	}


}