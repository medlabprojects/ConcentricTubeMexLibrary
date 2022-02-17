/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "SmoothStepFunction.h"
#include "BasicFunctions.h"

#ifndef STRONG_INLINE
#define STRONG_INLINE
#endif

namespace CTR
{
	namespace Functions
	{
		template <typename T> struct mollified_indicator_fun;
		template <typename T> struct mollified_indicator_fun_der;
	}

	namespace internal
	{
		struct mollified_indicator_der_maker;
	}

	namespace internal
	{
		template <typename T> struct traits;

		template <typename T>
		struct traits < Functions::mollified_indicator_fun<T> >
		{
			typedef mollified_indicator_der_maker	derivative_maker;
		};
	}

	namespace Functions
	{
		template <typename T>
		struct mollified_indicator_fun
		{

			mollified_indicator_fun() {}

			mollified_indicator_fun(double left, double right, double width) :
				m_stepLeft(static_cast<T>(0.0), static_cast<T>(1.0), left + width / 2.0, width),
				m_stepRight(static_cast<T>(0.0), static_cast<T>(-1.0), right - width / 2.0, width)
			{
			}

         STRONG_INLINE
			T operator()(double s) const
			{
#ifdef COUNT_OPS
				g_nadd += numel<T>::N;
#endif
				return m_stepLeft(s) + m_stepRight(s);
			}

			const smooth_step_function<T> &GetLeftStep() const
			{
				return m_stepLeft;
			}

			const smooth_step_function<T> &GetRightStep() const
			{
				return m_stepRight;
			}

			mollified_indicator_fun<T>& operator=(const mollified_indicator_fun<T>& other)
			{
				m_stepLeft = other.m_stepLeft;
				m_stepRight = other.m_stepRight;
				return *this;
			}

		private:
			smooth_step_function<T>	m_stepLeft;
			smooth_step_function<T> m_stepRight;
		};

		template <typename T>
		struct mollified_indicator_fun_der
		{
         mollified_indicator_fun_der()
         {
         }

			mollified_indicator_fun_der(const mollified_indicator_fun<T> &f) :
				m_derLeft(f.GetLeftStep()),
				m_derRight(f.GetRightStep())
			{
			}

         STRONG_INLINE
			T operator()(double s) const
			{
#ifdef COUNT_OPS
				g_nadd += numel<T>::N;
#endif
				return m_derLeft(s) + m_derRight(s);
			}

		private:
			smooth_step_fun_der<T> m_derLeft;
			smooth_step_fun_der<T> m_derRight;
		};
	}

	namespace internal
	{
		struct mollified_indicator_der_maker
		{
			template <typename T> struct result;

			template <typename T>
			struct result < mollified_indicator_der_maker(Functions::mollified_indicator_fun<T> const&) >
			{
				typedef Functions::mollified_indicator_fun_der<T>	type;
			};

			template <typename T>
			typename result< mollified_indicator_der_maker(Functions::mollified_indicator_fun<T> const&) >::type
				operator()(Functions::mollified_indicator_fun<T> const &ifun) const
			{
				return Functions::mollified_indicator_fun_der<T>(ifun);
			}
		};
	}
}
