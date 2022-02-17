/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "BasicFunctions.h"
#include "StandardLinearScale.h"

#ifndef STRONG_INLINE
#define STRONG_INLINE
#endif

namespace CTR
{
	namespace Functions
	{
		template <typename R>
		class quintic_poly_step;

		template <typename R>
		class quintic_poly_step_der;
	}

	namespace internal
	{
		struct poly_step_der_maker;
	}

	namespace internal
	{
		template <typename R>
		struct traits;

		//suport automatic differentiation
		template <typename R>
		struct traits < Functions::quintic_poly_step<R> >
		{
			typedef poly_step_der_maker	derivative_maker;
		};
	}

	namespace Functions
	{
		/*!
			\brief  A smooth step function with zero derivative at both ends

			This class implements a smoothed step function,
			[0,1]->[left,right], where the left and right values are passed
			in by the user
			*/
		template <typename R>
		class quintic_poly_step
		{
		public:
			typedef R result_type;
			quintic_poly_step() {}

			/*!
				\brief	Constructor which takes left and right values

				The constructor results in a quintic polynomial step which
				takes the domain [0,1] and maps it smoothly to the values in
				the range [left,right], with zero slope at the ends of the domain.
				*/
			quintic_poly_step(R const& left, R const& right) :
				m_scale(left, right)
			{
			}

			/*!
				\brief	Evaluate the polynomial step at s

				\param s input value
				\return The value of the step function at s
				*/
         STRONG_INLINE
			result_type operator()(double s) const
			{
#ifdef COUNT_OPS
				g_nadd += 1;
				g_nmult += 4;
#endif
				double s2 = s*s;
				return m_scale(s2*(3.0-2.0*s));
				//return m_scale(3.0*s*s - 2.0*s*s*s);
			}

			//! \brief Get the output value at the left endpoint
			const result_type & GetLeft() const
			{
				return m_scale.GetLeft();
			}

			//! \brief Get the output value at the right endpoint
			const result_type & GetRight() const
			{
				return m_scale.GetRight();
			}

			//! \brief Get the scaling function
			const standard_linear_scale<R> & GetScale() const
			{
				return m_scale;
			}

			quintic_poly_step<R>& operator=(const quintic_poly_step<R>& other)
			{
				m_scale = other.m_scale;
				return *this;
			}

		private:
			standard_linear_scale<R>	m_scale;
		};

		/*!
			\brief  Derivative of the quintic_poly_step class

			This class implements the derivative of a quintic_poly_step
			function.
			*/
		template <typename R>
		class quintic_poly_step_der
		{
		public:
			typedef R result_type;
         quintic_poly_step_der()
         {
         }
			/*!
				\brief Constructor taking left and right values
				*/
			quintic_poly_step_der(result_type const& left, result_type const& right) :
				m_scale_der(left, right)
			{}

			/*!
				\brief Constructor taking a quintic_poly_step object
				*/
			quintic_poly_step_der(quintic_poly_step<R> const& step_fun) :
				m_scale_der(step_fun.GetScale())
			{}

			/*!
				\brief Evaluates the derivative of the quintic_poly_step

				\param s input value
				\return  The value of the derivative
				*/
         STRONG_INLINE
			result_type operator()(double s) const
			{
				//Use the chain rule to evaluate the derivative
#ifdef COUNT_OPS
				g_nmult += 2 + numel<R>::N;
				g_nadd += 1;
#endif
				return m_scale_der(s)*(6.0*s*(1.0 - s));
			}

			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		private:
			const standard_linear_scale_der<R>	m_scale_der;
		};
	}

	namespace internal
	{
		struct poly_step_der_maker
		{
			template <typename T> struct result;

			template <typename T>
			struct result < poly_step_der_maker(Functions::quintic_poly_step<T> const&) >
			{
				typedef Functions::quintic_poly_step_der<T>	type;
			};

			template <typename T>
			Functions::quintic_poly_step_der<T>	operator()(Functions::quintic_poly_step<T> const &ps)
			{
				return Functions::quintic_poly_step_der<T>(ps);
			}

		};
	}

	namespace Functions
	{

		//Forward declarations
		template <typename R> class smooth_step_function;
		template <typename R> class smooth_step_fun_der;
	}

	namespace internal
	{
		struct smooth_step_der_maker;
	}

	namespace internal
	{
		template <typename T> struct traits;

		//The traits class for a function structure
		// needs to define the derivative_type for the
		// automatic derivative taking routine
		// get_derivative
		template <typename R>
		struct traits < Functions::smooth_step_function<R> >
		{
			typedef smooth_step_der_maker	derivative_maker;
		};
	}

	namespace Functions
	{
		/*!
			\brief A smoothed "step" function

			This class implements a smoothed "step" function, with user definable parameters
			such as the values at the left and right of the step, the "width" of the step, and
			the input value at which the center of the smoothed step should be.

			This function supports get_derivative.

			Example code:

			CTR::smooth_step_function<double> step_fun(0.0,10.0,5.0,1.0);
			double a = step_fun(0.0);   //a = 0.0
			double b = step_fun(5.0);   //b = 5.0
			double c = step_fun(10.0);  //c = 10.0
			*/
		template <typename R>
		class smooth_step_function
		{
		public:
			typedef R result_type;

			smooth_step_function() {}

			/*!
				\brief Constructor taking all user-definable parameters

				\param[in] left	The value to the left of the step
				\param[in] right The value to the right of the step
				\param[in] sw	The location in the domain where the center of the smoothed step should be
				\param[in] width The "width" of the smoothed step, so that [sw-width/2, sw+width/2] contains the smoothed step.
				*/
			smooth_step_function(result_type const& left, result_type const& right, double sw, double width) :
				m_left(left), m_right(right), m_sw(sw), m_halfwidth(width/2.0), m_poly(left, right), m_invWidth(1.0/width)
			{
			}

			/*!
				\brief Evaluates the smoothed step function

				\param[in] s The input value
				\return The appropriate value, according to the portion of the domain where the input is given
				If the input is to the left of the step, the value left is returned; if it is to the right of
				the step, the value right is returned. Otherwise, the smoothed quintic polynomial step is computed
				using quintic_poly_step
				*/
         STRONG_INLINE
			result_type operator()(double s) const
			{
				if (s < m_sw - m_halfwidth)
					return m_left;
				else if (s > m_sw + m_halfwidth)
					return m_right;
				else {
#ifdef COUNT_OPS
					g_nadd += 2;
					g_nmult += 1;
#endif
					return m_poly((s - m_sw + m_halfwidth) * m_invWidth);
			}
			}

			//! \brief Return the left value
			const R & GetLeft() const { return m_left; }
			//! \brief Return the right value
			const R & GetRight() const { return m_right; }
			//! \brief Return the location where the step occurs
			double GetSw() const { return m_sw; }
			//! \brief Return the Width of the smoothed step
			double GetWidth() const { return 2.0*m_halfwidth; }
			//! \brief Return the Half Width of the smoothed step
			double GetHalfWidth() const { return m_halfwidth; }
			//! \brief Return the quintic_poly_step<R> object which computes the smoothed step
			double GetInvWidth() const { return m_invWidth; }
			//! \brief Return the quintic_poly_step<R> object which computes the smoothed step
			const quintic_poly_step<R> & GetPoly() const { return m_poly; }

			smooth_step_function<R>& operator=(const smooth_step_function<R>& other)
			{
				m_left = other.m_left;
				m_right = other.m_right;
				m_sw = other.m_sw;
				m_halfwidth = other.m_halfwidth;
				m_invWidth = other.m_invWidth;
				m_poly = other.m_poly;
				return *this;
			}

		private:
			R m_left;
			R m_right;
			double m_sw;
			double m_halfwidth;
			double m_invWidth;
			quintic_poly_step<R>	m_poly;
		};

		/*!
			\brief Computes the derivative of a smooth_step_function

			This class implements the derivative of a smoothed "step" function.
			*/
		template <typename R>
		class smooth_step_fun_der
		{
		public:
			typedef R result_type;

         smooth_step_fun_der()
         {
         }

			/*!
				\brief Constructor taking a smooth_step_function<R> argument
				\param[in] step_fun	a smooth_step_function<R> object for which the derivative is needed
				*/
			smooth_step_fun_der(const smooth_step_function<R> &step_fun) :
				m_sw(step_fun.GetSw()), m_halfwidth(step_fun.GetHalfWidth()), m_invWidth(step_fun.GetInvWidth()),
				m_poly_der(get_derivative(step_fun.GetPoly())) //invoke copy constructor here
			{
			}

			/*!
				\brief Evaluate the derivative of the smooth step function
				\param[in] s The input value
				\return If the input is to the left or right of the step (i.e. less than of sw-width/2 or greater than sw+width/2),
				the return value is zero. Otherwise the value returned is computed by quintic_poly_step_der and the chain rule for scaling.

				This operator computes the derivative of the smoothed step function which was passed into
				the constructor at the value s.
				*/
         STRONG_INLINE
			result_type operator()(double s) const
			{
				if (s < m_sw - m_halfwidth)
					return zero_fun<result_type>()();
				else if (s > m_sw + m_halfwidth)
					return zero_fun<result_type>()();
				else {
#ifdef COUNT_OPS
					g_nmult += 2;
					g_nadd += 2;
#endif
					return m_invWidth*m_poly_der((s - m_sw + m_halfwidth) * m_invWidth ); //chain rule on smooth_step_fun formula
				}
			}

			double GetSw() const { return m_sw; }
			double GetWidth() const { return 2.0*m_halfwidth; }

			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		private:
			const double m_sw;
			const double m_halfwidth;
			const double m_invWidth;
			const quintic_poly_step_der<R>	m_poly_der;
		};
	}

	namespace internal
	{
		struct smooth_step_der_maker
		{
			template <typename T> struct result;

			template <typename T>
			struct result < smooth_step_der_maker( Functions::smooth_step_function<T> const&) >
			{
				typedef Functions::smooth_step_fun_der<T>	type;
			};

			template <typename T>
			Functions::smooth_step_fun_der<T>	operator()(Functions::smooth_step_function<T> const &stepfun)
			{
				return Functions::smooth_step_fun_der<T>(stepfun);
			}

		};
	}

}
