/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

namespace CTR {

	//Forward Declarations
	template <typename R>
	class standard_linear_scale;

	template <typename R>
	class standard_linear_scale_der;

	namespace internal
	{
		template <typename R>
		struct traits;

		//support automatic differentiation
		template <typename R>
		struct traits < standard_linear_scale<R> >
		{
			typedef standard_linear_scale_der<R> derivative_type;
		};
	}

	/**
	* \brief Performs the standard linear scaling
	*
	* This class performs a linear scaling of its input
	* to its output, so that the range [0,1] is mapped
	* to the the range [a,b] given in the constructor
	*
	* This class supports automatic differentiation through
	* CTR::get_derivative
	*/
	template <typename R>
	class standard_linear_scale
	{
	public:
		//! Typedef for the returned type
		typedef R result_type;

		standard_linear_scale() {}

		/*!
		\brief Constructs the standard linear scaling

		This class is a function object for the standard linear scaling
		between the input [0,1] and output [a,b], for any arbitrary type R
		which has operator*(double,R) and operator+(R,R) defined.
		*/
		standard_linear_scale(R const& a, R const& b) :
			m_a(a), m_b(b)
		{}

		/*!
		\brief Evaluates the standard linear scaling function
		\param s The input value
		\return The scaled output value with type result_type
		*/
		result_type operator()(double s) const
		{
#ifdef COUNT_OPS
			g_nmult += numel<R>::N * 2;
			g_nadd += 1 + numel<R>::N;
#endif
			return (1.0 - s)*m_a + (s)*m_b;
		}

		//! \brief Get the output value at input 0
		const result_type & GetLeft() const
		{
			return m_a;
		}

		//! \brief Get the output value at input 1
		const result_type & GetRight() const
		{
			return m_b;
		}

		standard_linear_scale<R>& operator=(const standard_linear_scale<R>& other)
		{
			m_a = other.m_a;
			m_b = other.m_b;
			return *this;
		}

	private:
		R m_a;
		R m_b;
	};

	/**
	* \brief Derivative of the standard linear scaling
	*
	* The derivative of a linear function is constant, so
	* this class is merely a convenience so that the 
	* automatic differentiation function CTR::get_derivative can
	* work appropriately for standard_linear_scale
	*/
	template <typename R>
	class standard_linear_scale_der
	{
	public:
		typedef R result_type;

      standard_linear_scale_der()
      {
      }

		/*!
			\brief Construct the derivative object from the left and right endpoint values
		*/
		standard_linear_scale_der(const R& left, const R& right) :
			m_slope(right - left)
		{
#ifdef COUNT_OPS
			g_nadd += numel<R>::N;
#endif
		}

		/*!
			\brief Construct the derivative object from the standard_linear_scale<R> object
		*/
		standard_linear_scale_der(const standard_linear_scale<R>& linear_scale) :
			m_slope(linear_scale.GetRight() - linear_scale.GetLeft())
		{
#ifdef COUNT_OPS
			g_nadd += numel<R>::N;
#endif
		}

		/*!
			\brief Get the derivative value

			Returns the constant slope of the linear scaling
		*/
		const result_type& operator()(double s = 0) const
		{
			return m_slope;
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	private:
		R m_slope;
	};
}