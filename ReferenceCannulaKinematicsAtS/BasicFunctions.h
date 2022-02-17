/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <Eigen/Dense>
#include <boost/utility.hpp>
#include <boost/utility/result_of.hpp>
#include "EvalPt.h"

namespace CTR
{
	/////////////////////////////////////////
	//                                     //
	//  Forward Declarations               //
	//                                     //
	/////////////////////////////////////////

	namespace Functions
	{
		template <typename R>
		struct constant_fun;
	}

	namespace internal
	{
		struct constant_fun_der;
	}

	namespace Functions
	{
		template <typename R>
		struct zero_fun;

		struct derivative_helper;
	}

	/////////////////////////////////////////
	//                                     //
	//  Constant Function                  //
	//                                     //
	/////////////////////////////////////////
	namespace internal {
		template <typename R>
		struct traits;

		template <typename R>
		struct traits < Functions::constant_fun<R> >
		{
			typedef constant_fun_der	derivative_maker;
		};
	}

	namespace Functions
	{
		/**
			\brief Models the unary constant function f(x) = c

			Retains a copy of the constant c, and evaluates to this constant
			for all input arguments.

			~~~{.cpp}
			constant_fun<double> f(4);
			cout << f(0.0) << endl << f(1.0) << endl << f(2.0) << endl;
			~~~

			~~~
			Output:
			4
			4
			4
			~~~
		*/
		template <typename R>
		struct constant_fun
		{
			typedef R result_type;
			constant_fun() {}
			constant_fun(const R& r) : m_r(r) {}

			result_type	operator()(double x = 0.0) const {
				return m_r;
			}

			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		private:
			R m_r;
		};
	}

	namespace internal
	{
		struct constant_fun_der
		{
			template <typename T>
			struct result;

			template <typename T>
			struct result < constant_fun_der(Functions::constant_fun<T> const&) >
			{
				typedef Functions::zero_fun<T>	type;
			};

			template <typename T>
			Functions::zero_fun<T>	operator()(Functions::constant_fun<T> const&)
			{
				return Functions::zero_fun<T>();
			}
		};
	}


	/////////////////////////////////////////
	//                                     //
	//  Zero Functions                     //
	//                                     //
	/////////////////////////////////////////
	namespace Functions
	{

		/*!
			\brief Implements the unary ``zero'' function 0(x) for an arbitrary type

			This class is a function object which returns the value zero
			for an arbitrary type, including integer, floating point,
			and Matrix types.

			Note that the function has a default argument provided, so that if you
			simply need the representation of zero for a particular type,

			Example Code:

				 double a = zero_fun<double>()(); //a = 0
				 Eigen::Vector3d b = zero_fun<Eigen::Vector3d>()(); //b = [0,0,0]
		*/
		template <typename T>
		struct zero_fun
		{
			typedef T result_type;

			zero_fun() {}

			/*!
				\brief Return the zero value
			*/
			result_type operator()(double s = 0) const
			{
				return static_cast<T>(0);
			}
		};

		template < typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
		struct zero_fun < Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> >
		{
			typedef Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> result_type;
			zero_fun() {}

			result_type operator()(double s = 0) const
			{
				return result_type::Zero();
			}
		};
	}

	namespace Functions
	{
		template <typename T>
		struct step_function;
	}

	namespace internal
	{
		struct step_fun_der;
	}

	namespace internal
	{
		template <typename T> struct traits;

		//
		// The derivative of the step function
		// is zero almost everywhere
		//
		template <typename T>
		struct traits < Functions::step_function<T> >
		{
			typedef step_fun_der	derivative_maker;
		};
	}

	namespace Functions
	{
		/**
			\brief Models the step function u(x)

			This function models the step function u(x) with user-selectable
			values for the left and right of the step and the location of the
			step.

			~~~{.cpp}
			step_function<double>	u(0.0, 1.0, 0.0);
			cout << u(-1.0e-6) << endl
				 << u(1.0e-6) << endl;
			~~~

			~~~
			Output:
			0
			1
			~~~

		*/
		template <typename T>
		struct step_function
		{
			step_function(T left_, T right_, double sw_) :
				sw(sw_),
				left(left_),
				right(right_)
			{}

			T operator()(double s) const
			{
				return (s > sw) ? right : left;
			}

			EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		private:
			double sw;
			T left;
			T right;
		};

	}

	namespace internal
	{
		struct step_fun_der
		{
			template <typename T>
			struct result;

			template <typename T>
			struct result < step_fun_der(T const&) >
			{
				typedef Functions::zero_fun<T>	type;
			};

			template <typename T>
			Functions::zero_fun<T> operator()(Functions::step_function<T> const&)
			{
				return Functions::zero_fun<T>();
			}
		};
	}

	namespace Functions
	{
		struct indicator_function;
	}

	namespace internal
	{
		struct indicator_fun_der;
	}

	namespace internal
	{
		template <typename T> struct traits;

		template <>
		struct traits < Functions::indicator_function >
		{
			typedef indicator_fun_der	derivative_maker;
		};
	}

	namespace Functions
	{
		/**
			\brief Models the indicator function for an interval

			This function models the indicator for an interval of the real line.
			The left and right values are the endpoints of the interval.
		*/
		struct indicator_function
		{
			indicator_function() {}

			indicator_function(double left_, double right_) :
				left(left_),
				right(right_)
			{}


			double operator()(double s) const
			{
				return (s >= left && s <= right) ? 1.0 : 0.0;
			}

			indicator_function& operator=(indicator_function const& other)
			{
				left = other.left;
				right = other.right;
				return *this;
			}

			double left;
			double right;
		};

	}

	namespace internal
	{
		struct indicator_fun_der
		{
			typedef Functions::zero_fun<double>	result_type;

			Functions::zero_fun<double>	operator()(Functions::indicator_function const&)
			{
				return Functions::zero_fun<double>();
			}
		};
	}

	/////////////////////////////////////////
	//                                     //
	//  Automatic Derivative Handling      //
	//                                     //
	/////////////////////////////////////////

	namespace internal {

		struct derivative_helper;

		/*!
		\brief This structure is used by CTR::get_derivative to deduce the return type

		For a function object type F to have a derivative computed automatically,
		there are two requirements. First, the structure CTR::internal::traits<F>
		must be defined for the function type F whose derivative is desired, and
		this structure must have a member typedef called derivative_type which
		defines the function object type of the derivative.

		Secondly, the derivative function object type needs to have a constructor
		which takes a constant reference to the original function object whose
		derivative is desired.

		If the function object type user_function has been defined already, and
		the derivative has been coded as the function object type derivative_type,
		then the following example code is sufficient to define the derivative type for
		get_derivative.

		The following code example shows the basics, defining a function type called
		user_function, which is the function user_function: \f$x \to m*x\f$, where m is
		passed into the constructor of the function object. The derivative function
		takes this slope from the function object reference and stores it to be
		returned as the derivative with respect to x.

		~~~{.cpp}
		template <typename R>
		struct user_function;

		template <typename R>
		struct user_function_derivative;

		namespace CTR
		{
		namespace internal
		{
		template <typename R> struct traits;

		template <typename R>
		struct traits< ::user_function<R> >
		{
		typedef ::user_function_derivative<R>	derivative_type;
		};
		}
		}

		template <typename R>
		struct user_function
		{
		user_function(R slope) : m(slope)

		R operator()(R x)
		{
		return m*x;
		}

		const R& GetSlope() const { return m; }
		private:
		R m;
		};

		template <typename R>
		struct user_function_derivative
		{
		///constructor
		user_function_derivative( const user_function<R>& fun ) :
		m(fun.GetSlope())
		{}

		R operator()(R x)
		{
		return m;
		}

		private:
		R m;
		};
		~~~
		*/
		struct derivative_helper
		{
			template <typename R>
			struct result;

			template <class T>
			struct result<derivative_helper(const T&)>
			{
				typedef typename traits<T>::derivative_maker	maker_type;
				typedef typename std::result_of< maker_type(const T&) >::type	type;
			};

			template <typename F>
			typename result< derivative_helper(const F&) >::type operator()(const F &fun) const
			{
				return typename traits<F>::derivative_maker()(fun);
			}
		};
	}

	/*!
	\brief Get the derivative of a function object which supports derivatives

	For most of the continous (i.e. at least C1) function types defined in this library,
	the derivative function object can be obtained through this helper. Whether this
	is supported will be noted in the documentation for the function object.

	Example Code:

	~~~{.cpp}
	//Create a smooth step function with double type, where the function
	// is 0 on the left, 10 on the right, switches at 5, and has a step
	// width of 1
	CTR::smooth_step_function<double> step_fun(0.0,10.0,5.0,1.0);
	double a = step_fun(5.0); //a = 5.0
	double b = get_derivative(step_fun)(5.0); //b = derivative of step_fun at input 5.0
	auto step_fun_der = get_derivative(step_fun); //step_fun_der is a function object representing the derivative
	~~~

	If you wish to write a new function object type which is compatible with get_derivative, see
	the documentation for CTR::internal::derivative_helper
	*/
	template <typename R>
	typename std::result_of< internal::derivative_helper(const R&) >::type
		get_derivative(const R& fun)
	{
		return internal::derivative_helper()(fun);
	}

	/*
		Get the number of elements in a scalar/vector
	*/
	template <typename T>
	struct numel
	{
		const static int N = T::SizeAtCompileTime;
	};

	//specialization for double/float
	template <>
	struct numel<double>
	{
		const static int N = 1;
	};

	template <>
	struct numel<float>
	{
		const static int N = 1;
	};

}