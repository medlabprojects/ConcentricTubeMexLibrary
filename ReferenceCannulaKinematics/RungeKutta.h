/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <algorithm>
#include <numeric>
#include "EvalPt.h"

namespace Mathematics
{

	template <class State, class Functor, class InIt, class OutIt>
	OutIt integrate(InIt first, InIt last, OutIt o_first, const State& y_init, const Functor& f);

	template <class State, class Functor, class TimeVector>
	State euler_final(Functor const& f, State const& y_init, TimeVector const& time);

	template <class State, class Functor, class TimeVector>
	State rk4_final(Functor const& f, State const& y_init, TimeVector const& time);

	namespace RK8_Coeffs
	{
		extern double c[13];
		extern double a[13][12];
		extern double b[12];
		extern double bhat[13];
	}

	template <class State, class Functor, class TimeVector>
	State rk8_final(Functor const& f, State const& y_init, TimeVector const& time);

	template <class State>
	struct StateTime
	{
		StateTime(const State& s, double t) : s_(s), t_(t) {}
		State s_;
		double t_;
	};

   namespace internal
   {
      struct empty_callback
      {
         template <typename State>
         void operator()( State &y ) const
         {
         }
      };

   }

	namespace internal
	{
		template <class Functor>
		struct euler_step_obj
		{
			template <typename State>
			struct result
			{
				typedef StateTime<State>	type;
			};

			euler_step_obj(Functor const& f) : f_(f) {}

			template <class State>
			StateTime<State> operator()(StateTime<State> const& prev, double next_time) const
			{
				double tn = prev.t_;
				double h = next_time - tn;
				State k1 = f_(eval_pt(tn, tn, next_time), prev.s_);
				return StateTime<State>(prev.s_ + h*k1, next_time);
			}

			Functor const& f_;
		};

		template <class Functor>
		struct rk4_step_obj
		{
			template <typename State>
			struct result
			{
				typedef StateTime<State>	type;
			};

			rk4_step_obj(Functor const& f) : f_(f) {}

			template <class State>
			StateTime<State> operator()(StateTime<State> const& prev, double next_time) const
			{
				double tn = prev.t_;
				double h = next_time - tn;
				State k1 = f_(eval_pt(tn,tn,next_time), prev.s_);
				State k2 = f_(eval_pt(tn + (h / 3.0),tn,next_time), prev.s_ + (h / 3.0)*k1);
				State k3 = f_(eval_pt(tn + (2.0 * h / 3.0),tn,next_time), prev.s_  + (-h / 3.0)*k1 + h*k2);
				State k4 = f_(eval_pt(tn + h,tn,next_time), prev.s_ + h*k1 - h*k2 + h*k3);
				return StateTime<State>(prev.s_ + h / 8.0*(k1 + 3.0*k2 + 3.0*k3 + k4), next_time);
			}

			Functor const& f_;
		};

		template <class Functor>
		struct rk4_richardson_step_obj
		{
			template <typename State>
			struct result
			{
				typedef StateTime<State>	type;
			};

			rk4_richardson_step_obj(Functor const& f) : f_(f) {}

			template <class State>
			StateTime<State> operator()(StateTime<State> const& prev, double next_time) const
			{
				rk4_step_obj<Functor> rk4_step(f_);

				StateTime<State> yhat1 = rk4_step(prev, next_time);
				StateTime<State> half = rk4_step(prev, (next_time - prev.t_) / 2.0 + prev.t_);
				StateTime<State> yhat2 = rk4_step(half, next_time);

				//combine yhat1 and yhat2 to get a more accurate solution estimate by
				// richardson extrapolation
				State	extrap = (1.0 / (31.0))*(32.0*yhat2.s_ - yhat1.s_);
				return StateTime<State>(extrap, next_time);
			}

			Functor const& f_;
		};


		template <class Functor>
		struct rk8_step_obj
		{
			template <typename State>
			struct result
			{
				typedef StateTime<State>	type;
			};

			rk8_step_obj(Functor const& f) : f_(f) {}

			template <class State>
			StateTime<State> operator()(StateTime<State> const& prev, double next_time) const
			{
				double tn = prev.t_;
				double h = next_time - tn;

				double eval_times[13];
				std::transform(std::begin(RK8_Coeffs::c), std::end(RK8_Coeffs::c), std::begin(eval_times),
				[tn,h](double c) {
					return tn + h*c;
				});

				State k[13];

				k[0] = f_(eval_pt(eval_times[0],tn,next_time), prev.s_);
				k[1] = f_(eval_pt(eval_times[1],tn,next_time),
					prev.s_
					+ h* RK8_Coeffs::a[1][0] * k[0]);

				k[2] = f_(eval_pt(eval_times[2],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[2][0] * k[0]
					+ h*RK8_Coeffs::a[2][1] * k[1]);

				k[3] = f_(eval_pt(eval_times[3],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[3][0] * k[0]
					+ h*RK8_Coeffs::a[3][1] * k[1]
					+ h*RK8_Coeffs::a[3][2] * k[2]);

				k[4] = f_(eval_pt(eval_times[4],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[4][0] * k[0]
					+ h*RK8_Coeffs::a[4][1] * k[1]
					+ h*RK8_Coeffs::a[4][2] * k[2]
					+ h*RK8_Coeffs::a[4][3] * k[3]);

				k[5] = f_(eval_pt(eval_times[5],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[5][0] * k[0]
					+ h*RK8_Coeffs::a[5][1] * k[1]
					+ h*RK8_Coeffs::a[5][2] * k[2]
					+ h*RK8_Coeffs::a[5][3] * k[3]
					+ h*RK8_Coeffs::a[5][4] * k[4]);
				
				k[6] = f_(eval_pt(eval_times[6],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[6][0] * k[0]
					+ h*RK8_Coeffs::a[6][1] * k[1]
					+ h*RK8_Coeffs::a[6][2] * k[2]
					+ h*RK8_Coeffs::a[6][3] * k[3]
					+ h*RK8_Coeffs::a[6][4] * k[4]
					+ h*RK8_Coeffs::a[6][5] * k[5]);

				k[7] = f_(eval_pt(eval_times[7],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[7][0] * k[0]
					+ h*RK8_Coeffs::a[7][1] * k[1]
					+ h*RK8_Coeffs::a[7][2] * k[2]
					+ h*RK8_Coeffs::a[7][3] * k[3]
					+ h*RK8_Coeffs::a[7][4] * k[4]
					+ h*RK8_Coeffs::a[7][5] * k[5]
					+ h*RK8_Coeffs::a[7][6] * k[6]);

				k[8] = f_(eval_pt(eval_times[8],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[8][0] * k[0]
					+ h*RK8_Coeffs::a[8][1] * k[1]
					+ h*RK8_Coeffs::a[8][2] * k[2]
					+ h*RK8_Coeffs::a[8][3] * k[3]
					+ h*RK8_Coeffs::a[8][4] * k[4]
					+ h*RK8_Coeffs::a[8][5] * k[5]
					+ h*RK8_Coeffs::a[8][6] * k[6]
					+ h*RK8_Coeffs::a[8][7] * k[7]);

				k[9] = f_(eval_pt(eval_times[9],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[9][0] * k[0]
					+ h*RK8_Coeffs::a[9][1] * k[1]
					+ h*RK8_Coeffs::a[9][2] * k[2]
					+ h*RK8_Coeffs::a[9][3] * k[3]
					+ h*RK8_Coeffs::a[9][4] * k[4]
					+ h*RK8_Coeffs::a[9][5] * k[5]
					+ h*RK8_Coeffs::a[9][6] * k[6]
					+ h*RK8_Coeffs::a[9][7] * k[7]
					+ h*RK8_Coeffs::a[9][8] * k[8]);

				k[10] = f_(eval_pt(eval_times[10],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[10][0] * k[0]
					+ h*RK8_Coeffs::a[10][1] * k[1]
					+ h*RK8_Coeffs::a[10][2] * k[2]
					+ h*RK8_Coeffs::a[10][3] * k[3]
					+ h*RK8_Coeffs::a[10][4] * k[4]
					+ h*RK8_Coeffs::a[10][5] * k[5]
					+ h*RK8_Coeffs::a[10][6] * k[6]
					+ h*RK8_Coeffs::a[10][7] * k[7]
					+ h*RK8_Coeffs::a[10][8] * k[8]
					+ h*RK8_Coeffs::a[10][9] * k[9]);

				k[11] = f_(eval_pt(eval_times[11],tn,next_time),
					prev.s_
					+ h*RK8_Coeffs::a[11][0] * k[0]
					+ h*RK8_Coeffs::a[11][1] * k[1]
					+ h*RK8_Coeffs::a[11][2] * k[2]
					+ h*RK8_Coeffs::a[11][3] * k[3]
					+ h*RK8_Coeffs::a[11][4] * k[4]
					+ h*RK8_Coeffs::a[11][5] * k[5]
					+ h*RK8_Coeffs::a[11][6] * k[6]
					+ h*RK8_Coeffs::a[11][7] * k[7]
					+ h*RK8_Coeffs::a[11][8] * k[8]
					+ h*RK8_Coeffs::a[11][9] * k[9]
					+ h*RK8_Coeffs::a[11][10] * k[10]);

				k[12] = f_(eval_pt(eval_times[12], tn, next_time),
					prev.s_
					+ h*RK8_Coeffs::a[12][0] * k[0]
					+ h*RK8_Coeffs::a[12][1] * k[1]
					+ h*RK8_Coeffs::a[12][2] * k[2]
					+ h*RK8_Coeffs::a[12][3] * k[3]
					+ h*RK8_Coeffs::a[12][4] * k[4]
					+ h*RK8_Coeffs::a[12][5] * k[5]
					+ h*RK8_Coeffs::a[12][6] * k[6]
					+ h*RK8_Coeffs::a[12][7] * k[7]
					+ h*RK8_Coeffs::a[12][8] * k[8]
					+ h*RK8_Coeffs::a[12][9] * k[9]
					+ h*RK8_Coeffs::a[12][10] * k[10]
					+ h*RK8_Coeffs::a[12][11] * k[11]);

				return StateTime<State>(prev.s_
					+ h*RK8_Coeffs::bhat[0] * k[0]
					+ h*RK8_Coeffs::bhat[1] * k[1]
					+ h*RK8_Coeffs::bhat[2] * k[2]
					+ h*RK8_Coeffs::bhat[3] * k[3]
					+ h*RK8_Coeffs::bhat[4] * k[4]
					+ h*RK8_Coeffs::bhat[5] * k[5]
					+ h*RK8_Coeffs::bhat[6] * k[6]
					+ h*RK8_Coeffs::bhat[7] * k[7]
					+ h*RK8_Coeffs::bhat[8] * k[8]
					+ h*RK8_Coeffs::bhat[9] * k[9]
					+ h*RK8_Coeffs::bhat[10] * k[10]
					+ h*RK8_Coeffs::bhat[11] * k[11]
					+ h*RK8_Coeffs::bhat[12] * k[12], next_time);
			}

			Functor const& f_;
		};
	}

	template <typename Functor>
	internal::rk4_step_obj<Functor> rk4_step(Functor const& F)
	{
		return internal::rk4_step_obj<Functor>(F);
	}

	template <typename Functor>
	internal::rk4_richardson_step_obj<Functor> rk4_richardson_step(Functor const& F)
	{
		return internal::rk4_richardson_step_obj<Functor>(F);
	}

	template <typename Functor>
	internal::rk8_step_obj<Functor> rk8_step(Functor const& F)
	{
		return internal::rk8_step_obj<Functor>(F);
	}

	template <typename Functor>
	internal::euler_step_obj<Functor>	euler_step(Functor const& F)
	{
		return internal::euler_step_obj<Functor>(F);
	}

	/**
	\brief Euler integrator template

	\param f	The callable object representing the first-order state equations
	\param y_init	The initial state vector
	\param time		The vector of integration times
	\return	The state at the end of the integration

	This template for the classical 1st order Euler method
	takes inputs of any types and performs the numerical integration.
	This template returns only the final state value.

	The requirements of the input types are

	Functor const& f must be a callable type, that is f must be able to
	be invoked

	S = f(t,y)

	for a time t and state y.

	The State type passed in as y_init must have the basic arithmetic operations

	operator+(const State& s1, const State& s2)
	operator*(double, const State& s)

	The time vector must be compatible with the standard template library
	algorithms and must have the const_iterator type defined for it, i.e.
	TimeVector::const_iterator must be an iterator which cannot change the
	contents of the TimeVector. Typically, the time vector should just be an
	std::vector type, possibly with a custom allocator if hard real-time
	support is needed.

	A single step forward of the Euler method is given by the following algorithm steps:
	~~~{.cpp}
	double h = next_time - current_time; //step length
	State k1 = f_(tn, prev);
	State next = prev + h*k1;
	~~~

	*/
	template <class State, class Functor, class TimeVector>
	State euler_final(Functor const& f, State const& y_init, TimeVector const& time)
	{
		typename TimeVector::const_iterator t = std::begin(time);
		typename TimeVector::const_iterator end_time = std::end(time);

		StateTime<State>	st(y_init, *t);
		++t;

		for (; t != end_time; ++t)
		{
			st = euler_step<Functor>(f)(st, *t);
		}

		return st.s_;
	}

	template <class State, class Functor, class Iter>
	State euler_final(Functor const& f, State const& y_init, Iter start_time, Iter const& end_time)
	{
		StateTime<State>	st(y_init, *start_time);
		++start_time;

		for (; start_time != end_time; ++start_time)
		{
			st = euler_step<Functor>(f)(st, *start_time);
		}

		return st.s_;
	}

	/**
		\brief Runge Kutta template

		\param f	The callable object representing the first-order state equations
		\param y_init	The initial state vector
		\param time		The vector of integration times
		\return	The state at the end of the integration

		This template for the classical 4th order Runge-Kutta method
		takes inputs of any types and performs the numerical integration.
		This template returns only the final state value.

		The requirements of the input types are

		Functor const& f must be a callable type, that is f must be able to
		be invoked

		S = f(t,y)

		for a time t and state y.

		The State type passed in as y_init must have the basic arithmetic operations

		operator+(const State& s1, const State& s2)
		operator*(double, const State& s)

		The time vector must be compatible with the standard template library
		algorithms and must have the const_iterator type defined for it, i.e.
		TimeVector::const_iterator must be an iterator which cannot change the
		contents of the TimeVector. Typically, the time vector should just be an
		std::vector type, possibly with a custom allocator if hard real-time
		support is needed.

		A single step forward of the Runge Kutta method is given by the following algorithm steps:

		~~~{.cpp}
		double h = next_time - current_time; //step length
		State k1 = f_(tn, prev);
		State k2 = f_(tn + h / 2.0, prev + h / 2.0*k1);
		State k3 = f_(tn + h / 2.0, prev + h / 2.0*k2);
		State k4 = f_(tn + h, prev + h*k3);
		State next = prev + h / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		~~~

		*/
	template <class State, class Functor, class TimeVector>
	State rk4_final(Functor const& f, State const& y_init, TimeVector const& time)
	{
		typename TimeVector::const_iterator t = std::begin(time);
		typename TimeVector::const_iterator end_time = std::end(time);

		StateTime<State>	st(y_init, *t);
		++t;

		for (; t != end_time; ++t)
		{
			st = rk4_step<Functor>(f)(st, *t);
		}

		return st.s_;
	}

	template <class State, class Functor, class Iter>
	State rk4_final(Functor const& f, State const& y_init, Iter start_time, Iter const& end_time)
	{
		StateTime<State>	st(y_init, *start_time);
		++start_time;

		for (; start_time != end_time; ++start_time)
		{
			st = rk4_step<Functor>(f)(st, *start_time);
		}

		return st.s_;
	}

   template <class State, class Functor, class TimeVector, class StateCallback>
   State rk4_final( Functor const& f, State const& y_init, TimeVector const& time, StateCallback const &callback )
   {
      typename TimeVector::const_iterator t = std::begin( time );
      typename TimeVector::const_iterator end_time = std::end( time );

      StateTime<State>	st( y_init, *t );
      ++t;

      for ( ; t != end_time; ++t )
      {
         st = rk4_step<Functor>( f )(st, *t);
         callback( st.s_ );
      }

      return st.s_;
   }

   template <class State, class Functor, class Iter, class StateCallback>
   State rk4_final( Functor const& f, State const& y_init, Iter start_time, Iter const& end_time, StateCallback const &callback )
   {
      StateTime<State>	st( y_init, *start_time );
      ++start_time;

      for ( ; start_time != end_time; ++start_time )
      {
         st = rk4_step<Functor>( f )(st, *start_time);
         callback( st.s_ );
      }

      return st.s_;
   }

	template <class State, class Functor, class TimeVector>
	State rk4_richardson_final(Functor const& f, State const& y_init, TimeVector const& time)
	{
		typename TimeVector::const_iterator t = std::begin(time);
		typename TimeVector::const_iterator end_time = std::end(time);

		StateTime<State>	st(y_init, *t);
		++t;

		for (; t != end_time; ++t)
		{
			st = rk4_richardson_step<Functor>(f)(st, *t);
		}

		return st.s_;
	}

	template <class State, class Functor, class Iter>
	State rk4_richardson_final(Functor const& f, State const& y_init, Iter start_time, Iter const& end_time)
	{
		StateTime<State>	st(y_init, *start_time);
		++start_time;

		for (; start_time != end_time; ++start_time)
		{
			st = rk4_richardson_step<Functor>(f)(st, *start_time);
		}

		return st.s_;
	}

	/**
		\brief Runge Kutta template

		\param f	The callable object representing the first-order state equations
		\param y_init	The initial state vector
		\param time		The vector of integration times
		\return	The state at the end of the integration

		This template for the classical 4th order Runge-Kutta method
		takes inputs of any types and performs the numerical integration.
		This template returns only the final state value.

		The requirements of the input types are

		Functor const& f must be a callable type, that is f must be able to
		be invoked

		S = f(t,y)

		for a time t and state y.

		The State type passed in as y_init must have the basic arithmetic operations

		operator+(const State& s1, const State& s2)
		operator*(double, const State& s)

		The time vector must be compatible with the standard template library
		algorithms and must have the const_iterator type defined for it, i.e.
		TimeVector::const_iterator must be an iterator which cannot change the
		contents of the TimeVector. Typically, the time vector should just be an
		std::vector type, possibly with a custom allocator if hard real-time
		support is needed.

		A single step forward of the Runge Kutta method is given by the following algorithm steps:

		~~~{.cpp}
		double h = next_time - current_time; //step length
		State k1 = f_(tn, prev);
		State k2 = f_(tn + h / 2.0, prev + h / 2.0*k1);
		State k3 = f_(tn + h / 2.0, prev + h / 2.0*k2);
		State k4 = f_(tn + h, prev + h*k3);
		State next = prev + h / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		~~~

	*/
	template <class State, class Functor, class TimeVector>
   State rk8_final( Functor const& f, State const& y_init, TimeVector const& time )
   {
      typename TimeVector::const_iterator t = std::begin( time );
      typename TimeVector::const_iterator end_time = std::end( time );

      StateTime<State>	st( y_init, *t );
      ++t;

      for ( ; t != end_time; ++t )
      {
         st = rk8_step<Functor>( f )(st, *t);
      }

      return st.s_;
   }
	

   template <class State, class Functor, class TimeVector, class StepCallback>
   State rk8_final( Functor const& f, State const& y_init, TimeVector const& time, StepCallback const& callback )
   {
      typename TimeVector::const_iterator t = std::begin( time );
      typename TimeVector::const_iterator end_time = std::end( time );

      StateTime<State>	st( y_init, *t );
      ++t;

      for ( ; t != end_time; ++t )
      {
         st = rk8_step<Functor>( f )(st, *t);
         callback( st.s_ );
      }

      return st.s_;
   }

	template <class State, class Functor, class Iter>
	State rk8_final(Functor const& f, State const& y_init, Iter start_time, Iter const& end_time)
	{
		StateTime<State>	st(y_init, *start_time);
		++start_time;

		for (; start_time != end_time; ++start_time)
		{
			st = rk8_step<Functor>(f)(st, *start_time);
		}

		return st.s_;
	}

   template <class State, class Functor, class Iter, class Callback>
   State rk8_final( Functor const& f, State const& y_init, Iter start_time, Iter const& end_time, Callback const &callback )
   {
      StateTime<State>	st( y_init, *start_time );
      ++start_time;

      for ( ; start_time != end_time; ++start_time )
      {
         st = rk8_step<Functor>( f )(st, *start_time);
         callback( st.s_ );
      }

      return st.s_;
   }

	/**
		\brief Runge Kutta template

		\param f	The callable object representing the first-order state equations
		\param y_init	The initial state vector
		\param time		The vector of integration times
		\return	The state at the end of the integration

		This template for the classical 4th order Runge-Kutta method
		takes inputs of any types and performs the numerical integration.
		This template returns only the final state value.

		The requirements of the input types are

		Functor const& f must be a callable type, that is f must be able to
		be invoked

		S = f(t,y)

		for a time t and state y.

		The State type passed in as y_init must have the basic arithmetic operations

		operator+(const State& s1, const State& s2)
		operator*(double, const State& s)

		The time vector must be compatible with the standard template library
		algorithms and must have the const_iterator type defined for it, i.e.
		TimeVector::const_iterator must be an iterator which cannot change the
		contents of the TimeVector. Typically, the time vector should just be an
		std::vector type, possibly with a custom allocator if hard real-time
		support is needed.

		A single step forward of the Runge Kutta method is given by the following algorithm steps:

		~~~{.cpp}
		double h = next_time - current_time; //step length
		State k1 = f_(tn, prev);
		State k2 = f_(tn + h / 2.0, prev + h / 2.0*k1);
		State k3 = f_(tn + h / 2.0, prev + h / 2.0*k2);
		State k4 = f_(tn + h, prev + h*k3);
		State next = prev + h / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		~~~

		*/
	template <class State, class Functor, class InIt, class OutIt, class Callback>
	OutIt integrate(InIt first, InIt last, OutIt o_first, const State& y_init, const Functor& f, Callback const &cback = Callback())
	{
		StateTime<State> st(y_init, *first);
		*o_first = y_init;
		while (++first != last)
		{
			st = f(st, *first);
			cback(st.s_);
			*++o_first = st.s_;
		}
		return o_first;
	}

	template <class State, class Functor, class InIt, class OutIt>
	OutIt integrate(InIt first, InIt last, OutIt o_first, const State& y_init, const Functor& f)
	{
		StateTime<State> st(y_init, *first);
		*o_first = y_init;
		while (++first != last)
		{
			st = f(st, *first);
			*++o_first = st.s_;
		}
		return o_first;
	}

}
