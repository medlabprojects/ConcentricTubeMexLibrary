#pragma once

namespace Mathematics
{

	/**
		\brief An evaluation point for the differential equation

		Since the differential equations involve discontinuous coefficients,
		we choose to evaluate the function at a combination of a single time (arc length)
		point and an interval of integration. In this way, the values of the discontinuous
		coefficients can depend on the interval as well as the actual time (arc length)
		variable. This struct is a simple wrapper, together with an implicit conversion
		to double so that an eval_pt can be passed to a function which takes only
		a single double argument.
	*/
	struct eval_pt
	{
		/**
			\param s_ the time value
			\param left_ the left of the evaluation interval
			\param right_ the right of the evaluation interval
		*/
		eval_pt(double s_, double left_, double right_) :
			s(s_),
			left(left_),
			right(right_)
		{}

		/**
			\brief Enable implicit conversion to double
			\return the time variable value
		*/
		operator double() const
		{
			return s;
		}

		double s;
		double left; double right;
	};

}