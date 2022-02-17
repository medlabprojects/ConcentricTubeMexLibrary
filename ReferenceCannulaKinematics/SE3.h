#pragma once

#include "CTRTypedefs.h"
#include <Eigen/Dense>
#include <tuple>
#include "Utility.h"

namespace SE3
{	
	struct transformation
	{
		transformation()
		{}

		transformation(Eigen::Vector3d const &p_, Eigen::Vector4d const &q_)
			:p(p_)
			,q(q_)
		{}
		
		Eigen::Vector3d	p;
		Eigen::Vector4d	q;

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	struct transformation_tie
	{
		transformation_tie(Eigen::Vector3d &p_, Eigen::Vector4d &q_)
			: p(p_)
			, q(q_)
		{}

		transformation_tie& operator=(transformation const &t)
		{
			p = t.p;
			q = t.q;
            return *this;
		}

		Eigen::Vector3d	&p;
		Eigen::Vector4d	&q;
	};

	transformation make_transformation(Eigen::Vector3d const &p, Eigen::Vector4d const &q)
	{
		return transformation(p,q);	
	}

	transformation_tie	tie_transformation(Eigen::Vector3d &p, Eigen::Vector4d &q)
	{
		return transformation_tie(p,q);
	}

	inline
	transformation transformation_inverse( Eigen::Vector3d const &p, Eigen::Vector4d const &q )
	{
		using Utility::quaternion_rotate;
		using Utility::quaternion_inv;
		transformation R;
		R.q = quaternion_inv(q);
		R.p = -quaternion_rotate(R.q, p);
		return R;
	}

	inline
	transformation operator*(transformation const &t1, transformation const &t2)
	{
		using Utility::quaternion_rotate;
		using Utility::quaternion_mult;
		return make_transformation(t1.p + quaternion_rotate(t1.q, t2.p), quaternion_mult(t1.q, t2.q) );
	}
}
