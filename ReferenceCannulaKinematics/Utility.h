/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <cstddef>

#include "Linspace.h"
#include "CTRTypedefs.h"
#include <Eigen/Dense>

namespace Utility
{
	inline Eigen::Vector4d	quaternion_mult(Eigen::Vector4d const& q1, Eigen::Vector4d const& q2)
	{
		Eigen::Vector4d 	q;
		q(0) = q1(0)*q2(0) - q1.bottomRows<3>().dot(q2.bottomRows<3>());
#ifdef COUNT_OPS
      g_nmult += 4;
      g_nadd += 1;
#endif
		q.bottomRows<3>() = q1.bottomRows<3>().cross(q2.bottomRows<3>()) + q1(0)*q2.bottomRows<3>() + q2(0)*q1.bottomRows<3>();
#ifdef COUNT_OPS
      g_nmult += 3 + 3 + 6;
      g_nadd += 3 + 3 + 3;
#endif
		return q;
	}

   inline Eigen::Vector4d	quaternion_mult_vec3( Eigen::Vector4d const& q1, Eigen::Vector3d const& w )
   {
      Eigen::Vector4d 	q;
      q( 0 ) = -q1.bottomRows<3>().dot( w ); //5 flops
      q.bottomRows<3>() = q1.bottomRows<3>().cross( w ) + q1( 0 )*w; //9 flops + 3 flops
      return q;
   }

	inline Eigen::Vector4d	quaternion_conj(Eigen::Vector4d const& q)
	{
		Eigen::Vector4d	qc;
		qc(0) = q(0);
		qc.bottomRows<3>() = -q.bottomRows<3>();
		return qc;
	}

   inline Eigen::Vector4d	quaternion_inv( Eigen::Vector4d const& q )
   {
      Eigen::Vector4d	qc;
      qc( 0 ) = q( 0 );
      qc.bottomRows<3>() = -q.bottomRows<3>();
#ifdef COUNT_OPS
      g_ndiv += 1;
      g_nnorm += 1;
      g_nmult += 4;
#endif
      return (1.0/q.norm())*qc;
   }

	inline Eigen::Vector3d	quaternion_rotate(Eigen::Vector4d const& q, const Eigen::Vector3d &p)
	{
		//Eigen::Vector4d	q_p = Eigen::Vector4d::Zero();
		//q_p.bottomRows<3>() = p;
		//return quaternion_mult(quaternion_mult(q, q_p), quaternion_inv(q)).bottomRows<3>();

#ifdef COUNT_OPS
      g_nmult += 3 + 6 + 6 + 3;
      g_nadd += 3 + 3 + 3 + 3;
#endif
      return p + 2.0 * q.bottomRows<3>().cross( q.bottomRows<3>().cross( p ) + q( 0 )*p );
	}

	inline double quaternion_angle(Eigen::Vector4d const& q)
	{
		double chalftheta = q(0);
		double shalftheta = q.bottomRows<3>().norm();
		double halftheta = atan2(shalftheta, chalftheta);
#ifdef COUNT_OPS
      g_nnorm += 1;
      g_ntrig += 1;
      g_nmult += 1;
#endif
		return 2. * halftheta;
	}

	inline Eigen::Vector3d	quaternion_zaxis(Eigen::Vector4d const& q)
	{
		Eigen::Vector3d	Z;
      //double n = q.stableNorm();
      //double nn = 1. / n;
	  double qr = q( 0 );
      double qi = q( 1 );
      double qj = q( 2 );
      double qk = q( 3 );

      //g_nnorm += 1;
      //g_ndiv += 1;
      //g_nmult += 4;


		Z(0) = 2.0*(qi*qk + qj*qr);
		Z(1) = 2.0*(qj*qk - qi*qr);
		Z(2) = 2.0*(qk*qk + qr*qr) - 1.0;

#ifdef COUNT_OPS
      g_nmult += 9;
      g_nadd += 4;
#endif

		return Z;
	}

   inline Eigen::Vector3d	quaternion_yaxis( Eigen::Vector4d const& q )
   {
      Eigen::Vector3d	Y;
      //double n = q.stableNorm();
      //double nn = 1. / n;
      double qr = q( 0 );
      double qi = q( 1 );
      double qj = q( 2 );
      double qk = q( 3 );

      //g_nnorm += 1;
      //g_ndiv += 1;
      //g_nmult += 4;

      Y( 0 ) = 2.0*(qi*qj - qk*qr);
      Y( 1 ) = 2.0*(qj*qj + qr*qr) - 1.0;
      Y( 2 ) = 2.0*(qj*qk + qi*qr);

#ifdef COUNT_OPS
      g_nmult += 9;
      g_nadd += 4;
#endif

      return Y;
   }

   inline Eigen::Vector3d	quaternion_xaxis( Eigen::Vector4d const& q )
   {
      Eigen::Vector3d	X;
      //double n = q.stableNorm();
      //double nn = 1. / n;
      double qr = q( 0 );
      double qi = q( 1 );
      double qj = q( 2 );
      double qk = q( 3 );

      X( 0 ) = 1.0 - 2.0*(qj*qj + qk*qk);
      X( 1 ) = 2.0*(qi*qj + qk*qr);
      X( 2 ) = 2.0*(qi*qk - qj*qr);

      //g_nnorm += 1;
#ifdef COUNT_OPS
      g_nmult += 9;
      g_nadd += 4;
#endif

      return X;
   }

	inline Eigen::Matrix3d	quaternion_to_rotation_matrix(const Eigen::Vector4d &q_)
	{
		Eigen::Matrix3d R;
		//Eigen::Vector4d q = q_ / q_.stableNorm();
#ifdef COUNT_OPS
		//g_nnorm += 1;
		//g_ndiv += 4;
#endif
		R.col(0) = quaternion_xaxis(q_);
		R.col(1) = quaternion_yaxis(q_);
		R.col(2) = quaternion_zaxis(q_);
		return R;
	}

	inline Eigen::Matrix3d	cpm3_2(const Eigen::Vector2d &u)
	{
		Eigen::Matrix3d	uhat;
		uhat(0, 0) = 0;     uhat(0, 1) = 0;     uhat(0, 2) = u(1);
		uhat(1, 0) = 0;     uhat(1, 1) = 0;     uhat(1, 2) = -u(0);
		uhat(2, 0) = -u(1); uhat(2, 1) = u(0);  uhat(2, 2) = 0;
		return uhat;
	}

	inline Eigen::Matrix3d	cpm3(const Eigen::Vector3d &u)
	{
		Eigen::Matrix3d	uhat;
		uhat(0, 0) = 0;     uhat(0, 1) = -u(2); uhat(0, 2) = u(1);
		uhat(1, 0) = u(2);  uhat(1, 1) = 0;     uhat(1, 2) = -u(0);
		uhat(2, 0) = -u(1); uhat(2, 1) = u(0);  uhat(2, 2) = 0;
		return uhat;
	}

	inline Eigen::Matrix<double, 6, 6>	Adjoint_p_q(const Eigen::Vector3d &p, const Eigen::Vector4d &q)
	{
		Eigen::Matrix3d R = quaternion_to_rotation_matrix(q);
		Eigen::Matrix<double, 6, 6>	Ad = Eigen::Matrix<double,6,6>::Zero();
		Ad.topLeftCorner<3, 3>() = R;
		Ad.bottomRightCorner<3, 3>() = R;
		Ad.topRightCorner<3, 3>() = cpm3(p)*R;
#ifdef COUNT_OPS
      g_nmult += 27;
      g_nadd += 18;
#endif
		return Ad;
	}

}
