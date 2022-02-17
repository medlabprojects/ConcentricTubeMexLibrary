/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include <type_traits>
#include "CTRTypedefs.h"
#include "TubePhysicalProps.h"
#include "BasicFunctions.h"
#include "IndicatorFunction.h"
#include "RungeKutta.h"
#include "EvalPt.h"

namespace CTR
{
	template <typename CurvFunType, typename KbFunType, typename CtFunType>
	class Tube;

	namespace internal
	{
		template <typename T>
		struct traits;

		template <typename CurvFunType, typename KbFunType, typename CtFunType>
		struct traits < Tube<CurvFunType, KbFunType, CtFunType> >
		{
		};
	}

	/**
		\brief The Tube class, which represents a single tube of a concentric tube robot

		This class completely encapsulates the details of a single tube of a concentric tube robot, including
		all necessary properties and the curvature function \f$\vec{\kappa}(s): \mathbb{R} \to \mathbb{R}^2\f$.
	*/
	template <typename CurvFunType, 
			  typename KbFunType = Functions::constant_fun<double>, 
			  typename CtFunType = Functions::constant_fun<double> >
	class Tube
	{
	public:
		Tube() {}

		Tube(const Tube& other) = default;

		/**
			\brief Construct a tube from parameters

			\param [in] L	The length of the tube (total length in meters)
			\param [in] k	The curvature function of the tube (must be a unary function returning a 2d vector)
			\param [in] kb	The bending stiffness of the tube
			\param [in] ct	The torsional compliance (inverse of torsional stiffness)
		*/
		Tube(double L, double tL, CurvFunType const& k, KbFunType const& kb, CtFunType const& ct) :
			m_L(L),
			m_tL(tL),
			m_ustar(k),
			m_kb(kb),
			m_ct(ct),
			m_chi(0.0,L),
			m_uChi(tL,L,1.0e-3),
         m_Dustar(get_derivative(m_ustar)),
         m_DuChi(get_derivative(m_uChi))
		{}

		/**
			\brief Get the precurvature of the tube at the point s

			The return type of this function is determined by the CurvFunType
			template parameter, but the only return type currently supported
			is they CTR::Vector<2>::type, which is an alias of Eigen::Vector2d or similar
		*/
		typename Vector<2>::type
			GetPrecurvature(double s) const
		{
			return m_ustar(s)*m_uChi(s);
		}

		typename Vector<2>::type
			GetDPrecurvature(double s) const
		{
			//Product rule for derivative
#ifdef COUNT_OPS
			g_nadd += 2;
			g_nmult += 4;
#endif
			//return get_derivative(m_ustar)(s)*m_uChi(s) + m_ustar(s)*get_derivative(m_uChi)(s);
			return m_Dustar(s)*m_uChi(s) + m_ustar(s)*m_DuChi(s);
		}

		/**
			\brief Get the total length
		*/
		double GetLength() const
		{
			return m_L;
		}

		/**
			\brief Get the straight transmission length
		*/
		double GetTransmissionLength() const
		{
			return m_tL;
		}

		/**
			\brief Get the bending stiffness
		*/
		double GetBendingStiffness(Mathematics::eval_pt s) const
		{
#ifdef COUNT_OPS
			g_nadd += 1;
			g_nmult += 1;
#endif
			double ival_center = (s.left + s.right) * 0.5;
			return m_kb(s.s)*m_chi(ival_center);
		}

		/**
			\brief Get the torsional compliance
		*/
		double GetTorsionalCompliance(Mathematics::eval_pt s) const
		{
			double ival_center = (s.left + s.right) * 0.5;
#ifdef COUNT_OPS
			g_nadd += 1;
			g_nmult += 1;
#endif
			return m_ct(s.s)*m_chi(ival_center);
		}

		/**
			\brief Get the derivative of bending stiffness
		*/
		double GetDKbDSigma(Mathematics::eval_pt	s) const
		{
			double ival_center = (s.left + s.right) * 0.5;
#ifdef COUNT_OPS
			g_nadd += 1;
			g_nmult += 1;
#endif
			return CTR::get_derivative(m_kb)(s.s)*m_chi(ival_center);
		}

		/**
			\brief Get the derivative of torsional compliance
		*/
		double GetDCtDSigma(Mathematics::eval_pt	s) const
		{
			double ival_center = (s.left + s.right) * 0.5;
#ifdef COUNT_OPS
			g_nadd += 1;
			g_nmult += 1;
#endif
			return CTR::get_derivative(m_ct)(s.s)*m_chi(ival_center);
		}

		const CurvFunType& GetCurvatureFun() const
		{
			return m_ustar;
		}

		const KbFunType& GetKbFun() const
		{
			return m_kb;
		}

		const CtFunType& GetCtFun() const
		{
			return m_ct;
		}

		const Functions::indicator_function& GetMaterialIndicator() const
		{
			return m_chi;
		}

		const Functions::mollified_indicator_fun<double>& GetCurvatureIndicator() const
		{
			return m_uChi;
		}

		Tube<CurvFunType, KbFunType, CtFunType>& operator=(Tube<CurvFunType, KbFunType, CtFunType>&& other)
		{
			m_L = other.m_L;
			m_tL = other.m_tL;
			m_ustar = other.m_ustar;
			m_kb = other.m_kb;
			m_ct = other.m_ct;
			m_chi = other.m_chi;
			m_uChi = other.m_uChi;
			return *this;
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	private:
		double m_L;
		double m_tL;
		CurvFunType	m_ustar;
		KbFunType m_kb;
		CtFunType m_ct;
		Functions::indicator_function				      m_chi;  //!< Used to mask the stiffnesses (discontinuous)
		Functions::mollified_indicator_fun<double>	m_uChi; //!< Used to mask the curvature (continuous)
      Functions::mollified_indicator_fun_der<double>  m_DuChi;
      typedef typename std::result_of< internal::derivative_helper( const CurvFunType& ) >::type   CurvDerType;
      CurvDerType m_Dustar;
	};
	
	/**
		\brief Makes a Tube from the given parameters (convenience function)

		\param[in] L	The total length of the tube in meters
		\param[in] OD	The outer diameter of the annular tube in meters
		\param[in] ID	The inner diameter of the annular tube in meters
		\param[in] k	The curvature function (unary function object returning a 2d vector)
		\param[in] E	The Young's Modulus of the material of the tube, so that E*I is the bending stiffness
		\param[in] G	The Shear Modulus of the material of the tube, so that G*J is the torsional stiffness

		\return The tube, with type Tube<CurvFunType>
	*/
	template <typename CurvFunType>
	Tube<CurvFunType>	
	make_annular_tube(double L, double Lt, double OD, double ID, const CurvFunType& k, double E, double G)
	{
		return Tube<CurvFunType>(L, 
								 Lt,
								 k, 
								 Functions::constant_fun<double>(AnnularBendingStiffness(OD, ID, E)), 
								 Functions::constant_fun<double>(AnnularTorsionalCompliance(OD, ID, G)) 
								 );
	}

	
}