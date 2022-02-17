/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once
#include <boost/mpl/bool.hpp>
#include <utility>
#include <Eigen/Geometry>
#include "CTRTypedefs.h"
#include "BasicFunctions.h"
#include "TupleTransform.h"
#include "TupleAccumulate.h"
#include "BinaryOperations.h"
#include "RungeKutta.h"
#include "EvalPt.h"
#include <limits>

#define NTUBES( X ) ((internal::traits<X>::N::value))
#define NTUBESC( X ) ((X::N::value))

namespace CTR
{
   namespace internal
   {
      template <typename State, typename Cannula>
      struct evaluation_context;
   }

   namespace internal {
      struct evaluate_precurvature
      {
         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef typename Vector<2>::type	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_precurvature( Tube, double )>::type operator()( const Tube& t, double sigma ) const
         {
            return t.GetPrecurvature( sigma );
         }

      private:
      };

      struct evaluate_dprecurvature
      {
         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef typename Vector<2>::type	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_dprecurvature( Tube, double )>::type operator()( const Tube& t, double sigma ) const
         {
            return t.GetDPrecurvature( sigma );
         }

      private:
      };

      struct evaluate_length
      {
         evaluate_length()
         {}

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube ) >
         {
            typedef double type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_length( Tube )>::type	operator()( const Tube& t ) const
         {
            return t.GetLength();
         }
      };

      struct evaluate_precurvature_at
      {
         evaluate_precurvature_at( double s ) : s_( s )
         {}

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube ) >
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result< evaluate_precurvature_at( Tube ) >::type  operator()( const Tube& t ) const
         {
            return t.GetPrecurvature( s_ );
         }

      private:
         double s_;
      };
   }

   namespace internal {
      struct evaluate_kb_at
      {
         evaluate_kb_at( Mathematics::eval_pt s ) : s_( s )
         {}

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_kb_at( Tube )>::type operator()( const Tube& t ) const
         {
            return t.GetBendingStiffness( s_ );
         }

      private:
         Mathematics::eval_pt s_;
      };

      struct evaluate_kb
      {

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_kb( Tube, double )>::type operator()( const Tube& t, Mathematics::eval_pt sigma ) const
         {
            return t.GetBendingStiffness( sigma );
         }
      };

      struct evaluate_dkb
      {

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_kb( Tube, double )>::type operator()( const Tube& t, Mathematics::eval_pt sigma ) const
         {
            return t.GetDKbDSigma( sigma );
         }
      };
   }

   namespace internal
   {
      struct evaluate_ct_at
      {
         evaluate_ct_at( Mathematics::eval_pt s ) : s_( s )
         {}

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_ct_at( Tube )>::type operator()( const Tube& t ) const
         {
            return t.GetTorsionalCompliance( s_ );
         }

      private:
         Mathematics::eval_pt s_;
      };

      struct evaluate_ct
      {

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_ct( Tube, double )>::type operator()( const Tube& t, Mathematics::eval_pt sigma ) const
         {
            return t.GetTorsionalCompliance( sigma );
         }

      };

      struct evaluate_dct
      {

         template <typename Tube>
         struct result;

         template <typename F, typename Tube>
         struct result < F( Tube, double ) >
         {
            typedef double	type;
         };

         template <typename Tube>
         STRONG_INLINE
            typename result<evaluate_ct( Tube, double )>::type operator()( const Tube& t, Mathematics::eval_pt sigma ) const
         {
            return t.GetDCtDSigma( sigma );
         }

      };
   }

   namespace internal
   {
      struct rotation_about_z
      {
         template <class>
         struct result;

         template <class F, class Vector, typename Scalar>
         struct result<F( Vector, Scalar )>
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <class F, class Vector>
         struct result<F( Vector, Vector )>
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <typename Vector>
         STRONG_INLINE
            typename result<rotation_about_z( Vector, double )>::type operator()( Vector const& v, double a ) const
         {
#ifdef COUNT_OPS
            g_ntrig += 2;
            g_nmult += 4;
            g_nadd += 2;
#endif
            double sa = sin( a ); double ca = cos( a );
            typename CTR::Vector<2>::type	r;
            r( 0 ) = ca*v( 0 ) - sa*v( 1 );
            r( 1 ) = sa*v( 0 ) + ca*v( 1 );
            return r;
         }

         template <typename Vector1, typename Vector2>
         STRONG_INLINE
            typename result<rotation_about_z( Vector1, Vector2 )>::type operator()( Vector1 const& v, Vector2 const& sincos ) const
         {
            //g_ntrig += 2;
#ifdef COUNT_OPS
            g_nmult += 4;
            g_nadd += 2;
#endif
            double sa = sincos( 0 ); double ca = sincos( 1 );
            typename CTR::Vector<2>::type	r;
            r( 0 ) = ca*v( 0 ) - sa*v( 1 );
            r( 1 ) = sa*v( 0 ) + ca*v( 1 );
            return r;
         }

      };
   }

   namespace internal
   {
      struct d_rotation_about_z
      {
         template <class>
         struct result;

         template <class F, class Vector, typename Scalar>
         struct result<F( Vector, Scalar )>
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <typename Vector>
         STRONG_INLINE
            typename result<d_rotation_about_z( Vector, double )>::type operator()( Vector const& v, double a ) const
         {
#ifdef COUNT_OPS
            g_ntrig += 2;
            g_nmult += 4;
            g_nadd += 2;
#endif
            double dsa = cos( a ); double dca = -sin( a );
            typename CTR::Vector<2>::type	dr;
            dr( 0 ) = dca*v( 0 ) - dsa*v( 1 );
            dr( 1 ) = dsa*v( 0 ) + dca*v( 1 );
            return dr;
         }

         template <typename Vector1, typename Vector2>
         STRONG_INLINE
            typename result<d_rotation_about_z( Vector1, Vector2 )>::type operator()( Vector1 const& v, Vector2 const& sincos ) const
         {
            //g_ntrig += 2;
#ifdef COUNT_OPS
            g_nmult += 4;
            g_nadd += 2;
#endif
            double dsa = sincos( 1 ); double dca = -sincos( 0 );
            typename CTR::Vector<2>::type	dr;
            dr( 0 ) = dca*v( 0 ) - dsa*v( 1 );
            dr( 1 ) = dsa*v( 0 ) + dca*v( 1 );
            return dr;
         }
      };
   }

   namespace internal
   {
      struct dd_rotation_about_z
      {
         template <class>
         struct result;

         template <class F, class Vector, typename Scalar>
         struct result<F( Vector, Scalar )>
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <class F, class Vector>
         struct result<F( Vector, Vector )>
         {
            typedef typename CTR::Vector<2>::type	type;
         };

         template <typename Vector>
         STRONG_INLINE
            typename result<d_rotation_about_z( Vector, double )>::type operator()( Vector const& v, double a ) const
         {
#ifdef COUNT_OPS
            g_ntrig += 2;
            g_nmult += 4;
            g_nadd += 2;
#endif
            double ddsa = -sin( a ); double ddca = -cos( a );
            typename CTR::Vector<2>::type	ddr;
            ddr( 0 ) = ddca*v( 0 ) - ddsa*v( 1 );
            ddr( 1 ) = ddsa*v( 0 ) + ddca*v( 1 );
            return ddr;
         }

         template <typename Vector1, typename Vector2>
         STRONG_INLINE
            typename result<d_rotation_about_z( Vector1, Vector2 )>::type operator()( Vector1 const& v, Vector2 const& sincos ) const
         {
            //g_ntrig += 2;
#ifdef COUNT_OPS
            g_nmult += 4;
            g_nadd += 2;
#endif
            double ddsa = -sincos( 0 ); double ddca = -sincos( 1 );
            typename CTR::Vector<2>::type	ddr;
            ddr( 0 ) = ddca*v( 0 ) - ddsa*v( 1 );
            ddr( 1 ) = ddsa*v( 0 ) + ddca*v( 1 );
            return ddr;
         }
      };
   }

   namespace internal
   {
      template <typename Vector1>
      struct dot_with_vector
      {
         dot_with_vector( Vector1 const& v_ ) : v( v_ )
         {}

         template <typename T>
         struct result;

         template <typename F, typename Vector2>
         struct result<F( Vector2 )>
         {
            typedef decltype(std::declval<Vector1>().dot( std::declval<Vector2>() ))	type;
         };

         template <typename Vector2>
         STRONG_INLINE
            typename result < dot_with_vector<Vector1>( Vector2 ) >::type
            operator()( Vector2 const& v2 ) const
         {
#ifdef COUNT_OPS
            g_nmult += Vector1::RowsAtCompileTime;
            g_nadd += Vector1::RowsAtCompileTime - 1;
#endif
            return v.dot( v2 );
         }

      private:
         Vector1 const& v;

      };
   }

   namespace internal
   {
      template <typename Context>
      STRONG_INLINE
         typename Matrix< 2, Context::N::value >::type
         compute_sines_and_cosines( const Context &c )
      {
         typename Matrix<2, Context::N::value >::type	sincos;
         sincos.template topRows<1>() = c.y.Psi.array().sin();
         sincos.template bottomRows<1>() = c.y.Psi.array().cos();
#ifdef COUNT_OPS
         g_ntrig += 2 * c.y.Psi.size();
#endif
         return sincos;
      }

      template <typename Context>
      STRONG_INLINE
         typename Matrix< 2, NTUBESC( Context ) >::type
         get_precurvatures( const Context &c )
      {
         using namespace tuple_transform;
         typename Matrix<2, NTUBESC( Context ) >::type	u_star;
			auto t = transform(c.C, c.y.Sigma, internal::evaluate_precurvature());
         copy( t, columns_of( u_star ) );
         return u_star;
      }

      template <typename Context>
      typename Matrix< 2, NTUBESC( Context ) >::type
         get_F_matrix( const Context &c )
      {
         using namespace tuple_transform;
         typename Matrix<2, NTUBESC( Context ) >::type F;
         copy( transform( columns_of( c.u_star ), columns_of( c.sincos ), internal::rotation_about_z() ), columns_of( F ) );
         return F;
      }

      template <typename Context>
      typename Matrix< 2, NTUBESC( Context ) >::type
         get_G_matrix( const Context &c )
      {
         using namespace tuple_transform;
         typename Matrix<2, NTUBESC( Context ) >::type G;
         copy( transform( columns_of( c.z_context.du_star ), columns_of( c.sincos ), internal::rotation_about_z() ), columns_of( G ) );
         return G;
      }

      template <typename Context>
      typename Matrix< 2, NTUBESC( Context ) >::type
         get_S_matrix( const Context &c )
      {
         using namespace tuple_transform;
         typename Matrix<2, NTUBESC( Context ) >::type S;
         copy( transform( columns_of( c.z_context.du_star ), columns_of( c.sincos ), internal::d_rotation_about_z() ), columns_of( S ) );
         return S;
      }


      template <typename Context>
      typename Matrix< 2, NTUBESC( Context ) >::type
         get_M_matrix( const Context &c )
      {
         using namespace tuple_transform;
         typename Matrix<2, NTUBESC( Context ) >::type	M;
         copy( transform( columns_of( c.u_star ), columns_of( c.sincos ), internal::d_rotation_about_z() ), columns_of( M ) );
         return M;
      }

      struct make_sigma_eval_pt
      {
         make_sigma_eval_pt( std::tuple<double, double> const& offsets_ ) :
            offsets( offsets_ )
         {}

         Mathematics::eval_pt	operator()( double sigma ) const
         {
#ifdef COUNT_OPS
            g_nadd += 2;
#endif
            return Mathematics::eval_pt( sigma, sigma + std::get<0>( offsets ), sigma + std::get<1>( offsets ) );
         }

         std::tuple<double, double>	offsets;
      };

      template <typename State, typename Cannula>
      typename Vector<NTUBES( State ) >::type
         get_bending_stiffnesses( const Mathematics::eval_pt &s, const State &y, const Cannula &C )
      {
         using namespace tuple_transform;
#ifdef COUNT_OPS
         g_nadd += 2;
#endif
         std::tuple<double, double> offsets = std::make_tuple( s.left - s.s, s.right - s.s );
         auto sigma_pts = transform( y.Sigma, make_sigma_eval_pt( offsets ) );
         auto kb_list = transform( C, sigma_pts, internal::evaluate_kb() );
         typename Vector<NTUBES( State ) >::type	kb;
         copy( kb_list, kb );
         return kb;
      }

      template < typename State, typename Cannula >
      typename Vector<NTUBES( State )>::type
         get_d_bending_stiffnesses( const Mathematics::eval_pt &s, const State &y, const Cannula &C )
      {
         using namespace tuple_transform;
#ifdef COUNT_OPS
         g_nadd += 2;
#endif
         std::tuple<double, double> offsets = std::make_tuple( s.left - s.s, s.right - s.s );
         auto sigma_pts = transform( y.Sigma, make_sigma_eval_pt( offsets ) );
         auto dkb_list = transform( C, sigma_pts, internal::evaluate_dkb() );
         typename Vector<NTUBES( State )>::type	dkb;
         copy( dkb_list, dkb );
         return dkb;
      }

      template <typename Context >
      typename Vector<NTUBESC( Context ) >::type
         get_torsional_compliances( const Context &c )
      {
         using namespace tuple_transform;
#ifdef COUNT_OPS
         g_nadd += 2;
#endif
         std::tuple<double, double> offsets = std::make_tuple( c.s.left - c.s.s, c.s.right - c.s.s );
         auto sigma_pts = transform( c.y.Sigma, make_sigma_eval_pt( offsets ) );
         auto ct_list = transform( c.C, sigma_pts, internal::evaluate_ct() );
         typename Vector<NTUBESC( Context ) >::type	ct;
         copy( ct_list, ct );
         return ct;
      }

      template <typename Context >
      typename Vector<NTUBESC( Context ) >::type
         get_d_torsional_compliances( const Context &c )
      {
         using namespace tuple_transform;
#ifdef COUNT_OPS
         g_nadd += 2;
#endif
         std::tuple<double, double> offsets = std::make_tuple( c.s.left - c.s.s, c.s.right - c.s.s );
         auto sigma_pts = transform( c.y.Sigma, make_sigma_eval_pt( offsets ) );
         auto dct_list = transform( c.C, sigma_pts, internal::evaluate_dct() );
         typename Vector<NTUBESC( Context ) >::type	dct;
         copy( dct_list, dct );
         return dct;
      }


      template <typename Context>
      typename Vector<2>::type get_bishop_curvature( const Context &c )
      {
         typedef typename internal::traits< typename Context::State >::ExternalLoads	LoadOption;

         return get_bishop_curvature_impl( c, LoadOption() );
      }

      template <typename Context>
      typename Vector<2>::type get_bishop_curvature_impl( const Context &c, boost::mpl::false_ no_loads )
      {
         auto &kb = c.kb;
         auto &F = c.F;

         //The bishop curvature is a weighted sum of the columns of F, which consists of
         // the precurvature vectors of each tube, rotated into the common Bishop frame
#ifdef COUNT_OPS
         g_nmult += F.rows() * kb.cols() * F.cols();
         g_nadd += F.rows() * kb.cols() * (F.cols() - 1);
         g_nmult += F.rows();
#endif
         return c.kb_tot_inv*(F*kb);
      }

      template <typename Context>
      typename Vector<2>::type get_bishop_curvature_impl( const Context &c, boost::mpl::true_ has_loads )
      {
         using namespace tuple_transform;

         auto &kb = c.kb;
         auto &F = c.F;

         typename Vector<3>::type   Mb = Utility::quaternion_rotate( Utility::quaternion_conj( c.g_context.q_renorm ), c.y.Mt );

         //The bishop curvature is a weighted sum of the columns of F, which consists of
         // the precurvature vectors of each tube, rotated into the common Bishop frame
         // plus the effect of the body frame moment carried by the collection of tubes
#ifdef COUNT_OPS
         g_nmult += F.rows() * kb.cols() * F.cols();
         g_nadd += F.rows() * kb.cols() * (F.cols() - 1) + F.rows();
         g_nmult += 2 * F.rows();
#endif
         return c.kb_tot_inv*(F*kb) + c.kb_tot_inv * Mb.topRows<2>();
      }

      template <typename Context>
      typename Vector< NTUBESC( Context ) >::type
         get_angle_derivatives( const Context &c )
      {
         //coefficient-wise product between 
         //torsional compliances and torsional moments gives
         //the angular rate of change for each tube
         typename Vector<NTUBESC( Context ) >::type ct = internal::get_torsional_compliances( c );
#ifdef COUNT_OPS
         g_nmult += NTUBESC( Context );
#endif
         return c.y.Mz.cwiseProduct( ct );
      }

      template <typename Context>
      typename Vector< NTUBESC( Context ) >::type	get_moment_derivatives( const Context &c )
      {
         using namespace tuple_transform;
         const static int N = NTUBESC( Context );

         //Define cross-product matrix of e3 (top left corner)
         typename Matrix<2, 2>::type	e3hat = Matrix<2, 2>::type::Zero();
         e3hat( 0, 1 ) = -1.0; e3hat( 1, 0 ) = 1.0;

         //Compute the derivatives of the torsional moments (vectorized)
#ifdef COUNT_OPS
         g_nmult += 4 + 3 * NTUBESC( Context );
         g_nadd += 2 + N;
#endif
         return (c.kb.asDiagonal()*(c.F.transpose()*(e3hat*c.uB)));

      }

      template <typename Context, typename State>
      void xi_fun_impl_t( const Context &c, State& yprime )
      {
         //compute the sigma derivatives
         yprime.Sigma = Vector< NTUBESC( Context ) >::type::Ones();

         //compute the angle derivatives	
         yprime.Psi = get_angle_derivatives( c );

         //compute the moment derivatives
         yprime.Mz = get_moment_derivatives( c );
      }

      template <typename Context, typename State>
      void xi_fun_impl_g_impl( const Context &c, State& yprime, boost::mpl::true_ )
      {
         //Get the derivative of the position (body-frame)
         yprime.p = Utility::quaternion_zaxis( c.g_context.q_renorm );

         //Get the derivative of the orientation (body-frame)
         typename Vector<3>::type    qUb = Functions::zero_fun< typename Vector<3>::type >()();
         qUb.topRows<2>() = c.uB;
         yprime.q = 0.5 * Utility::quaternion_mult_vec3( c.g_context.q_renorm, qUb );

#ifdef COUNT_OPS
         g_nmult += 4;
#endif

      }

      template <typename Context, typename State>
      void xi_fun_impl_g_impl( const Context& c, State& yprime, boost::mpl::false_ )
      {}

      template <typename Context, typename State>
      void xi_fun_impl_g( Context const& c, State& yprime )
      {
         //Figure out if we actually need to compute the geometry, and if so
         // delegate it to the appropriate function call below
         typedef typename internal::traits< typename Context::State >::Options	OType;
         typedef typename Option::has_option<OType, Option::ComputeGeometry>::type	ComputeGeometry;

         //Dispatch via function overload to either an empty implementation
         // or the function which performs the computation for dp/ds and dq/ds
         xi_fun_impl_g_impl( c, yprime, ComputeGeometry() );
      }

      template <typename Context>
      typename Vector<3>::type xi_n( const Context &c )
      {
         //Fill out 3D Bishop curvature vector
         //typename Vector<3>::type	uB = typename Vector<3>::type::Zero();
         //uB.topRows<2>() = c.uB;

         //return -uB.cross(c.y.Nb);
         return Functions::zero_fun< typename Vector<3>::type >()();
      }

      template <typename Context>
      typename Vector<3>::type   xi_m( const Context &c )
      {
         //Get the moment vector by pulling the x- and y- components
         // from the state and the z component from the sum of mz
         //typename Vector<3>::type	Mb = typename Vector<3>::type::Zero();
         //Mt.topRows<2>() = c.y.Mb;
         //M(2) = c.y.Mz.sum();

         ////Define third standard basis vector
         //typename Vector<3>::type	e3 = typename Vector<3>::type::UnitZ();

         ////Fill out 3D Bishop curvature vector
         //typename Vector<3>::type uB = typename Vector<3>::type::Zero();
         //uB.topRows<2>() = c.uB;

         //return (-uB.cross(Mb) - e3.cross(c.y.Nb)).topRows<2>();

         typename Vector<3>::type pdot = Utility::quaternion_zaxis( c.g_context.q_renorm );
#ifdef COUNT_OPS
         g_nmult += 6;
         g_nadd += 3;
#endif
         return c.y.Nt.cross( pdot );
      }

      template <typename Context, typename State>
      void xi_fun_impl_l_impl( const Context &c, State &yprime, boost::mpl::true_ has_loads )
      {
         yprime.Nt = xi_n( c );
         yprime.Mt = xi_m( c );
      }

      template <typename Context, typename State>
      void xi_fun_impl_l_impl( const Context &c, State &yprime, boost::mpl::false_ no_loads )
      {}

      template <typename Context, typename State>
      void xi_fun_impl_l( const Context &c, State& yprime )
      {
         typedef typename internal::traits<State>::Options	OType;
         typedef typename Option::has_option<OType, Option::ExternalLoads>::type	ExternalLoads;

         xi_fun_impl_l_impl( c, yprime, ExternalLoads() );
      }



      template <typename Context>
      typename Matrix< NTUBESC( Context ), NTUBESC( Context ) >::type
         d_xi_psi_d_sigma( const Context &c )
      {
         typename Vector<NTUBESC( Context ) >::type	dct = internal::get_d_torsional_compliances( c );
#ifdef COUNT_OPS
         g_nmult += NTUBESC( Context );
#endif
         return (dct.cwiseProduct( c.y.Mz )).asDiagonal();
      }

      template <typename Context>
      typename Matrix< NTUBESC( Context ), NTUBESC( Context ) >::type
         d_xi_psi_d_mz( const Context &c )
      {
         auto ct_v = internal::get_torsional_compliances( c );
         return ct_v.asDiagonal();
      }

      template <typename Context>
      typename Matrix<2, NTUBESC( Context ) >::type
         d_ub_d_psi( const Context &c )
      {
         using namespace tuple_transform;

         typename Matrix<2, 2>::type    J = Matrix<2, 2>::type::Zero();
         J( 0, 1 ) = -1.0; J( 1, 0 ) = 1.0;
#ifdef COUNT_OPS
         g_nmult += 6 * NTUBESC( Context );
         g_nadd += 2 * NTUBESC( Context );
#endif
         return  J*c.F * (c.kb_tot_inv*c.kb).asDiagonal();
      }

      template <typename Context>
      typename Matrix<2, 3>::type
         d_ub_d_q( const Context &c )
      {
         typename Matrix<3, 3>::type R = Utility::quaternion_to_rotation_matrix( c.g_context.q_renorm );
#ifdef COUNT_OPS
         g_nmult += 4 + 9;
         g_nadd += 6;
#endif
         return (c.kb_tot_inv*Utility::cpm3( R.transpose()*c.y.Mt )).template topRows<2>();
      }

      template <typename Context>
      typename Matrix< 2, NTUBESC( Context ) >::type
         d_ub_d_sigma( const Context &c )
      {
         const int N = NTUBESC( Context );

         using namespace tuple_transform;
         typename Vector<N>::type	const &kb = c.kb;
         typename Vector<N>::type	const &dkb = c.z_context.dkb;

         typename Matrix<2, N>::type	const &F = c.F;
         typename Matrix<2, N>::type G = internal::get_G_matrix( c );
         typename Vector<2>::type	const &uB = c.uB;
#ifdef COUNT_OPS
         g_nmult += 2 * NTUBESC( Context ) * 1 + F.size() + G.size();
#endif
         return -uB*(c.kb_tot_inv*dkb).transpose()
            + F*(c.kb_tot_inv*dkb).asDiagonal()
            + G*(c.kb_tot_inv*kb).asDiagonal();
      }

      template <typename Context>
      typename Matrix< 2, 3 >::type
         d_ub_d_mt( const Context &c )
      {
         typedef typename Matrix<2, 2>::type	MType;

         typename Matrix<3, 3>::type R = Utility::quaternion_to_rotation_matrix( c.g_context.q_renorm );
#ifdef COUNT_OPS
         g_nmult += 6;
#endif
         return c.kb_tot_inv * R.transpose().topRows<2>();
      }

      template <typename Context>
      typename Matrix< NTUBESC( Context ), 3 >::type
         d_xi_mz_d_mt( const Context &c )
      {
         using namespace tuple_transform;
         const static int N = NTUBESC( Context );

         typename Matrix<2, 2>::type	e3hat = Matrix<2, 2>::type::Zero();
         e3hat( 0, 1 ) = -1.0; e3hat( 1, 0 ) = 1.0;
#ifdef COUNT_OPS
         g_nmult += 2 * 2 * 3 + 2 * NTUBESC( Context ) + NTUBESC( Context ) * 2 * 3;
         g_nadd += 3 * NTUBESC( Context ) + 6;
#endif
         return (c.kb.asDiagonal()*c.F.transpose())*(e3hat*d_ub_d_mt( c ));
      }



      template <typename Context>
      typename Matrix< NTUBESC( Context ), NTUBESC( Context ) >::type
         g_d_sigma( const Context &c )
      {
         using namespace tuple_transform;
         const int N = NTUBESC( Context );

         typename Vector<N>::type	const &kb = c.kb;
         typename Vector<N>::type	const &dkb = c.z_context.dkb;
         typename Vector<2>::type	const &uB = c.uB;

         typename Matrix<2, N>::type   G = internal::get_G_matrix( c );
         typename Matrix<2, 2>::type	J = Matrix<2, 2>::type::Zero();
         J( 0, 1 ) = -1.0; J( 1, 0 ) = 1.0;

         typename Matrix<2, N>::type	d_ub = c.z_context.dub_d_sigma;
#ifdef COUNT_OPS
         g_nadd += 2 + NTUBESC( Context );
         g_nmult += 3 * NTUBESC( Context ) + 4;
#endif
         typename Matrix<N, N>::type	A = dkb.cwiseProduct( c.F.transpose()*(J*uB) ).asDiagonal();

#ifdef COUNT_OPS
         g_nadd += NTUBESC( Context )*NTUBESC( Context ) + 2 * NTUBESC( Context );
         g_nmult += 3 * NTUBESC( Context )*NTUBESC( Context ) + 4 * NTUBESC( Context );
#endif
         typename Matrix<N, N>::type	B = (kb.asDiagonal()) * (c.F.transpose()*(J*d_ub));

#ifdef COUNT_OPS
         g_nadd += 2 + NTUBESC( Context );
         g_nmult += 3 * NTUBESC( Context ) + 4;
#endif
         typename Matrix<N, N>::type	D = kb.cwiseProduct( G.transpose()*J*uB ).asDiagonal();

#ifdef COUNT_OPS
         g_nadd += A.size() + B.size();
#endif

         return (A + B + D);
      }

      template <typename Context>
      typename Matrix< NTUBESC( Context ), NTUBESC( Context ) >::type
         g_d_psi( const Context &c )
      {
         using namespace tuple_transform;
         const int N = NTUBESC( Context );

         typename Vector<N>::type const &kb = c.kb;

         typename Matrix<2, 2>::type J = Matrix<2, 2>::type::Zero();
         J( 0, 1 ) = -1.0;
         J( 1, 0 ) = 1.0;

         typename Vector<2>::type	const &uB = c.uB;
         typename Matrix<2, N>::type	d_ub = c.z_context.dub_d_psi;
#ifdef COUNT_OPS
         g_nmult += 4 * NTUBESC( Context )*NTUBESC( Context ) + 7 * NTUBESC( Context );
         g_nadd += 3 * NTUBESC( Context ) *NTUBESC( Context ) + NTUBESC( Context );
#endif
         return (kb.asDiagonal()*c.F.transpose())*(J*d_ub) + (kb.cwiseProduct( c.F.transpose()*uB ).asDiagonal().toDenseMatrix());
      }

      template <typename Context>
      typename Matrix<2 * NTUBESC( Context ), NTUBESC( Context )>::type
         xi_t_sigma( const Context &c )
      {
         //Layout of xi_ts
         // [dxi_psi/dsigma]
         // [dxi_mz/dsigma ]
         const int N = NTUBESC( Context );
         typename Matrix<2 * N, N>::type	xi_t_s = Matrix<2 * N, N>::type::Zero();

         //Set the block dxi_psi/dsigma
         xi_t_s.template block<N, N>( 0, 0 ) = d_xi_psi_d_sigma( c );

         //Set the block dxi_mz/dsigma
         xi_t_s.template block<N, N>( N, 0 ) = g_d_sigma( c );

         return xi_t_s;
      }

      template <typename Context>
      typename Matrix<2 * NTUBESC( Context ), 2 * NTUBESC( Context )>::type
         xi_t_t( const Context &c )
      {
         //Layout of xi_tt
         // [dxi_psi/dpsi   dxi_psi/dmz ]
         // [dxi_mz/dpsi    dxi_mz/dmz  ]
         const int N = NTUBESC( Context );
         typename Matrix<2 * N, 2 * N>::type	xi_t_t = Matrix<2 * N, 2 * N>::type::Zero();

         //Set the block dxi_psi/dmz
         xi_t_t.template block<N, N>( 0, N ) = d_xi_psi_d_mz( c );

         //Set the block dxi_mz/dpsi
         xi_t_t.template block<N, N>( N, 0 ) = g_d_psi( c );

         return xi_t_t;
      }

      template <typename Context>
      typename Matrix<2 * NTUBESC( Context ), 6>::type
         xi_t_g( const Context &c )
      {
         //Layout of xi_tg
         // [dxi_psi/dp    dxi_psi/dR ]
         // [dxi_mz/dp     dxi_mz/dR  ]
         const int N = NTUBESC( Context );
         typename Matrix<2 * N, 6>::type	xi_tg = Matrix<2 * N, 6>::type::Zero();

         //Set the block dxi_mz/dR
         typename Matrix<2, 2>::type	e3hat = Matrix<2, 2>::type::Zero();
         e3hat( 0, 1 ) = -1.0; e3hat( 1, 0 ) = 1.0;
#ifdef COUNT_OPS
         g_nadd += 3 * NTUBESC( Context ) + 6;
         g_nmult += 12 + 9 * NTUBESC( Context );
#endif
         xi_tg.template bottomRightCorner<N, 3>() = c.kb.asDiagonal() * c.F.transpose() * e3hat * d_ub_d_q( c );

         return xi_tg;
      }

      template <typename Context>
      typename Matrix<6, NTUBESC( Context ) >::type
         xi_g_sigma( const Context &c )
      {
         //Layout of xi_gs
         // [dxi_p/dsigma]
         // [dxi_R/dsigma]
         const int N = NTUBESC( Context );
         typedef typename Matrix<6, NTUBESC( Context ) >::type	MType;

         MType xi_gs = Functions::zero_fun< MType >()();
         //Set the block dxi_R/dsigma
         xi_gs.template block<2, N>( 3, 0 ) = c.z_context.dub_d_sigma;
         return xi_gs;
      }

      template <typename Context>
      typename Matrix<6, 2 * NTUBESC( Context ) >::type
         xi_g_t( const Context &c )
      {
         //Layout of xi_gt
         // [dxi_p/dpsi   dxi_p/dmz ]
         // [dxi_R/dpsi   dxi_R/dmz ]
         //
         const int N = NTUBESC( Context );
         typedef typename Matrix<6, 2 * NTUBESC( Context ) >::type	MType;

         MType xi_gt = Functions::zero_fun< MType >()();

         //Set the block dxi_R/dpsi
         xi_gt.template block<2, N>( 3, 0 ) = c.z_context.dub_d_psi;

         return xi_gt;
      }

      template <typename Context>
      typename Matrix<6, 6 >::type
         xi_g_l( const Context &c )
      {
         //Layout of xi_gl
         //  [dxi_p/dnt     dxi_p/dmt ]
         //  [dxi_R/dnt     dxi_R/dmt ]
         //
         const int N = NTUBESC( Context );
         typedef typename Matrix<6, 6 >::type	MType;

         MType xi_gl = Functions::zero_fun< MType >()();

         //Set the block dxi_R/dmxy
         xi_gl.block<2, 3>( 3, 3 ) = d_ub_d_mt( c );

         return xi_gl;
      }

      template <typename Context>
      typename Matrix<6, 6 >::type
         xi_g_g( const Context &c )
      {
         //Layout of xi_gl
         //  [dxi_p/dp     dxi_p/dR ]
         //  [dxi_R/dp     dxi_R/dR ]
         //
         const int N = NTUBESC( Context );
         typedef typename Matrix<6, 6 >::type	MType;

         MType xi_gg = Functions::zero_fun< MType >()();

         //Set the block dxi_R/dR
         xi_gg.block<2, 3>( 3, 3 ) = d_ub_d_q( c );

         return xi_gg;
      }

      template <typename Context>
      typename Matrix<2 * NTUBESC( Context ), 6 >::type
         xi_t_l( const Context &c )
      {
         //Layout of xi_tl
         //  [dxi_psi/dnb    dxi_psi/dmxy]
         //  [dxi_mz/dnb     dxi_mz/dmxy ]
         //
         const int N = NTUBESC( Context );
         typedef typename Matrix<2 * NTUBESC( Context ), 6 >::type	MType;

         MType	xi_tl = Functions::zero_fun< MType >()();

         //Set the block dxi_mz/dmxy
         xi_tl.template block<N, 3>( N, 3 ) = d_xi_mz_d_mt( c );

         return xi_tl;
      }

      template <typename Context>
      typename Matrix<6, 6>::type
         xi_l_g( const Context &c )
      {
         //Layout of xi_lg
         // [dxi_nt/dp   dxi_mt/dp ]
         // [dxi_nt/dR   dxi_mt/dR ]
         typedef typename Matrix<6, 6>::type	MType;
         MType xi_lg = Functions::zero_fun< MType >()();

         typename Matrix<3, 3>::type NHat = Utility::cpm3( c.y.Nt );
#ifdef COUNT_OPS
         g_nadd += 12;
         g_nmult += 18;
#endif
         //only block here is dxi_mt/dR
         xi_lg.bottomRightCorner<3, 3>().col( 0 ) = -NHat*Utility::quaternion_yaxis( c.g_context.q_renorm );
         xi_lg.bottomRightCorner<3, 3>().col( 1 ) = NHat*Utility::quaternion_xaxis( c.g_context.q_renorm );

         return xi_lg;
      }

      template <typename Context>
      typename Matrix< 6, 6>::type
         xi_l_l( const Context &c )
      {
         //Layout of xi_ll
         //  [dxi_nt/dnt      dxi_nt/dmt]
         //  [dxi_mt/dnt      dxi_mt/dmt ]
         //
         const int N = NTUBESC( Context );
         typedef typename Matrix<6, 6>::type	MType;

         MType xi_ll = Functions::zero_fun<MType >()();

         xi_ll.bottomLeftCorner<3, 3>() = -Utility::cpm3( Utility::quaternion_zaxis( c.g_context.q_renorm ) );

         return xi_ll;
      }
   }

   namespace internal
   {
      template <typename Context>
      typename Matrix< 6, 6 >::type	ad_SE3( const Context &c )
      {
         typename Matrix<6, 6>::type	ad = Functions::zero_fun< typename Matrix<6, 6>::type >()();
         ad.block<3, 3>( 0, 0 ) = Utility::cpm3_2( c.uB );
         ad.block<3, 3>( 3, 3 ) = Utility::cpm3_2( c.uB );
         typename Vector<3>::type	e3;
         e3( 0 ) = 0; e3( 1 ) = 0; e3( 2 ) = 1;
         ad.template block<3, 3>( 0, 3 ) = Utility::cpm3( e3 );
         return ad;
      }
   }

   namespace internal
   {
      template <typename Context>
      typename Matrix< NTUBESC( Context ), 2 * NTUBESC( Context ) >::type
         Z_sa( Context const& )
      {
         const int N = NTUBESC( Context );
         typedef typename Matrix< NTUBESC( Context ), 2 * NTUBESC( Context ) >::type	MType;
         MType Zsa;
         Zsa.template block<N, N>( 0, 0 ) = Eigen::Matrix<double, N, N>::Zero();
         Zsa.template block<N, N>( 0, N ) = -Eigen::Matrix<double, N, N>::Identity();
         return Zsa;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zta_impl( const Context &c, State& yprime, boost::mpl::true_, boost::mpl::false_ )
      {
         //Differential equation for Zta when Zla is not present
         const int N = NTUBESC( Context );

         typename Matrix<2 * N, 2 * N>::type	xi_tt = xi_t_t( c );
         typename Matrix<2 * N, N>::type		xi_ts = xi_t_sigma( c );
#ifdef COUNT_OPS
         g_nadd += xi_ts.rows()*Z_sa( c ).cols()*(xi_ts.cols() - 1) + xi_tt.rows()*c.y.Zta.cols()*(xi_tt.cols() - 1);
         g_nmult += xi_ts.rows()*Z_sa( c ).rows()*Z_sa( c ).cols() + xi_tt.rows()*c.y.Zta.rows()*c.y.Zta.cols();
#endif
         yprime.Zta = xi_ts*Z_sa( c ) + xi_tt*c.y.Zta;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zta_impl( const Context &c, State& yprime, boost::mpl::true_, boost::mpl::true_ )
      {
         //differential equation for Zta when Zla is present
         const int N = NTUBESC( Context );

         typename Matrix<2 * N, 2 * N>::type	xi_tt = xi_t_t( c );
         typename Matrix<2 * N, N>::type		xi_ts = xi_t_sigma( c );
         typename Matrix<2 * N, 6>::type		xi_tl = xi_t_l( c );
         typename Matrix<2 * N, 6>::type     xi_tg = xi_t_g( c );
#ifdef COUNT_OPS
         g_nmult += xi_ts.rows()*Z_sa( c ).cols()*Z_sa( c ).rows() + xi_tt.rows()*c.y.Zta.rows()*c.y.Zta.cols();
         g_nmult += xi_tl.rows()*c.y.Zla.cols()*c.y.Zla.rows() + xi_tg.rows()*xi_tg.cols()*c.y.Zga.cols();

         g_nadd += xi_ts.rows()*Z_sa( c ).cols()*(Z_sa( c ).rows() - 1) + xi_tt.rows()*c.y.Zta.cols()*(c.y.Zta.rows() - 1);
         g_nadd += xi_tl.rows()*c.y.Zla.cols()*(c.y.Zla.rows() - 1) + xi_tg.rows()*c.y.Zga.cols()*(xi_tg.cols() - 1);
#endif
         yprime.Zta = xi_ts*Z_sa( c ) + xi_tt*c.y.Zta + xi_tl*c.y.Zla + xi_tg*c.y.Zga;
      }

      template <typename Context, typename State, typename DoesNotMatter>
      void xi_fun_impl_Zta_impl( const Context &c, State& yprime, boost::mpl::false_, DoesNotMatter )
      {

      }

      template <typename Context, typename State>
      void xi_fun_impl_Zga_impl( const Context &c, State& yprime, boost::mpl::true_ has_Zga, boost::mpl::true_ has_Zla )
      {
         typename Matrix<6, 2 * NTUBES( State ) >::type	xi_gt = xi_g_t( c );
         typename Matrix<6, NTUBES( State ) >::type		xi_gs = xi_g_sigma( c );
         typename Matrix<6, 6 >::type					   xi_gl = xi_g_l( c );
         typename Matrix<6, 6>::type                  xi_gg = xi_g_g( c );
#ifdef COUNT_OPS
      //TODO: Finish this op counting
#endif
         yprime.Zga = xi_gs*Z_sa( c ) + xi_gt*c.y.Zta + xi_gl*c.y.Zla - internal::ad_SE3( c )*c.y.Zga + xi_gg*c.y.Zga;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zga_impl( const Context &c, State& yprime, boost::mpl::true_ has_Zga, boost::mpl::false_ no_Zla )
      {
         typename Matrix<6, 2 * NTUBES( State ) >::type	xi_gt = xi_g_t( c );
         typename Matrix<6, NTUBES( State ) >::type		xi_gs = xi_g_sigma( c );
         //yprime.Zga = xi_gs * Z_sa(c) + xi_gt*c.y.Zta - internal::ad_SE3(c)*c.y.Zga;
#ifdef COUNT_OPS
      //TODO: FINISH THIS OP COUNTING
#endif

         yprime.Zga = xi_gt*c.y.Zta - internal::ad_SE3( c )*c.y.Zga;
         yprime.Zga.template rightCols< NTUBES( State ) >() -= xi_gs;
      }

      template <typename Context, typename State, typename DoesNotMatter>
      void xi_fun_impl_Zga_impl( const Context &c, State& yprime, boost::mpl::false_ no_Zga, DoesNotMatter )
      {}

      template <typename Context, typename State>
      void xi_fun_impl_Zga( const Context &c, State& yprime )
      {
         typedef typename internal::traits<State>::Options	OType;
         xi_fun_impl_Zga_impl( c, yprime, Option::has_option<OType, Option::ComputeZga>(), Option::has_option<OType, Option::ComputeZla>() );
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zla( const Context &c, State &yprime, boost::mpl::false_ )
      {

      }

      template <typename Context, typename State>
      void xi_fun_impl_Zla( const Context &c, State &yprime, boost::mpl::true_ )
      {
         typename Matrix<6, 6>::type						   xi_ll = xi_l_l( c );
         typename Matrix<6, 6>::type		               xi_lg = xi_l_g( c );
#ifdef COUNT_OPS
         //TODO: FINISH THIS OP COUNTING
#endif
         yprime.Zla = xi_ll*c.y.Zla + xi_lg*c.y.Zga;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Ztf( const Context &c, State &yprime, boost::mpl::true_ )
      {
         typename Matrix<2 * NTUBESC( Context ), 2 * NTUBESC( Context ) >::type   xi_tt = xi_t_t( c );
         typename Matrix<2 * NTUBESC( Context ), 6 >::type                        xi_tl = xi_t_l( c );
         typename Matrix<2 * NTUBESC( Context ), 6 >::type                        xi_tg = xi_t_g( c );

#ifdef COUNT_OPS
       //TODO: FINISH THIS OP COUNTING
#endif
         yprime.Ztf = xi_tt*c.y.Ztf + xi_tl*c.y.Zlf + xi_tg*c.y.Zgf;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Ztf( const Context &c, State &yprime, boost::mpl::false_ )
      {

      }

      template <typename Context, typename State>
      void xi_fun_impl_Zlf( const Context &c, State &yprime, boost::mpl::true_ )
      {
         typename Matrix<6, 6>::type  xi_ll = xi_l_l( c );
         typename Matrix<6, 6>::type  xi_lg = xi_l_g( c );

#ifdef COUNT_OPS
       //TODO: FINISH THIS OP COUNTING
#endif
         yprime.Zlf = xi_ll*c.y.Zlf + xi_lg*c.y.Zgf;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zlf( const Context &c, State &yprime, boost::mpl::false_ )
      {

      }

      template <typename Context, typename State>
      void xi_fun_impl_Zgf( const Context &c, State &yprime, boost::mpl::true_ )
      {
         typename Matrix<6, 2 * NTUBESC( Context ) >::type   xi_gt = xi_g_t( c );
         typename Matrix<6, 6>::type                         xi_gl = xi_g_l( c );
         typename Matrix<6, 6>::type						 xi_gg = xi_g_g( c );
         typename Matrix<6, 6>::type                         ad = internal::ad_SE3( c );

#ifdef COUNT_OPS
       //TODO: FINISH THIS OP COUNTING
#endif
         yprime.Zgf = xi_gt*c.y.Ztf + xi_gl*c.y.Zlf + xi_gg*c.y.Zgf - ad*c.y.Zgf;
      }

      template <typename Context, typename State>
      void xi_fun_impl_Zgf( const Context &c, State &yprime, boost::mpl::false_ )
      {

      }

      template <typename Context, typename State>
      void xi_fun_impl_Z( const Context &c, State& yprime )
      {
         xi_fun_impl_Zta_impl( c, yprime, typename internal::traits<State>::ComputeZta(), typename internal::traits<State>::ComputeZla() );
         xi_fun_impl_Zga( c, yprime );
         xi_fun_impl_Zla( c, yprime, typename internal::traits<State>::ComputeZla() );

         xi_fun_impl_Ztf( c, yprime, typename internal::traits<State>::ComputeZtf() );
         xi_fun_impl_Zlf( c, yprime, typename internal::traits<State>::ComputeZlf() );
         xi_fun_impl_Zgf( c, yprime, typename internal::traits<State>::ComputeZgf() );
      }
   }

   namespace internal
   {
      template <typename State_, typename Cannula_>
      struct evaluation_context;

      template <typename State_, typename Cannula_>
      struct no_evaluation_context_Z
      {
         typedef typename traits<State_>::N	N;
         typedef State_		State;
         typedef Cannula_	Cannula;

         no_evaluation_context_Z( Mathematics::eval_pt const&, const State&, const Cannula& )
         {
            //nothing to do here
         }

         void construct( evaluation_context<State_, Cannula_> const& eval_context )
         {}
      };

      template <typename State_, typename Cannula_>
      struct evaluation_context_Z
      {
         typedef typename traits<State_>::N	N;
         typedef State_		State;
         typedef Cannula_	Cannula;

         evaluation_context_Z( Mathematics::eval_pt const& s, const State &y, const Cannula &C ) :
            dkb( internal::get_d_bending_stiffnesses( s, y, C ) )
         {
            using namespace tuple_transform;
            //get the derivative precurvature list
            copy( transform( C, y.Sigma, internal::evaluate_dprecurvature() ),
                  columns_of( du_star ) );

               //Can't fill out dub_d_sigma yet since we need the rest of the eval_context,
               //which isn't available in this constructor
         }

         void construct( evaluation_context<State_, Cannula_> const& eval_context )
         {
            dub_d_sigma = internal::d_ub_d_sigma( eval_context );
            dub_d_psi = internal::d_ub_d_psi( eval_context );
         }

         typename Matrix<2, N::value>::type	      du_star;
         typename Vector<N::value>::type		      dkb;
         typename Matrix<2, N::value>::type        dub_d_sigma;
         typename Matrix<2, N::value>::type        dub_d_psi;
      };

      template <typename State_, typename Cannula_>
      struct evaluation_context_G
      {
         typedef typename traits<State_>::N	N;
         typedef State_		State;
         typedef Cannula_	Cannula;

         evaluation_context_G( const State &y, const Cannula &C ) :
            q_renorm( (1.0 / y.q.norm()) * y.q )
         {

         }

         typename Vector<4>::type	q_renorm;
      };

      template <typename State_, typename Cannula_>
      struct no_evaluation_context_G
      {
         typedef typename traits<State_>::N	N;
         typedef State_		State;
         typedef Cannula_	Cannula;

         no_evaluation_context_G( const State &y, const Cannula &C )
         {

         }

      };

       /**
          \brief Evaluation context for differential equations

          This structure was created out of the observation via profiling that
          the compiler was unable to eliminate many common subexpressions. These
          subexpressions would be computed many times in a single call to the
          d.e. right hand side. To avoid this, the evaluation_context holds
          these subexpressions in a structure. Depending on whether derivative
          matrices are computed, some extra subexpressions may be cached as well.
       */
      template <typename State_, typename Cannula_>
      struct evaluation_context
      {
         typedef State_		State;
         typedef Cannula_	Cannula;

         typedef typename boost::mpl::if_< typename traits<State_>::AnyJac,
            evaluation_context_Z<State_, Cannula_>,
            no_evaluation_context_Z<State_, Cannula_> >::type	Z_context;

         typedef typename boost::mpl::if_< typename traits<State_>::ComputeGeometry,
            evaluation_context_G<State_, Cannula_>,
            no_evaluation_context_G<State_, Cannula_> >::type	G_context;

         typedef typename traits<State>::N	N;
         evaluation_context( Mathematics::eval_pt const& s_, const State &y_, const Cannula &C_ ) :
            s( s_ ),
            y( y_ ),
            C( C_ ),
            sincos( internal::compute_sines_and_cosines( *this ) ),
            kb( internal::get_bending_stiffnesses( s, y, C ) ),
            kb_tot_inv( 1.0 / kb.sum() ), //this is guaranteed to be nonzero by the stiffness discontinuity formulation
            u_star( internal::get_precurvatures( *this ) ),
            z_context( s_, y_, C_ ),
            g_context( y_, C_ )
         {
#ifdef COUNT_OPS
            g_nnorm += 1;
            g_ndiv += 5;
            g_nadd += kb.size();
#endif
            using namespace tuple_transform;
            //pretend that there is an "infinitely stiff" tube which reacts to all external
            //loads and also straightens the collection for all s < 0

#ifdef COUNT_OPS
            g_nadd += 1;
#endif
            if ((s.left + s.right) / 2.0 < 0) {
               kb_tot_inv = 0;
            }

            //get the precurvature list
            //copy(transform(C, y.Sigma, internal::evaluate_precurvature()),
            //   columns_of(u_star));

            //Get the sines and cosines
            //sincos = internal::compute_sines_and_cosines(*this);

            //get the F matrix
            F = internal::get_F_matrix( *this );

            //OK, just needs u_star and kb which are already initialized
            uB = internal::get_bishop_curvature( *this );

            //OK, needs ustar
            //M = internal::get_M_matrix(*this);

            z_context.construct( *this );
         }

         Mathematics::eval_pt const& s;
         State const& y;
         Cannula const& C;

         //Variables which appear commonly in equations
         typename Matrix<2, N::value>::type	sincos; //sine and cosine storage for angles
         typename Vector<N::value>::type		kb;
         double								      kb_tot_inv;
         typename Matrix<2, N::value>::type	u_star; //precurvature list
         typename Matrix<2, N::value>::type	F; //always needed for uB
         typename Vector<2>::type			   uB;//needed for most equations
         //typename Vector<4>::type            q_renorm; //renormalized version of quaternion

         //Z context
         Z_context                           z_context;
         G_context							 g_context;
      };
   }

   /**
      \brief The state equation for the concentric tube robot

      The state obeys a differential
      state equation, where if \f$y\f$ is a state, we have
      \f[
      \frac{ \partial }{\partial s}y = (L_y)_* \xi(s, y)
      \f]

      The entire right hand side of this equation is modeled by this function.
      This includes the evolution of all states along the robot, including the
      Jacobian states, internal load states, and geometric states if they are
      present for the problem being solved.

      Note that the evaluation here
      takes place for a Mathematics:eval_pt rather than at a single arc length.
      This is due to the presence of discontinuous functions representing the
      bending stiffness and torsional compliance of the individual tubes. In
      order to correctly evaluate this functions, the interval of integration is
      needed in addition to the value of the arc length. If the stiffnesses
      were assumed functions of arc length only, they would have to be modeled
      as either left continuous or right continuous, either of which would pose
      a problem for the integrator. To alleviate these problems, the stiffnesses
      are assumed to have multiple values at the single point of discontinuity.
      Note that under integration, the value of a function at a point does not
      matter (in Lebesgue integration theory) so making this choice helps the
      convergence of the integrator without changing the problem. It is a "hack"
      for numerical purposes.

      The function also takes a const reference to the cannula, which is needed
      to evaluate the precurvature and stiffness functions for each tube. In the
      CTR::Kinematics routine, this argument is bound to xi_fun before calling
      the integrator so that the integrator sees this function as only depending
      on the parameters s and y.

      \param s The evaluation point (arc length, left interval endpoint, right interval endpoint)
      \param y A state element
      \param C A cannula object (tuple of CTR::Tube objects of any types)
   */
   template <typename State, typename Cannula>
   State xi_fun( const Mathematics::eval_pt &s, const State &y, const Cannula &C )
   {
      //Declare the return value
      State yprime;

      //Create the evaluation context, which will evaluate commonly
      // used subexpressions first
      internal::evaluation_context<State, Cannula>	context( s, y, C );

      //compute the torsional state equations
      internal::xi_fun_impl_t( context, yprime );

      //Compute the geometric state equations (if applicable)
      internal::xi_fun_impl_g( context, yprime );

      //Compute the loading state equations (if applicable)
      internal::xi_fun_impl_l( context, yprime );

      //Compute the jacobian state equations (if applicable)
      internal::xi_fun_impl_Z( context, yprime );

	  //std::cout << "(s, uB) = (" << s << ", " << context.uB << ")\n";

      return yprime;
   }

}

#undef NTUBES
#undef NTUBESC
