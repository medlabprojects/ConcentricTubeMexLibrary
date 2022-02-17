/*******************************************************
*
* Copyright 2015 Vanderbilt University
* Author: Hunter B. Gilbert
*
*
********************************************************/

#pragma once

#include "Tube.h"
#include "SmoothStepFunction.h"
#include "IndicatorFunction.h"
#include "CTRTypedefs.h"
#include "SystemState.h"
#include "Cannula.h"
#include "Utility.h"
#include "TaggedInterval.h"
#include "FunctionIntervals.h"
#include "FunctionDiscontinuities.h"
#include "RungeKutta.h"
#include "Options.h"
#include "CodeTimer.h"
#include "TupleZip.h"
#include "Linspace.h"
#include "TupleTransform.h"
#include "TupleAccumulate.h"
#include "StateEquations.h"
#include <boost/container/static_vector.hpp>
#include <iostream>
#include <iterator>
#include <functional>
#include <algorithm>

//This is complete overkill, but better safe than sorry 
// here.  In some cases, such as for a processor with very little
// or no cache available, it could be wise to lower this number
// to limit the amount of required stack-space. The program will
// crash if not enough space is available.
#define MAX_INTEGRATION_POINTS	2048

namespace CTR
{
   enum IntegrationEventType
   {
      TUBE_DISCONTINUITY,
      REACTION_LOAD,
	  EVALUATE_AT_S
   };

   struct Discontinuity
   {
      Discontinuity( double l, int t, IntegrationEventType ev_ ) :
         loc( l ),
         tube_responsible( t ),
         ev( ev_ )
      {}
      double					loc;
      int						tube_responsible;
      IntegrationEventType	ev;
   };
   typedef boost::container::static_vector< Discontinuity, 50 >	DiscontinuityList;

   std::ostream& operator<<( std::ostream& str, const Discontinuity& disc )
   {
      str << disc.loc << ": " << disc.tube_responsible;
      return str;
   }
}

namespace CTR
{

   namespace internal
   {
      //Simple binder to turn the three-argument xi_fun function
      // into the required two-argument form for the numerical
      // IVP integrator
      template <class Cannula>
      struct	xi_bind
      {
         xi_bind( const Cannula &C ) : C_( C )
         {}

         template <class State>
         State operator()( const Mathematics::eval_pt &s, State const& y ) const
         {
            return CTR::xi_fun( s, y, C_ );
         }

      private:
         const Cannula &C_;
      };

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_t_sigma( State& y, const Configuration& q, const Cannula& C )
      {
         using namespace tuple_transform;
         typedef typename std::tuple_size<Cannula>::type N;
         typedef typename Vector<N::value>::type	VectorN;
         VectorN	L;
         copy( transform( C, internal::evaluate_length() ), L );
         VectorN S = q.Beta + L;
         double s_max = S.maxCoeff();

         y.Sigma = s_max*VectorN::Ones() - q.Beta;
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_t_loaded_impl( State &y, const Configuration &q, const Cannula &C, boost::mpl::true_ )
      {
         y.Mz( 0 ) = q.Ttip( 2 ); //add the tip torque!
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_t_loaded_impl( State &y, const Configuration &q, const Cannula &C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_t( State& y, const Configuration& q, const Cannula& C )
      {
         //Initialize psi and Mz states at the tip
         // of the robot.
         y.Psi = q.PsiL;
         y.Mz = Functions::zero_fun< typename Vector< internal::traits<State>::N::value >::type >()();

         //initialize Mz1 if a tip torque is present
         init_state_t_loaded_impl( y, q, C, typename internal::traits<State>::ExternalLoads() );

         //Initialize sigma states
         init_state_impl_t_sigma( y, q, C );
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_g_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         y.p << 0.0, 0.0, 0.0;
         y.q << 1.0, 0.0, 0.0, 0.0;
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_g_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {

      }

      //Internal load states always start at zero, even when 
      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_l_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         y.Mt = q.Ttip;
         y.Nt = q.Ftip;
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_l_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {

      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_g( State& y, const Configuration& q, const Cannula& C )
      {
         //Get the option type for the state
         typedef typename internal::traits<State>::Options	OType;
         init_state_impl_g_impl( y, q, C, typename Option::has_option<OType, Option::ComputeGeometry>::type() );
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_l( State& y, const Configuration& q, const Cannula& C )
      {
         typedef typename internal::traits<State>::Options	OType;
         init_state_impl_l_impl( y, q, C, typename Option::has_option<OType, Option::ExternalLoads>::type() );
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zta_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         const static int N = internal::traits<State>::N::value;
         y.Zta = Matrix<2 * N, 2 * N>::type::Zero();

         //Derivative block of psi with respect to alpha (psi(L) = alpha, so this is the identity matrix)
         y.Zta.template block<N, N>( 0, 0 ) = Matrix<N, N>::type::Identity();
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zta_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zga_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         const static int N = internal::traits<State>::N::value;
         y.Zga = Matrix<6, 2 * N>::type::Zero();
         //y.Zga( 2, N ) = -1.0;
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zga_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zgf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         y.Zgf = Matrix<6, 6>::type::Zero();
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zgf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zla_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         const static int N = internal::traits<State>::N::value;
         y.Zla = Matrix<6, 2 * N>::type::Zero();
         //This initial condition for Zla is based on the e3 cross n term in the derivative
         // of m with respect to s
         y.Zla( 3, N ) = -q.Ftip( 1 );
         y.Zla( 4, N ) = q.Ftip( 0 );

      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zla_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Ztf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         const static int N = internal::traits<State>::N::value;
         y.Ztf = Matrix<2 * N, 6>::type::Zero();
         y.Ztf( N, 5 ) = 1.0; //z-moment is directly carried by innermost tube
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Ztf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zlf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::true_ )
      {
         y.Zlf = Matrix<6, 6>::type::Identity();
      }

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Zlf_impl( State& y, const Configuration& q, const Cannula& C, boost::mpl::false_ )
      {}

      template <typename State, typename Configuration, typename Cannula>
      void init_state_impl_Z( State& y, const Configuration& q, const Cannula& C )
      {
         //TODO: Implement state initialization for
         init_state_impl_Zta_impl( y, q, C, typename internal::traits<State>::ComputeZta() );
         init_state_impl_Zga_impl( y, q, C, typename internal::traits<State>::ComputeZga() );
         init_state_impl_Zgf_impl( y, q, C, typename internal::traits<State>::ComputeZgf() );
         init_state_impl_Zla_impl( y, q, C, typename internal::traits<State>::ComputeZla() );
         init_state_impl_Ztf_impl( y, q, C, typename internal::traits<State>::ComputeZtf() );
         init_state_impl_Zlf_impl( y, q, C, typename internal::traits<State>::ComputeZlf() );
      }
   }

   template <typename Cannula, typename Configuration, typename Options>
   State< std::tuple_size<Cannula>::type::value, Options >
      init_state( const Cannula& C, const Configuration& q, Options )
   {
      typedef typename std::tuple_size<Cannula>::type	N;
      State<N::value, Options> y;
      internal::init_state_impl_t( y, q, C ); //Initialize the torsional states
      internal::init_state_impl_g( y, q, C ); //Initialize the geometry states
      internal::init_state_impl_l( y, q, C ); //Initialize the internal load states
      internal::init_state_impl_Z( y, q, C ); //Initialize the Jacobian matrices

      return y;
   }

   namespace internal
   {
      struct get_discontinuities
      {
         template <typename T>
         struct result;

         template <typename Tube, typename F>
         struct result < F( Tube ) >
         {
            typedef DiscontinuityList&	type;
         };

         template <typename Tube>
         DiscontinuityList& operator()( DiscontinuityList& list, std::tuple< const Tube&, const int > data ) const
         {
            using namespace boost::container;
            const Tube& t = std::get<0>( data );
            int Number = std::get<1>( data );

            static_vector<double, 50>	temp_list;
            Functions::add_discontinuities( t.GetMaterialIndicator(), temp_list );

            struct make_discontinuity
            {
               Discontinuity operator()( double d, int t )
               {
                  return Discontinuity( d, t, TUBE_DISCONTINUITY );
               }
            };

            transform( temp_list.begin(), temp_list.end(), std::back_inserter( list ),
                       std::bind( make_discontinuity(), std::placeholders::_1, Number ) );

            return list;
         }
      };

      template <typename VectorType>
      void shift_discontinuities( DiscontinuityList& list, VectorType const &beta )
      {
         std::for_each( list.begin(), list.end(),
                        [&beta]( Discontinuity &v ) {
            v.loc += beta( v.tube_responsible );
         } );
      }

      struct tubes_to_intervals
      {
         template <typename T>
         struct result;

         template <typename Tube, typename F>
         struct result < F( Tube ) >
         {
            typedef TInterval::IntervalList&	type;
         };

         struct shift_interval
         {
            shift_interval( double d ) : shift( d )
            {}

            TInterval::TaggedInterval operator()( TInterval::TaggedInterval const& interval ) const
            {
               return TInterval::TaggedInterval( TInterval::left( interval ) + shift,
                                                 TInterval::right( interval ) + shift,
                                                 TInterval::density( interval ) );
            }
         private:
            double shift;
         };

         template <typename Tube>
         TInterval::IntervalList&	operator()( TInterval::IntervalList& il, std::tuple<const Tube&, const double&> data ) const
         {
            TInterval::IntervalList temp;
            const Tube& t = std::get<0>( data );
            add_intervals( t.GetCurvatureFun(), temp ); //intervals in curvature function
            add_intervals( t.GetKbFun(), temp );        //intervals in stiffness function
            add_intervals( t.GetCtFun(), temp );        //intervals in torsional compliance function
            add_intervals( t.GetMaterialIndicator(), temp );  //interval for tube itsel
            add_intervals( t.GetCurvatureIndicator(), temp ); //intervals for smooth steps on curvature function
            shift_interval si( std::get<1>( data ) );   //the interval shifter (models actuation in beta variables)
            std::transform( temp.begin(), temp.end(), back_inserter( il ), si );
            //Add the transmission interval (note that if this has zero length, it will be ignored by 
            // the resolve_to_disjoint function).
            il.push_back( si( TInterval::TaggedInterval( 0.0, t.GetTransmissionLength(), TInterval::SPARSE ) ) );
            //Add the regular interval
            il.push_back( si( TInterval::TaggedInterval( t.GetTransmissionLength(), t.GetLength(), TInterval::NORMAL ) ) );
            return il;
         }
      };
   }

   namespace internal
   {
      //Implementation of return structure which provides
      // tip location and orientation for convenience
      struct KinRet_impl_g
      {
         Vector<3>::type	pTip;
         Vector<4>::type	qTip;
      };

      //Empty struct in case no geometry is requested
      struct KinRet_impl_no_g
      {};

      //The return implementation for computing the tip position
      // and orientation
      template <typename Ret, typename Configuration>
      void return_impl_g( Ret &r, Configuration const& q, boost::mpl::true_ )
      {
         using namespace Utility;
         //compute the tip position
         r.pTip = (-quaternion_rotate( quaternion_conj( r.y_final.q ), r.y_final.p ) + q.Beta( 0 )*Eigen::Vector3d::UnitZ());
         //compute the tip orientation
         r.qTip = quaternion_conj( r.y_final.q );
      }

      template <typename Ret, typename Configuration>
      void return_impl_g( Ret &r, Configuration const& q, boost::mpl::false_ )
      {

      }
   }

   template <typename State>
   struct KinRet :
      public boost::mpl::if_< typename internal::traits<State>::ComputeGeometry, internal::KinRet_impl_g, internal::KinRet_impl_no_g >::type
   {
      State y_final;
      int normal_points_used;
   };

   template <typename State>
   struct KinRetDense :
      public boost::mpl::if_< typename internal::traits<State>::ComputeGeometry, internal::KinRet_impl_g, internal::KinRet_impl_no_g >::type
   {
      State y_final;
      std::vector<State, Eigen::aligned_allocator<State> > dense_state_output;
      std::vector<double> arc_length_points;
      int normal_points_used;

   };

   namespace internal
   {
      template <typename T>
      struct traits;

      template <typename State_>
      struct traits < KinRet<State_> >
      {
         typedef State_ State;
         typedef typename traits<State>::N  N;
      };
   }

   namespace internal
   {
      bool discontinuity_compare( Discontinuity const& disc1, Discontinuity const& disc2 )
      {
         return (disc1.loc < disc2.loc);
      }

      template <typename Cannula, typename Configuration>
      DiscontinuityList	get_discontinuity_list( Cannula const& C, Configuration const &q )
      {
         DiscontinuityList discs;
         //zip_with_index passes each element along with its index in the tuple
         tuple_transform::accumulate_in_place( tuple_transform::zip_with_index( C ), internal::get_discontinuities(), discs );
         internal::shift_discontinuities( discs, q.Beta );

         ////push back the "zero arc length" discontinuity where reaction loads take place
         discs.push_back( Discontinuity( 0.0, -1, REACTION_LOAD ) );

         std::sort( discs.begin(), discs.end(), discontinuity_compare );

         return discs;
      }

      template <typename Cannula, typename Configuration >
      void add_loaded_interval( TInterval::IntervalList &ilist, Cannula const &C, Configuration const &q, boost::mpl::true_ )
      {
         double end = q.Beta( 0 ) + std::get<0>( C ).GetLength();
         ilist.push_back( TInterval::TaggedInterval( 0.0, end, TInterval::NORMAL ) );
      }

      template <typename Cannula, typename Configuration >
      void add_loaded_interval( TInterval::IntervalList &ilist, Cannula const &C, Configuration const &q, boost::mpl::false_ )
      {}


      template <typename Cannula, typename Configuration, typename OType>
      TInterval::IntervalList	get_integration_intervals( Cannula const& C, Configuration const &q, DiscontinuityList const& discs, OType )
      {
         using namespace TInterval;
         IntervalList intervals;
         tuple_transform::accumulate_in_place( tuple_transform::zip( C, q.Beta ), tubes_to_intervals(), intervals );
         //Resolve the integration intervals to a union of adjacent intervals
         std::transform( discs.begin(), discs.end(), back_inserter( intervals ),
         []( Discontinuity const &d ) {
            return TaggedInterval( d.loc, d.loc, DISCONTINUITY );
         } );

         add_loaded_interval( intervals, C, q, CTR::Option::has_option<OType, Option::ExternalLoads>() );

         intervals = resolve_to_disjoint( intervals );

         return intervals;
      }

      double get_normal_length( TInterval::IntervalList const& ilist )
      {
         using namespace TInterval;
         double L = 0;
         for (auto it = ilist.begin(); it != ilist.end(); ++it) {
            if (has_tag( NORMAL )(*it)) {
               L += TInterval::length( *it );
            }
         }
         return L;
      }
   }

   namespace internal
   {
      typedef boost::container::static_vector<double, MAX_INTEGRATION_POINTS>	IntegrationPoints;
      void make_normal_points( TInterval::TaggedInterval const& ival, double normal_density, IntegrationPoints& pts )
      {
         using namespace TInterval;
         double L = length( ival );

         //2+ is for the endpoints of the interval
         int NPts = 2 + (int)floor( L * normal_density );
         NPts = (NPts > 2) ? NPts : 2;

         Utility::linspace<double> space( left( ival ), right( ival ), NPts );
         auto it_end = space.end(); it_end-=1; //don't add the last point
         std::copy( space.begin(), it_end, std::back_inserter( pts ) );
      }

      void make_dense_points( TInterval::TaggedInterval const& ival, double dense_density, IntegrationPoints& pts )
      {
         using namespace TInterval;
         double L = length( ival );

         //2+ is for the endpoints of the interval
         int NPts = 2 + (int)ceil( L * dense_density );
         NPts = (NPts > 2) ? NPts : 2;

         Utility::linspace<double> space( left( ival ), right( ival ), NPts );
         auto it_end = space.end(); it_end-=1; //don't add the last point
         std::copy( space.begin(), it_end, std::back_inserter( pts ) );
      }

      void make_sparse_points( TInterval::TaggedInterval const& ival, IntegrationPoints& pts )
      {
         pts.push_back( TInterval::left( ival ) );
      }

   }

   namespace internal
   {
      void make_integration_points( TInterval::TaggedInterval const& ival, double normal_density, double dense_density,
                                    internal::IntegrationPoints &pts )
      {
         using namespace boost::container;

         switch (TInterval::density( ival )) {
         case TInterval::NORMAL:
            internal::make_normal_points( ival, normal_density, pts );
            break;
         case TInterval::DENSE:
            internal::make_dense_points( ival, dense_density, pts );
            break;
         case TInterval::SPARSE:
            internal::make_sparse_points( ival, pts );
            break;
         }

         pts.push_back( right( ival ) );
      }
   }

   namespace internal
   {
      template <typename State, typename Context>
      void transition_impl_Zga( State &y, Context const& context_before, Context const& context_after, Discontinuity const &disc, boost::mpl::true_ )
      {
         const int N = traits<State>::N::value;

         typename Vector<2>::type term_before = context_before.uB;
         typename Vector<2>::type term_after = context_after.uB;

         y.Zga.col( N + disc.tube_responsible ).template middleRows<2>( 3 ) += (term_before - term_after);
      }

      template <typename State, typename Context>
      void transition_impl_Zga( State &y, Context const& eval_before, Context const& eval_after, Discontinuity const &disc, boost::mpl::false_ )
      {}

      template <typename State, typename Context>
      void transition_impl_Zta( State &y, Context const& context_before, Context const& context_after, Discontinuity const &disc, boost::mpl::true_ )
      {
         const static int N = traits<State>::N::value;

         auto g_before = get_moment_derivatives( context_before );
         auto g_after = get_moment_derivatives( context_after );

         auto xi_psi_before = get_angle_derivatives( context_before );
         auto xi_psi_after = get_angle_derivatives( context_after );

         y.Zta.col( N + disc.tube_responsible ).template topRows<N>() += (xi_psi_before - xi_psi_after);
         y.Zta.col( N + disc.tube_responsible ).template bottomRows<N>() += (g_before - g_after);
      }

      template <typename State, typename Context>
      void transition_impl_Zta( State &y, Context const& eval_before, Context const& eval_after, Discontinuity const &disc, boost::mpl::false_ )
      {

      }

      template <typename State>
      void transition_impl_l( State &y, Mathematics::eval_pt const& eval_before, Mathematics::eval_pt const& eval_after, Discontinuity const& d, boost::mpl::false_ )
      {

      }

      template <typename State>
      void transition_impl_l( State &y, Mathematics::eval_pt const& eval_before, Mathematics::eval_pt const& eval_after, Discontinuity const& d, boost::mpl::true_ )
      {
         if (d.ev == REACTION_LOAD) {
            y.Nt = Vector<3>::type::Zero();
            y.Mt = Vector<3>::type::Zero();
         }
      }

      template <typename State, typename Context>
      void transition_impl_Zla( State &y, Context const& eval_before, Context const& eval_after, Discontinuity const& d, boost::mpl::false_ )
      {

      }

      template <typename State, typename Context>
      void transition_impl_Zla( State &y, Context const& context_before, Context const& context_after, Discontinuity const& d, boost::mpl::true_ )
      {
         const int N = traits<State>::N::value;
         int tube_id = d.tube_responsible;

         typename Vector<3>::type   xi_m_before = internal::xi_m( context_before );
         typename Vector<3>::type   xi_m_after = internal::xi_m( context_after );

         typename Vector<3>::type   xi_n_before = internal::xi_n( context_before );
         typename Vector<3>::type   xi_n_after = internal::xi_n( context_after );

         y.Zla.col( N + tube_id ).template topRows<3>() += xi_n_before - xi_n_after;
         y.Zla.col( N + tube_id ).template bottomRows<3>() += xi_m_before - xi_m_after;
      }

      template <typename State, typename Context>
      void transition_impl_Zlf( State &y, Context const& context_before, Context const& context_after, Discontinuity const& d, boost::mpl::false_ )
      {

      }

      template <typename State, typename Context>
      void transition_impl_Zlf( State &y, Context const& context_before, Context const& context_after, Discontinuity const& d, boost::mpl::true_ )
      {
         typedef typename Matrix<6, 6>::type	MType;
         y.Zlf = Functions::zero_fun<MType>()();
      }

      template <typename State>
      void transition_impl_y( State &y, Mathematics::eval_pt const &eval_before, Mathematics::eval_pt const &eval_after, Discontinuity const &d )
      {
         transition_impl_l( y, eval_before, eval_after, d, typename traits<State>::ExternalLoads() );
      }

      template <typename State, typename Cannula>
      void transition_impl( State &y, Cannula const &C, Mathematics::eval_pt const& eval_before, Mathematics::eval_pt const& eval_after, Discontinuity const &disc )
      {
         ////transition state variables if Mb, Nb are discontinuous here (and if they exist)
         //transition_impl_y( y_after, eval_before, eval_after, disc ); 

         //Get the evaluation contexts for the interval and state on each side of the discontinuity
         internal::evaluation_context<State, Cannula>  context_before( eval_before, y, C );
         internal::evaluation_context<State, Cannula>  context_after( eval_after, y, C );

         const int N = traits<State>::N::value;

         //Do the transitions on Zga, Zta, and Zla if due to a tube discontinuity
         if (disc.ev == TUBE_DISCONTINUITY) {
            transition_impl_Zla( y, context_before, context_after, disc, typename traits<State>::ComputeZla() );
            transition_impl_Zga( y, context_before, context_after, disc, typename traits<State>::ComputeZga() );
            transition_impl_Zta( y, context_before, context_after, disc, typename traits<State>::ComputeZta() );
         }
      }

      template <typename State, typename Cannula>
      void inner_transition( State& y, Cannula const& C, TInterval::TaggedInterval const& before, TInterval::TaggedInterval const& after, Discontinuity const& disc )
      {
         double s = TInterval::left( before );
         Mathematics::eval_pt eval_before( s, left( before ), right( before ) );
         Mathematics::eval_pt eval_after( s, left( after ), right( after ) );

         transition_impl( y, C, eval_before, eval_after, disc );
      }

      struct state_callback
      {
         template <typename State>
         void callback_impl( State &y, boost::mpl::false_ ) const
         {}

         template <typename State>
         void callback_impl( State &y, boost::mpl::true_ ) const
         {
#ifdef COUNT_OPS
            g_nmult += 4;
            g_ndiv += 1;
            g_nnorm += 1;
#endif
            y.q = (1.0 / y.q.norm()) * y.q;
         }

         template <typename State>
         void operator()( State &y ) const
         {
            typedef typename internal::traits<State>::ComputeGeometry   GeomOpt;
            callback_impl( y, GeomOpt() );
         }
      };

      template <typename State, typename Cannula>
      State integrate_impl( Cannula const& C, State const& init_state, TInterval::IntervalList const& ivals, DiscontinuityList const& discs,
                            int n_normal_points, double dense_density, int *normal_points_used = 0 )
      {
         using namespace boost::container;
         double normal_density = static_cast<double>(n_normal_points) / get_normal_length( ivals );
         State y_final = init_state;

         int normal_pts_used = 0;

         //TODO: Fix this. The comment below is incorrect - the last discontinuity must be integrated over.
         //Note: first and last intervals are always discontinuities, so no need to 
         // integrate over these
         for (auto it = ivals.rbegin() + 1; it != ivals.rend() - 1; ++it) {
            switch (density( *it )) {
            case TInterval::DENSE:
            {
               if (TInterval::length( *it ) > 10 * std::numeric_limits<double>::epsilon()) {
                  IntegrationPoints	pts;
                  make_integration_points( *it, normal_density, dense_density, pts );
                  y_final = Mathematics::rk8_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend(), state_callback() );
               }
               break;
            }
            case TInterval::SPARSE:
            {
               IntegrationPoints pts;
               make_integration_points( *it, normal_density, dense_density, pts );

               y_final = Mathematics::euler_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend() );
               break;
            }
            case TInterval::NORMAL:
            {
               if (TInterval::length( *it ) > 10 * std::numeric_limits<double>::epsilon()) {
                  IntegrationPoints pts;
                  make_integration_points( *it, normal_density, dense_density, pts );
                  normal_pts_used += static_cast<int>(pts.size()) - 2;

                  y_final = Mathematics::rk8_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend(), state_callback() );
               }
               break;
            }
            case TInterval::DISCONTINUITY:
            {
               double loc = TInterval::left( *it );
               //This really should have a better implementation, but it does work and the call to 
               // find_if is not a speed penalty (has been profiled)
               auto disc_it = std::find_if( discs.begin(), discs.end(), [loc]( Discontinuity const& d ) {
                  return (d.loc == loc);
               } );

               //transition across the discontinuity
               internal::inner_transition( y_final, C, *(it - 1), *(it + 1), *disc_it );
               break;
            }
            }
         }
         if (normal_points_used != 0)
            *normal_points_used = normal_pts_used;
         return y_final;
      }

      template <typename State, typename Cannula, typename OutputContainer, typename OutputArcLengths>
      State integrate_impl_dense( Cannula const& C, State const& init_state, TInterval::IntervalList const& ivals, DiscontinuityList const& discs,
                                  int n_normal_points, double dense_density, OutputContainer &out, OutputArcLengths &outS, int *normal_points_used = 0 )
      {
         using namespace boost::container;
         double normal_density = static_cast<double>(n_normal_points) / get_normal_length( ivals );
         State y_final = init_state;

         int normal_pts_used = 0;

         //Note: first and last intervals are always discontinuities, so no need to 
         // integrate over these
         for (auto it = ivals.rbegin() + 1; it != ivals.rend() - 1; ++it) {
            switch (density( *it )) {
            case TInterval::DENSE:
            {
               if (TInterval::length( *it ) > 10 * std::numeric_limits<double>::epsilon()) {
                  IntegrationPoints	pts;
                  make_integration_points( *it, normal_density, dense_density, pts );

                  //y_final = Mathematics::rk8_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend(), state_callback() );
                  Mathematics::integrate( pts.rbegin(), pts.rend(), std::back_inserter( out ), y_final, Mathematics::rk8_step( xi_bind<Cannula>( C ) ), state_callback() );
                  y_final = out.back();
                  std::copy( pts.rbegin(), pts.rend(), std::back_inserter( outS ) );

                  outS.pop_back();
                  out.pop_back();
               }
               break;
            }
            case TInterval::SPARSE:
            {
               IntegrationPoints pts;
               make_integration_points( *it, normal_density, dense_density, pts );

               //y_final = Mathematics::euler_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend() );
               Mathematics::integrate( pts.rbegin(), pts.rend(), std::back_inserter( out ), y_final, Mathematics::euler_step( xi_bind<Cannula>( C ) ) );
               y_final = out.back();
               std::copy( pts.rbegin(), pts.rend(), std::back_inserter( outS ) );

               outS.pop_back();
               out.pop_back();
               break;
            }
            case TInterval::NORMAL:
            {
               if (TInterval::length( *it ) > 10 * std::numeric_limits<double>::epsilon()) {
                  IntegrationPoints pts;
                  make_integration_points( *it, normal_density, dense_density, pts );
                  normal_pts_used += static_cast<int>(pts.size()) - 2;

                  //y_final = Mathematics::rk8_final( xi_bind<Cannula>( C ), y_final, pts.rbegin(), pts.rend(), state_callback() );
                  Mathematics::integrate( pts.rbegin(), pts.rend(), std::back_inserter( out ), y_final, Mathematics::rk4_step( xi_bind<Cannula>( C ) ), state_callback() );
                  y_final = out.back();
                  std::copy( pts.rbegin(), pts.rend(), std::back_inserter( outS ) );

                  outS.pop_back();
                  out.pop_back();
               }
               break;
            }
            case TInterval::DISCONTINUITY:
            {
               double loc = TInterval::left( *it );
               //This really should have a better implementation, but it does work and the call to 
               // find_if is not really a speed penalty (has been profiled)
               auto disc_it = std::find_if( discs.begin(), discs.end(), [loc]( Discontinuity const& d ) {
                  return (d.loc == loc);
               } );

               //transition across the discontinuity
               internal::inner_transition( y_final, C, *(it - 1), *(it + 1), *disc_it );
               break;
            }
            }
         }
         if (normal_points_used != 0)
            *normal_points_used = normal_pts_used;

         outS.push_back( left( ivals.front() ) );
         out.push_back( y_final );

         return y_final;
      }
   }

   //static LARGE_INTEGER t1, t2;
   //static LARGE_INTEGER total_integration_time;
   //static LARGE_INTEGER total_setup_time;

   /**
      \brief Concentric Tube Robot Kinematics

      \param C The concentric tube robot, or "Cannula", as a std::tuple
      \param q The configuration of the robot, as a struct. This struct has the following format
      ~~~{.cpp}
      struct Configuration
      {
      Eigen::Matrix<double,N,1>	Psi;
      Eigen::Matrix<double,N,1>	Beta;
      };
      ~~~
      \param Options An options object, which is discussed in the documentation for DeclarOptions
      \param n_normal_points The number of additional integrator points to add in normal intervals
      \param dense_density The minimum allowable density of points in dense intervals

      \return a KinRet<State>	object with the result of the forward kinematics
      */
   template < typename Cannula, typename Configuration, typename Options >
   KinRet< State< std::tuple_size<Cannula>::type::value, Options> >
      Kinematics( const Cannula& C, const Configuration& q, const Options o,
                  int n_normal_points = 0, double dense_density = 0.0 )
   {
      using namespace TInterval;
      using namespace boost::container;
      typedef typename std::tuple_size<Cannula>::type	 N;
      typedef typename Vector<N::value>::type	VectorN;

      //Get the list of discontinuities
      //QueryPerformanceCounter( &t1 );
      DiscontinuityList discs = internal::get_discontinuity_list( C, q );

      //Get the integration intervals
      TInterval::IntervalList	intervals = internal::get_integration_intervals( C, q, discs, Options() );

      //get the initial state
      State<N::value, Options> y_init = init_state( C, q, Options() );

      //Perform the (backwards) integration of the IVP
      KinRet< State<N::value, Options> >	ret;
      //QueryPerformanceCounter( &t2 );
      //total_setup_time.QuadPart += t2.QuadPart - t1.QuadPart;

      //QueryPerformanceCounter( &t1 );
      ret.y_final = internal::integrate_impl( C, y_init, intervals, discs, n_normal_points, dense_density, &ret.normal_points_used );

      //If the geometry is specified, compute the tip location and orientation with respect to the
      // base frame at p(beta1) = [0 0 beta1], q(beta1) = [1 0 0 0]
      internal::return_impl_g( ret, q, Option::has_option<Options, Option::ComputeGeometry>() );
      //QueryPerformanceCounter( &t2 );
      //total_integration_time.QuadPart += t2.QuadPart - t1.QuadPart;

      return ret;
   }

   template < typename Cannula, typename Configuration, typename Options >
   KinRetDense< State< std::tuple_size<Cannula>::type::value, Options> >
      Kinematics_with_dense_output( const Cannula& C, const Configuration& q, const Options o,
                                    int n_normal_points = 0, double dense_density = 0.0, double sEval = 0.0 )
   {
      using namespace TInterval;
      using namespace boost::container;
      typedef typename std::tuple_size<Cannula>::type	 N;
      typedef typename Vector<N::value>::type	VectorN;

      //Get the list of discontinuities
      DiscontinuityList discs = internal::get_discontinuity_list( C, q );

	  /*for(int i = 0; i < 10; i++)
	  std::cout << discs[i] << '\n';

	  std::cout << '\n' << "break" << '\n' << '\n';*/

	  discs.push_back(Discontinuity(sEval, -1, EVALUATE_AT_S));

	  /*for (int i = 0; i < 10; i++)
	  std::cout << discs[i] << '\n';*/

      //Get the integration intervals
      TInterval::IntervalList	intervals = internal::get_integration_intervals( C, q, discs, Options() );

      //get the initial state
      State<N::value, Options> y_init = init_state( C, q, Options() );

      //Perform the (backwards) integration of the IVP
      KinRetDense< State<N::value, Options> >	ret;

      ret.y_final = internal::integrate_impl_dense( C, y_init, intervals, discs, n_normal_points, dense_density, ret.dense_state_output, ret.arc_length_points, &ret.normal_points_used );

      //If the geometry is specified, compute the tip location and orientation with respect to the
      // base frame at p(beta1) = [0 0 beta1], q(beta1) = [1 0 0 0]
      internal::return_impl_g( ret, q, Option::has_option<Options, Option::ComputeGeometry>() );

      return ret;
   }

   /**
      \brief Get the Jacobian matrix relating tip motion to q

      \param y The final state vector of a kinematics run which has done the Jacobian matrix computations

      \return The Body-frame Jacobian matrix at the tip of the robot
      */
   template <typename State>
   typename Matrix< 6, 2 * internal::traits<State>::N::value >::type
      GetTipJacobianForTube1( const State &y )
   {
      const int N = internal::traits<State>::N::value;
      Eigen::Matrix<double, 6, 6> Ad = Utility::Adjoint_p_q( y.p, y.q );
      Eigen::Matrix<double, 6, 2 * N > J;
      J = -Ad*y.Zga;
      J( 2, N ) += 1; //arc-length attach
      return J;
   }

   /**
   \brief Get the Jacobian matrix relating tip motion to the external forces

   \param y The final state vector of a kinematics run which has done the Jacobian matrix computations

   \return The Body-frame Compliance matrix at the tip of the robot
   */
   template <typename State>
   typename Matrix< 6, 2 * internal::traits<State>::N::value >::type
      GetTipComplianceForTube1( const State &y )
   {
      Eigen::Matrix<double, 6, 6> Ad = Utility::Adjoint_p_q( y.p, y.q );

      return -Ad*y.Zgf;
   }

   /**
      \brief Get the solution stability from the final state

      Returns the solution stability, as computed by
      \f[ \text{stability} = \text{det}(Z_{\psi,\alpha}(\beta_1)) \f]
      */
   template <typename State>
   double GetStability( const State &y )
   {
      const int N = internal::traits<State>::N::value;
      return ((y.Zta).template block<N, N>( 0, 0 )).determinant();
   }

   namespace internal
   {
      template <typename Configuration>
      struct offset_psi
      {
         template <typename VType>
         Configuration	operator()( const Configuration &q, const VType &dir ) const
         {
            Configuration qnew = q;
            qnew.PsiL += dir;
            return qnew;
         }
      };

      template <typename Configuration>
      struct offset_beta
      {
         template <typename VType>
         Configuration	operator()( const Configuration &q, const VType &dir ) const
         {
            Configuration qnew = q;
            qnew.Beta += dir;
            return qnew;
         }
      };

      template <typename Cannula, typename Configuration, typename KinRet, typename OffsetFun>
      struct finite_differencer_obj
      {
         finite_differencer_obj( const Cannula &C_, const Configuration &q_, const KinRet &y_init, const OffsetFun &f_ ) :
            C( C_ ),
            q( q_ ),
            y_mid( y_init ),
            f( f_ )
         {}

         template <typename F> struct result;

         template <typename VType>
         struct result < finite_differencer_obj( const VType& ) >
         {
            typedef typename Vector<6>::type	type;
         };

         template <typename VType>
         typename Vector<6>::type	operator()( const VType &dir ) const
         {
            Configuration qpos = f( q, dir );
            Configuration qneg = f( q, -dir );

            typedef typename CTR::DeclareOptions< Option::ComputeGeometry >::options	OType;
            auto kinret_pos = Kinematics( C, qpos, OType(), 50, 2.0 / 0.001 );
            auto kinret_neg = Kinematics( C, qneg, OType(), 50, 2.0 / 0.001 );

            Eigen::Vector4d dq = Utility::quaternion_mult( Utility::quaternion_conj( y_mid.qTip ), kinret_pos.qTip - kinret_neg.qTip );
            Eigen::Vector3d dp = kinret_pos.pTip - kinret_neg.pTip;

            typename Vector<6>::type	v;
            v.topRows<3>() = Utility::quaternion_rotate( Utility::quaternion_conj( y_mid.qTip ), dp );
            v.bottomRows<3>() = 2.0*dq.bottomRows<3>();

            return v;
         }

      private:
         const KinRet &y_mid;
         const Configuration &q;
         const Cannula &C;
         OffsetFun	  f;
      };

      template <typename Cannula, typename Configuration, typename KinRet, typename OffsetFun>
      struct loaded_finite_differencer_obj
      {
         loaded_finite_differencer_obj( const Cannula &C_, const Configuration &q_, const KinRet &y_init, const OffsetFun &f_ ) :
            C( C_ ),
            q( q_ ),
            y_mid( y_init ),
            f( f_ )
         {}

         template <typename F> struct result;

         template <typename VType>
         struct result < loaded_finite_differencer_obj( const VType& ) >
         {
            typedef typename Vector<6>::type	type;
         };

         template <typename VType>
         typename Vector<6>::type	operator()( const VType &dir ) const
         {
            Configuration qpos = f( q, dir );
            Configuration qneg = f( q, -dir );

            typedef typename CTR::DeclareOptions< Option::ComputeGeometry, Option::ExternalLoads >::options	OType;
            auto kinret_pos = Kinematics( C, qpos, OType(), 300, 10. / 0.001 );
            auto kinret_neg = Kinematics( C, qneg, OType(), 300, 10. / 0.001 );

            Eigen::Vector4d dq = Utility::quaternion_mult( Utility::quaternion_conj( y_mid.qTip ), kinret_pos.qTip - kinret_neg.qTip );
            Eigen::Vector3d dp = kinret_pos.pTip - kinret_neg.pTip;

            typename Vector<6>::type	v;
            v.topRows<3>() = Utility::quaternion_rotate( Utility::quaternion_conj( y_mid.qTip ), dp );
            v.bottomRows<3>() = 2.0*dq.bottomRows<3>();

            return v;
         }

      private:
         const KinRet &y_mid;
         const Configuration &q;
         const Cannula &C;
         OffsetFun	  f;
      };

      template <typename Cannula, typename Configuration, typename KinRet>
      finite_differencer_obj<Cannula, Configuration, KinRet, offset_psi<Configuration> >
         fd_psi( const Cannula &C, const Configuration &q, const KinRet &y_mid )
      {
         return finite_differencer_obj<Cannula, Configuration, KinRet, offset_psi<Configuration> >(
            C,
            q,
            y_mid,
            offset_psi<Configuration>()
            );
      }

      template <typename Cannula, typename Configuration, typename KinRet>
      loaded_finite_differencer_obj<Cannula, Configuration, KinRet, offset_psi<Configuration> >
         fd_psi_loaded( const Cannula &C, const Configuration &q, const KinRet &y_mid )
      {
         return loaded_finite_differencer_obj<Cannula, Configuration, KinRet, offset_psi<Configuration> >(
            C,
            q,
            y_mid,
            offset_psi<Configuration>()
            );
      }

      template <typename Cannula, typename Configuration, typename KinRet>
      finite_differencer_obj<Cannula, Configuration, KinRet, offset_beta<Configuration> >
         fd_beta( const Cannula &C, const Configuration &q, const KinRet &y_mid )
      {
         return finite_differencer_obj<Cannula, Configuration, KinRet, offset_beta<Configuration> >(
            C,
            q,
            y_mid,
            offset_beta<Configuration>()
            );
      }

      template <typename Cannula, typename Configuration, typename KinRet>
      loaded_finite_differencer_obj<Cannula, Configuration, KinRet, offset_beta<Configuration> >
         fd_beta_loaded( const Cannula &C, const Configuration &q, const KinRet &y_mid )
      {
         return loaded_finite_differencer_obj<Cannula, Configuration, KinRet, offset_beta<Configuration> >(
            C,
            q,
            y_mid,
            offset_beta<Configuration>()
            );
      }

   }

   template <typename Cannula, typename Configuration>
   typename Matrix< 6, 2 * std::tuple_size<Cannula>::value >::type
      GetFDJacobianForTube1( const Cannula &C, const Configuration &q, const double step_size_psi = 1e-5, const double step_size_beta = 1e-5 )
   {
      const int NTUBES = std::tuple_size<Cannula>::value;

      typedef typename CTR::DeclareOptions< Option::ComputeGeometry >::options	OType;
      typedef State<NTUBES, OType >	SType;

      typename Matrix<6, NTUBES >::type	JbPsi;
      typename Matrix<6, NTUBES >::type	JbBeta;
      typename Matrix<NTUBES, NTUBES>::type	DirsPsi = step_size_psi*Matrix<NTUBES, NTUBES>::type::Identity();
      typename Matrix<NTUBES, NTUBES>::type	DirsBeta = step_size_beta*Matrix<NTUBES, NTUBES>::type::Identity();

      auto y_init = Kinematics( C, q, OType(), 50, 2. / 0.001 );

      auto cols_psi = tuple_transform::transform( tuple_transform::columns_of( DirsPsi ), internal::fd_psi( C, q, y_init ) );
      auto cols_beta = tuple_transform::transform( tuple_transform::columns_of( DirsBeta ), internal::fd_beta( C, q, y_init ) );

      tuple_transform::copy( cols_psi, tuple_transform::columns_of( JbPsi ) );
      tuple_transform::copy( cols_beta, tuple_transform::columns_of( JbBeta ) );

      typename Matrix<6, 2 * NTUBES >::type	Jb;
      Jb.template block<6, NTUBES >( 0, 0 ) = (1.0 / 2.0 / step_size_psi)*JbPsi;
      Jb.template block<6, NTUBES >( 0, NTUBES ) = (1.0 / 2.0 / step_size_beta)*JbBeta;

      return Jb;
   }

   template <typename Cannula, typename Configuration>
   typename Matrix< 6, 2 * std::tuple_size<Cannula>::value >::type
      GetLoadedFDJacobianForTube1( const Cannula &C, const Configuration &q, const double step_size_psi = 1e-5, const double step_size_beta = 1e-5 )
   {
      const int NTUBES = std::tuple_size<Cannula>::value;

      typedef typename CTR::DeclareOptions< Option::ComputeGeometry, Option::ExternalLoads >::options	OType;
      typedef State<NTUBES, OType >	SType;

      typename Matrix<6, NTUBES >::type	JbPsi;
      typename Matrix<6, NTUBES >::type	JbBeta;
      typename Matrix<NTUBES, NTUBES>::type	DirsPsi = step_size_psi*Matrix<NTUBES, NTUBES>::type::Identity();
      typename Matrix<NTUBES, NTUBES>::type	DirsBeta = step_size_beta*Matrix<NTUBES, NTUBES>::type::Identity();

      auto y_init = Kinematics( C, q, OType(), 300, 10. / 0.001 );

      auto cols_psi = tuple_transform::transform( tuple_transform::columns_of( DirsPsi ), internal::fd_psi_loaded( C, q, y_init ) );
      auto cols_beta = tuple_transform::transform( tuple_transform::columns_of( DirsBeta ), internal::fd_beta_loaded( C, q, y_init ) );

      tuple_transform::copy( cols_psi, tuple_transform::columns_of( JbPsi ) );
      tuple_transform::copy( cols_beta, tuple_transform::columns_of( JbBeta ) );

      typename Matrix<6, 2 * NTUBES >::type	Jb;
      Jb.template block<6, NTUBES >( 0, 0 ) = (1.0 / 2.0 / step_size_psi)*JbPsi;
      Jb.template block<6, NTUBES >( 0, NTUBES ) = (1.0 / 2.0 / step_size_beta)*JbBeta;

      return Jb;
   }

}
