unsigned long long g_ntrig = 0;
unsigned long long g_nadd = 0;
unsigned long long g_nmult = 0;
unsigned long long g_nsqrt = 0;
unsigned long long g_ndiv = 0;
unsigned long long g_nnorm = 0;

#include <iostream>
#include <Eigen/Dense>
#include "Tube.h"
#include "SmoothStepFunction.h"
#include "IndicatorFunction.h"
#include "CTRTypedefs.h"
#include "SystemState.h"
#include "Cannula.h"
#include "Utility.h"
#include "TaggedInterval.h"
#include "Kinematics.h"
#include <tuple>
#include "RungeKutta.h"
#include "CodeTimer.h"
#include <memory>
#include <algorithm>
#include <vector>
#include <fstream>
#include "TupleTransform.h"
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <iomanip>
#include "Options.h"
#include "FunctionIntervals.h"
#include "BumpFunction.h"
#include <boost/io/ios_state.hpp>
#include <random>
#include <math.h>

#ifdef WIN32
#include "Processthreadsapi.h"
#endif

#undef max
#undef min

using namespace std;
using namespace CTR;

struct test_derivative
{
   double operator()( double t, double y ) const
   {
      return y;
   }
};

template <typename R>
void do_test( R( *f ) )
{
   //try {
   //	(*f)();
   //	cout << "OK" << endl;
   //}
   //catch (std::exception& e) {
   //	cout << "caught exception: ";
   //	cout << e.what() << endl;
   //}
   (*f)();
}

using std::tuple;
using namespace CTR::Functions;
using Eigen::Vector2d;
typedef tuple < Tube< constant_fun< Vector2d > >,
   Tube< constant_fun< Vector2d > >,
   Tube< constant_fun< Vector2d > >,
   Tube< constant_fun< Vector2d > >> Cannula3;

struct Configuration3
{
   Eigen::Vector3d	PsiL;
   Eigen::Vector3d	Beta;
};

struct LoadedConfiguration3
{
   Eigen::Vector3d	PsiL;
   Eigen::Vector3d	Beta;
   Eigen::Vector3d   Ftip;
   Eigen::Vector3d	Ttip;
};

Configuration3	get_random_config_for_cannula( Cannula3 const& C, default_random_engine& eng )
{
   uniform_real_distribution<double> beta1_dist( -get<0>( C ).GetLength(), 0.0e-3 );
   double beta1 = beta1_dist( eng );

   double minBeta2 = std::max( -get<1>( C ).GetLength(), beta1 );
   double maxBeta2 = std::min( beta1 + get<0>( C ).GetLength() - get<1>( C ).GetLength(), 0.0 );
   uniform_real_distribution<double> beta2_dist( minBeta2, maxBeta2 );
   double beta2 = beta2_dist( eng );

   double minBeta3 = std::max( -get<2>( C ).GetLength(), beta2 );
   double maxBeta3 = std::min( beta2 + get<1>( C ).GetLength() - get<2>( C ).GetLength(), 0.0 );
   uniform_real_distribution<double> beta3_dist( minBeta3, maxBeta3 );
   double beta3 = beta3_dist( eng );

   Configuration3 q;
   q.Beta << beta1, beta2, beta3;

   uniform_real_distribution<double> psi_dist( -M_PI, M_PI );
   q.PsiL << psi_dist( eng ), psi_dist( eng ), psi_dist( eng );

   return q;
}
//
//LoadedConfiguration3	get_random_config_for_cannula_with_loads( Cannula3 const& C, default_random_engine& eng )
//{
//   uniform_real_distribution<double> beta1_dist( -get<0>( C ).GetLength(), 0.0e-3 );
//   double beta1 = beta1_dist( eng );
//
//   double minBeta2 = std::max( -get<1>( C ).GetLength(), beta1 );
//   double maxBeta2 = std::min( beta1 + get<0>( C ).GetLength() - get<1>( C ).GetLength(), 0.0 );
//   uniform_real_distribution<double> beta2_dist( minBeta2, maxBeta2 );
//   double beta2 = beta2_dist( eng );
//
//   double minBeta3 = std::max( -get<2>( C ).GetLength(), beta2 );
//   double maxBeta3 = std::min( beta2 + get<1>( C ).GetLength() - get<2>( C ).GetLength(), 0.0 );
//   uniform_real_distribution<double> beta3_dist( minBeta3, maxBeta3 );
//   double beta3 = beta3_dist( eng );
//
//   LoadedConfiguration3 q;
//   q.Beta << beta1, beta2, beta3;
//
//   uniform_real_distribution<double> psi_dist( -M_PI, M_PI );
//   q.PsiL << psi_dist( eng ), psi_dist( eng ), psi_dist( eng );
//
//   //normal distribution for forces
//   uniform_real_distribution<double> F_dist( -2.0, 2.0 );
//   q.Ftip( 0 ) = F_dist( eng );
//   q.Ftip( 1 ) = F_dist( eng );
//   q.Ftip( 2 ) = F_dist( eng );
//
//   q.Ttip = Eigen::Vector3d::Zero();
//
//   return q;
//}
//
//void test_smooth_step_function()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*        Smooth Step Function Test          *\n"
//      "*********************************************";
//   cout << msg << endl;
//   double left = 1; double right = 2;
//   double sw = 1; double width = 1;
//
//   smooth_step_function<double> step_fun( left, right, sw, width );
//   cout << step_fun( 1.5 ) << endl;
//   auto der_fun = get_derivative( step_fun );
//
//   Utility::linspace<double> x( 0.0, 2.0, 2000 );
//
//   vector<double> y, z;
//   transform( begin( x ), end( x ), back_inserter( y ), step_fun );
//   transform( begin( x ), end( x ), back_inserter( z ), der_fun );
//
//   ofstream fx( "testX.txt" );
//   fx.precision( 10 );
//   copy( begin( x ), end( x ), ostream_iterator<double>( fx, "\n" ) );
//
//   ofstream f( "test.txt" );
//   f.precision( 10 );
//   copy( begin( y ), end( y ), ostream_iterator<double>( f, "\n" ) );
//   ofstream f2( "test2.txt" );
//   f2.precision( 10 );
//   copy( begin( z ), end( z ), ostream_iterator<double>( f2, "\n" ) );
//}
//
//void test_identity_zero_functions()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*        Basic Function Tests               *\n"
//      "*********************************************";
//   cout << msg << endl;
//   cout << "testing identity function" << endl;
//   //Check identity function 
//   constant_fun<double>	id( 4.0 );
//   std::cout << "four = " << id() << endl; //default argument
//   constant_fun<Eigen::Vector2d>	id_v( Eigen::Vector2d::UnitX() + 4.0*Eigen::Vector2d::UnitY() );
//   cout << "id_v = " << endl << id_v() << endl;
//
//   cout << "testing zero function" << endl;
//   //Check zero function
//   zero_fun<double> zero = get_derivative( id );
//   std::cout << "Zero: " << zero( 0 ) << endl;
//   zero_fun<Eigen::Vector2d> zero_v = get_derivative( id_v );
//   cout << "zero_v = " << endl << zero_v() << endl;
//}
//
//void test_mollified_indicator_function()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*        Mollified Indicator Fun Test       *\n"
//      "*********************************************";
//   cout << msg << endl;
//   cout << "testing mollified indicator function" << endl;
//
//   using namespace CTR;
//   using namespace Utility;
//   using namespace std;
//   mollified_indicator_fun<double> f( 0.0, 5.0, 0.1 );
//
//   auto v = linspace<double>( 0.0, 5.0, 1000 );
//   std::vector<double> r;
//   std::transform( std::begin( v ), std::end( v ), std::back_inserter( r ), f );
//
//   ofstream file1( "testX_chi.txt" );
//   std::copy( std::begin( v ), std::end( v ), ostream_iterator<double>( file1, "\n" ) );
//
//   ofstream file2( "testY_chi.txt" );
//   std::copy( std::begin( r ), std::end( r ), ostream_iterator<double>( file2, "\n" ) );
//
//   auto df = get_derivative( f );
//
//   std::vector<double> dr;
//   transform( begin( v ), end( v ), back_inserter( dr ), df );
//
//   ofstream file3( "testdY_chi.txt" );
//   std::copy( begin( dr ), end( dr ), ostream_iterator<double>( file3, "\n" ) );
//
//   TInterval::IntervalList ilist;
//   Functions::add_intervals( f, ilist );
//   std::cout << "Intervals: " << std::endl;
//   std::cout << ilist << std::endl;
//}
//
//void test_tube_class()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*        Tube Class Tests                   *\n"
//      "*********************************************";
//   cout << msg << endl;
//   double sw = 1; double width = 1;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T1_type;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T2_type;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T3_type;
//   Vector<2>::type	left_k = Vector<2>::type::UnitX();
//   Vector<2>::type right_k = Vector<2>::type::UnitY();
//   smooth_step_function< Vector<2>::type > step_fun_k( left_k, right_k, sw, width );
//   T1_type T1 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto v = make_tuple( T1, T2, T3 );
//   cout << std::get<0>( v ).GetBendingStiffness( Mathematics::eval_pt( 0.0, -1.0, 1.0 ) ) * std::get<0>( v ).GetTorsionalCompliance( Mathematics::eval_pt( 0.0, -1.0, 1.0 ) ) << endl;
//}
//
//void test_get_precurvature()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*        Get Precurvature Test              *\n"
//      "*********************************************";
//   cout << msg << endl;
//
//   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;
//
//   constant_fun< Vector<2>::type >	k_fun( Eigen::Vector2d::UnitX() );
//
//   T1_type T1( 0.4, 0.1, k_fun, 1.0, 1.0 / 1.33 );
//   T2_type T2( 0.2, 0.1, k_fun, 1.0, 1.0 / 1.33 );
//   T3_type T3( 0.1, 0.05, k_fun, 1.0, 1.0 / 1.33 );
//
//   auto v = make_tuple( T1, T2, T3 );
//   cout << std::get<0>( v ).GetPrecurvature( 1.0 ) << endl;
//}
//
//void test_tagged_intervals()
//{
//   const char *msg = "*********************************************\n"
//      "*        Tagged Interval Tests              *\n"
//      "*********************************************";
//   cout << msg << endl;
//   TInterval::TaggedInterval ti1( 0.0, 1.0, TInterval::DENSE );
//   TInterval::TaggedInterval ti2( 0.5, 2, TInterval::NORMAL );
//   TInterval::TaggedInterval ti3( 1.5, 3.0, TInterval::DENSE );
//   TInterval::TaggedInterval ti4( 2.0, 10.0, TInterval::SPARSE );
//
//   TInterval::IntervalList IL;
//   IL.push_back( ti1 );
//   IL.push_back( ti2 );
//   IL.push_back( ti3 );
//   IL.push_back( ti4 );
//   TInterval::IntervalList DL = TInterval::resolve_to_disjoint( IL );
//   std::cout << DL << std::endl;
//}
//
//void test_rk4_integration()
//{
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           RK4 Integration Tests           *\n"
//      "*********************************************";
//   boost::io::ios_all_saver saver( cout );
//   cout << msg << endl;
//   Utility::CodeTimer ct;
//
//   long N = 5;
//
//   Utility::linspace<double> times( 0.0, 1.0, N );
//   std::vector<double >	R;
//   R.reserve( times.size() );
//   ct.start();
//   Mathematics::integrate( times.begin(),
//                           times.end(),
//                           back_inserter( R ),
//                           1.0,
//                           Mathematics::rk4_step( test_derivative() )
//                           );
//   ct.stop();
//
//   std::copy( R.begin(), R.end(), std::ostream_iterator<double>( cout, "\n" ) );
//   std::cout << "Elapsed time per RK step: " << ct.elapsed() / static_cast<double>(N) << std::endl;
//
//   std::cout << "RK8_final test" << std::endl;
//   ct.start();
//   double e = Mathematics::rk8_final( test_derivative(), 1.0, Utility::linspace<double>( 0.0, 1.0, N ) );
//   ct.stop();
//
//   cout.precision( 20 );
//   cout << e << endl;
//   cout << "error: " << exp( 1.0 ) - e << std::endl;
//   std::cout << "Elapsed time per RK step: " << ct.elapsed() / static_cast<double>(N) << std::endl;
//}
//
//struct times_two
//{
//   template <typename T>
//   auto operator()( T const& t ) const -> decltype(t)
//   {
//      return 2.0*t;
//   }
//};
//
//void test_tuple_transform()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           Tuple Transform Tests           *\n"
//      "*********************************************";
//   std::cout << msg << endl;
//   double sw = 1; double width = 1;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T1_type;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T2_type;
//   typedef Tube< smooth_step_function< Vector<2>::type> >	T3_type;
//   Vector<2>::type	left_k = Vector<2>::type::UnitX();
//   Vector<2>::type right_k = Vector<2>::type::UnitY();
//   smooth_step_function< Vector<2>::type > step_fun_k( left_k, right_k, sw, width );
//   T1_type T1 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 1.0, 0.1, 1.0e-3, 0.8e-3, step_fun_k, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto v = std::make_tuple( T1, T2, T3 );
//
//   //this provides a view to the precurvature of each tube
//   auto view_kb = tuple_transform::transform( v, internal::evaluate_kb_at( Mathematics::eval_pt( 0.0, -1.0, 1.0 ) ) );
//   auto view_ct = tuple_transform::transform( v, internal::evaluate_ct_at( Mathematics::eval_pt( 0.0, -1.0, 1.0 ) ) );
//   auto view_u = tuple_transform::transform( v, internal::evaluate_precurvature_at( 0.0 ) );
//
//   std::tuple<double, double, double>	psi = std::make_tuple( 0.0, M_PI / 4.0, M_PI / 2.0 );
//   auto view_ru = tuple_transform::transform( view_u, psi, internal::rotation_about_z() );
//   std::cout << tuple_transform::at<0>( view_kb ) * tuple_transform::at<1>( view_ct ) << std::endl;
//   std::cout << "view_u: " << std::endl;
//   std::cout << tuple_transform::at<1>( view_u ) << std::endl;
//   std::cout << "unrotated: " << std::endl;
//   std::cout << tuple_transform::at<0>( view_ru ) << std::endl;
//   std::cout << "pi/4: " << std::endl;
//   std::cout << tuple_transform::at<1>( view_ru ) << std::endl;
//   std::cout << "pi/2: " << std::endl;
//   std::cout << tuple_transform::at<2>( view_ru ) << std::endl;
//
//   Eigen::Vector3d	kb_vec;
//   tuple_transform::copy( view_kb, kb_vec );
//   std::cout << "kb_vector: " << std::endl << kb_vec << std::endl;
//
//   Eigen::Matrix<double, 2, 3> u_mat;
//
//   tuple_transform::copy( view_ru, tuple_transform::columns_of( u_mat ) );
//
//   std::cout << "u_mat = " << std::endl;
//   std::cout << u_mat << std::endl;
//}
//
void test_kinematics_geom_only()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (geom)        *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun( 5.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.5, 0.02, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.3, 0.1, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
   } q;

   q.Beta << -0.150, -0.05, -0.025;
   q.PsiL << 0, -M_PI * 0.5, M_PI * 0.5;

   typedef CTR::DeclareOptions< Option::ComputeGeometry >::options	OType;

   //total_integration_time.QuadPart = 0;
   //total_setup_time.QuadPart = 0;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret;
   long N = 100000;
   for (int i = 0; i < N; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 50, 0 );
   }
   ct.stop();
   std::cout << "kinematics rate: " << static_cast<double>(N) / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;
   std::cout << "nmult = " << g_nmult << std::endl;
   std::cout << "nadd = " << g_nadd << std::endl;
   std::cout << "ndiv = " << g_ndiv << std::endl;
   std::cout << "ntrig = " << g_ntrig << std::endl;
   std::cout << "nsqrt = " << g_nsqrt << std::endl;
   std::cout << "nnorm = " << g_nnorm << std::endl;
   //LARGE_INTEGER freq;
   //QueryPerformanceFrequency( &freq );
   //std::cout << "ratio of setup time to integration time:    1:" << double( total_integration_time.QuadPart ) / double( total_setup_time.QuadPart ) << endl;
}

void test_kinematics_torsion_only()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (torsion)     *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun( 5.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.5, 0.2, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.3, 0.1, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
   } q;

   q.Beta << -0.150, -0.05, -0.025;
   q.PsiL << 0, -M_PI * 0.5, M_PI * 0.5;

   typedef CTR::DeclareOptions< Option::NoOption >::options	OType;

   //total_integration_time.QuadPart = 0;
   //total_setup_time.QuadPart = 0;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret;
   long N = 100000;
   for (int i = 0; i < N; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 0, 0 );
   }
   ct.stop();
   std::cout << "kinematics rate: " << static_cast<double>(N) / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;
   std::cout << "nmult = " << g_nmult << std::endl;
   std::cout << "nadd = " << g_nadd << std::endl;
   std::cout << "ndiv = " << g_ndiv << std::endl;
   std::cout << "ntrig = " << g_ntrig << std::endl;
   std::cout << "nsqrt = " << g_nsqrt << std::endl;
   std::cout << "nnorm = " << g_nnorm << std::endl;
   //LARGE_INTEGER freq;
   //QueryPerformanceFrequency( &freq );
   //std::cout << "ratio of setup time to integration time:    1:" << double( total_integration_time.QuadPart ) / double( total_setup_time.QuadPart ) << endl;
}

void test_kinematics_Zga()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (Zga)         *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun( 5.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.5, 0.2, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.3, 0.1, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
   } q;

   q.Beta << -0.150, -0.05, -0.025;
   q.PsiL << 0, -M_PI * 0.5, M_PI * 0.5;

   typedef CTR::DeclareOptions< Option::ComputeZga >::options	OType;

   //total_integration_time.QuadPart = 0;
   //total_setup_time.QuadPart = 0;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret;
   long N = 100000;
   for (int i = 0; i < N; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 0, 0 );
   }
   ct.stop();
   std::cout << "kinematics rate: " << static_cast<double>(N) / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;
   std::cout << "nmult = " << g_nmult << std::endl;
   std::cout << "nadd = " << g_nadd << std::endl;
   std::cout << "ndiv = " << g_ndiv << std::endl;
   std::cout << "ntrig = " << g_ntrig << std::endl;
   std::cout << "nsqrt = " << g_nsqrt << std::endl;
   std::cout << "nnorm = " << g_nnorm << std::endl;
   //LARGE_INTEGER freq;
   //QueryPerformanceFrequency( &freq );
   //std::cout << "ratio of setup time to integration time:    1:" << double( total_integration_time.QuadPart ) / double( total_setup_time.QuadPart ) << endl;
}

void test_kinematics_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (loads)       *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun( 5.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.5, 0.2, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.3, 0.1, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
      Vector<3>::type   Ftip;
      Vector<3>::type   Ttip;
   } q;

   q.Beta << -0.150, -0.05, -0.025;
   q.PsiL << 0, -M_PI * 0.5, M_PI * 0.5;
   q.Ftip << -0.25, 0.25, 0.3;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ExternalLoads >::options	OType;

   //total_integration_time.QuadPart = 0;
   //total_setup_time.QuadPart = 0;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret;
   long N = 100000;
   for (int i = 0; i < N; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 0, 0 );
   }
   ct.stop();
   std::cout << "kinematics rate: " << static_cast<double>(N) / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;
   std::cout << "nmult = " << g_nmult << std::endl;
   std::cout << "nadd = " << g_nadd << std::endl;
   std::cout << "ndiv = " << g_ndiv << std::endl;
   std::cout << "ntrig = " << g_ntrig << std::endl;
   std::cout << "nsqrt = " << g_nsqrt << std::endl;
   std::cout << "nnorm = " << g_nnorm << std::endl;
   //LARGE_INTEGER freq;
   //QueryPerformanceFrequency( &freq );
   //std::cout << "ratio of setup time to integration time:    1:" << double( total_integration_time.QuadPart ) / double( total_setup_time.QuadPart ) << endl;
}

void test_kinematics_jacobians_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (jac/loads)   *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun( 5.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.5, 0.2, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.3, 0.1, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
      Vector<3>::type   Ftip;
      Vector<3>::type   Ttip;
   } q;

   q.Beta << -0.150, -0.05, -0.025;
   q.PsiL << 0, -M_PI * 0.5, M_PI * 0.5;
   q.Ftip << -0.25, 0.25, 0.3;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeJacobian, Option::ComputeCompliance >::options	OType;

   //total_integration_time.QuadPart = 0;
   //total_setup_time.QuadPart = 0;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret;
   long N = 100000;
   for (int i = 0; i < N; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 0, 0 );
   }
   ct.stop();
   std::cout << "kinematics rate: " << static_cast<double>(N) / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;
   std::cout << "nmult = " << g_nmult << std::endl;
   std::cout << "nadd = " << g_nadd << std::endl;
   std::cout << "ndiv = " << g_ndiv << std::endl;
   std::cout << "ntrig = " << g_ntrig << std::endl;
   std::cout << "nsqrt = " << g_nsqrt << std::endl;
   std::cout << "nnorm = " << g_nnorm << std::endl;
   //LARGE_INTEGER freq;
   //QueryPerformanceFrequency( &freq );
   //std::cout << "ratio of setup time to integration time:    1:" << double( total_integration_time.QuadPart ) / double( total_setup_time.QuadPart ) << endl;
}


void test_Zta_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Test Zta with Loads   (Zta)     *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;


   q.Beta << -0.2, -0.16, -0.075;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 1.0, 1.0, 1.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeZta, Option::ExternalLoads >::options	OType;
   auto kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0.0 / 0.001 );

   std::cout << "Zta (integrated): " << std::endl;
   std::cout << kinRet.y_final.Zta << std::endl << std::endl;

   Eigen::Matrix<double, 6, 6>  Zta = Eigen::Matrix<double, 6, 6>::Zero();

   double h = 1e-7;
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 0 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 0 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 1 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 1 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 2 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 2 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }

   h = 1e-6;
   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 3 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 3 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
      Zta( 0, 3 ) = Zta( 0, 3 ) - kinRet_p.y_final.Mz( 0 ) * std::get<0>( cannula ).GetTorsionalCompliance( Mathematics::eval_pt( 0.0, 0.0, 0.001 ) );
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 4 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 4 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zta.block<3, 1>( 3, 5 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Zta.block<3, 1>( 0, 5 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }

   std::cout << "Zta (fd): " << std::endl;
   std::cout << Zta << std::endl << std::endl;

   std::cout << "Error: " << std::endl;
   std::cout << Zta - kinRet.y_final.Zta;
   std::cout << std::endl << std::endl;

}

void test_Ztf_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Test Ztf with Loads   (Ztf)     *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;


   q.Beta << -0.2, -0.16, -0.075;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 1.0, 1.0, 1.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeZtf, Option::ExternalLoads >::options	OType;
   auto kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0.0 / 0.001 );

   std::cout << "Ztf (integrated): " << std::endl;
   std::cout << kinRet.y_final.Ztf << std::endl << std::endl;

   Eigen::Matrix<double, 6, 6>  Ztf = Eigen::Matrix<double, 6, 6>::Zero();

   double h = 1e-7;
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 0 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 0 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 1 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 1 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 2 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 2 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }

   h = 1e-6;
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 3 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 3 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
      //Zta( 0, 3 ) = Zta( 0, 3 ) - kinRet_p.y_final.Mz( 0 ) * std::get<0>( cannula ).GetTorsionalCompliance( Mathematics::eval_pt( 0.0, 0.0, 0.001 ) );
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 4 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 4 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Ztf.block<3, 1>( 3, 5 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Mz - kinRet_n.y_final.Mz);
      Ztf.block<3, 1>( 0, 5 ) = 1.0 / h / 2.0*(kinRet_p.y_final.Psi - kinRet_n.y_final.Psi);
   }

   std::cout << "Ztf (fd): " << std::endl;
   std::cout << Ztf << std::endl << std::endl;

   std::cout << "Error: " << std::endl;
   std::cout << Ztf - kinRet.y_final.Ztf;
   std::cout << std::endl << std::endl;

}

void test_Zlf_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests  (Zlf)         *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;

   q.Beta << -0.2, -0.16, -0.075;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 1.0, 1.0, 1.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeZta, Option::ExternalLoads, Option::ComputeZlf >::options	OType;
   auto kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0 );

   std::cout << "Zlf (integrated): " << std::endl;
   std::cout << kinRet.y_final.Zlf << std::endl << std::endl;

   Eigen::Matrix<double, 6, 6>  Zlf = Eigen::Matrix<double, 6, 6>::Zero();

   cout << kinRet.y_final.p.cross( q.Ftip ) << endl;

   double h = 1e-7;
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zlf.block<3, 1>( 0, 0 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 0 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zlf.block<3, 1>( 0, 1 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 1 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zlf.block<3, 1>( 0, 2 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 2 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }

   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zlf.block<3, 1>( 0, 3 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 3 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zlf.block<3, 1>( 0, 4 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 4 ) = 1.0 / h / 2.0* (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 0, 0 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 0, 0 / 0.001 );

      Zlf.block<3, 1>( 0, 5 ) = 1.0 / h / 2.0 * (kinRet_p.y_final.Nt - kinRet_n.y_final.Nt);
      Zlf.block<3, 1>( 3, 5 ) = 1.0 / h / 2.0 * (kinRet_p.y_final.Mt - kinRet_n.y_final.Mt);
   }

   std::cout << "Zlf (fd): " << std::endl;
   std::cout << Zlf << std::endl << std::endl;

   std::cout << "Error: " << std::endl;
   std::cout << Zlf - kinRet.y_final.Zlf;
   std::cout << std::endl << std::endl;

}

void test_Zgf_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests (Zgf)         *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;

   q.Beta << -0.2, -0.16, -0.075;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 1.0, 1.0, 1.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeZta, Option::ExternalLoads, Option::ComputeZgf >::options	OType;
   auto kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0.0 / 0.001 );

   std::cout << "Zgf (integrated): " << std::endl;
   Eigen::Matrix<double, 6, 6>  Ad = Utility::Adjoint_p_q( kinRet.y_final.p, kinRet.y_final.q );
   std::cout << -Ad*kinRet.y_final.Zgf << std::endl << std::endl;

   Eigen::Matrix<double, 6, 6>  Zgf = Eigen::Matrix<double, 6, 6>::Zero();
   using Utility::quaternion_conj;
   using Utility::quaternion_mult;
   using Utility::quaternion_rotate;

   double h = 1e-8;
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zgf.block<3, 1>( 0, 0 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 0 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zgf.block<3, 1>( 0, 1 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 1 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ftip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ftip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zgf.block<3, 1>( 0, 2 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 2 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }

   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zgf.block<3, 1>( 0, 3 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 3 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zgf.block<3, 1>( 0, 4 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 4 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Ttip += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 0, 0 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Ttip -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 0, 0 / 0.001 );

      Zgf.block<3, 1>( 0, 5 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zgf.block<3, 1>( 3, 5 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }

   std::cout << "Zgf (fd): " << std::endl;
   std::cout << Zgf << std::endl << std::endl;

   std::cout << "Error: " << std::endl;
   std::cout << -Zgf - Ad*kinRet.y_final.Zgf;
   std::cout << std::endl << std::endl;

}

void test_Zga_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests (Zga w/loads) *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;

   q.Beta << -0.2, -0.16, -0.075;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 1.0, 1.0, 1.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeZta, Option::ExternalLoads, Option::ComputeZgf >::options	OType;
   auto kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0.0 / 0.001 );

   std::cout << "Zga (integrated): " << std::endl;
   Eigen::Matrix<double, 6, 6>  Ad = Utility::Adjoint_p_q( kinRet.y_final.p, kinRet.y_final.q );
   std::cout << -Ad*kinRet.y_final.Zga << std::endl << std::endl;

   Eigen::Matrix<double, 6, 6>  Zga = Eigen::Matrix<double, 6, 6>::Zero();
   using Utility::quaternion_conj;
   using Utility::quaternion_mult;
   using Utility::quaternion_rotate;

   double h = 1e-8;
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zga.block<3, 1>( 0, 0 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 0 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zga.block<3, 1>( 0, 1 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 1 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.PsiL += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zga.block<3, 1>( 0, 2 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 2 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }

   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitX();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitX();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zga.block<3, 1>( 0, 3 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 3 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitY();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 50, 2 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitY();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 50, 2 / 0.001 );

      Zga.block<3, 1>( 0, 4 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 4 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }
   {
      LoadedConfiguration3 qP = q;
      qP.Beta += h*Eigen::Vector3d::UnitZ();
      auto kinRet_p = CTR::Kinematics( cannula, qP, OType(), 0, 0 / 0.001 );
      LoadedConfiguration3 qN = q;
      qN.Beta -= h*Eigen::Vector3d::UnitZ();
      auto kinRet_n = CTR::Kinematics( cannula, qN, OType(), 0, 0 / 0.001 );

      Zga.block<3, 1>( 0, 5 ) = 1.0 / h / 2.0* quaternion_rotate( quaternion_conj( kinRet.qTip ), (kinRet_p.pTip - kinRet_n.pTip) );
      Zga.block<3, 1>( 3, 5 ) = 2.0 / h / 2.0* quaternion_mult( quaternion_conj( kinRet.qTip ), (kinRet_p.qTip - kinRet_n.qTip) ).bottomRows<3>();
   }

   std::cout << "Zga (fd): " << std::endl;
   std::cout << Zga << std::endl << std::endl;

   std::cout << "Error: " << std::endl;
   std::cout << -Zga - Ad*kinRet.y_final.Zga;
   std::cout << std::endl << std::endl;

}

void test_kinematics_with_loads()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests                *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10 * Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.35, 0.3, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.25, 0.2, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.15, 0.1, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   LoadedConfiguration3 q;

   q.Beta << -0.2, -0.15, -0.075;
   q.PsiL << 0.0, 0.0, 0.0;
   //q.PsiL << 0.0, 0.0, 0.0;
   q.Ftip << 2.0, -2.0, 0.0;
   q.Ttip << 0.0, 0.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian, Option::ExternalLoads, Option::ComputeCompliance >::options	OType;

   std::default_random_engine	eng;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret, ret2;

   for (int i = 0; i < 1000; ++i) {
      ret = CTR::Kinematics( cannula, q, OType(), 0 );
   }
   ct.stop();
   ret2 = CTR::Kinematics( cannula, q, OType(), 200, 10.0 / 0.001 );

   std::cout << "position error: " << (ret.pTip - ret2.pTip).norm() << std::endl;
   std::cout << "quaternion mag error: " << ret.qTip.norm() - 1 << std::endl;
   std::cout << "quaternion angle err: " << Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( ret.qTip ), ret2.qTip ) ) << std::endl;

   //Eigen::Matrix<double, 6, 6>	JJ = CTR::GetTipJacobianForTube1(ret2.y_final);
   std::cout << "kinematics rate: " << 1000.0 / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret2.y_final << std::endl;

   using namespace Utility;
   std::cout << "position: " << ret2.pTip.transpose() << std::endl;
   std::cout << "orientation: " << std::endl << quaternion_to_rotation_matrix( ret2.qTip ) << std::endl;

   Eigen::Matrix<double, 6, 6> Zta = Eigen::Matrix<double, 6, 6>::Zero();
   Eigen::Matrix<double, 6, 6> Zla = Eigen::Matrix<double, 6, 6>::Zero();
   double h = 1e-7;
   {
      LoadedConfiguration3 qP = q;
      Eigen::Vector3d dir = Eigen::Vector3d::UnitX();
      qP.PsiL += h * dir;
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h * dir;

      auto retP = CTR::Kinematics( cannula, qP, OType(), 100, 10.0 / 0.001 );
      auto retN = CTR::Kinematics( cannula, qN, OType(), 100, 10.0 / 0.001 );

      Zta.col( 0 ).topRows<3>() = (retP.y_final.Psi - retN.y_final.Psi) / 2.0 / h;
      Zta.col( 0 ).bottomRows<3>() = (retP.y_final.Mz - retN.y_final.Mz) / 2.0 / h;

      Zla.col( 0 ).topRows<3>() = (retP.y_final.Nt - retN.y_final.Nt) / 2.0 / h;
      Zla.col( 0 ).bottomRows<3>() = (retP.y_final.Mt - retN.y_final.Mt) / 2.0 / h;
   }
   {
      LoadedConfiguration3 qP = q;
      Eigen::Vector3d dir = Eigen::Vector3d::UnitY();
      qP.PsiL += h * dir;
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h * dir;

      auto retP = CTR::Kinematics( cannula, qP, OType(), 100, 10.0 / 0.001 );
      auto retN = CTR::Kinematics( cannula, qN, OType(), 100, 10.0 / 0.001 );

      Zta.col( 1 ).topRows<3>() = (retP.y_final.Psi - retN.y_final.Psi) / 2.0 / h;
      Zta.col( 1 ).bottomRows<3>() = (retP.y_final.Mz - retN.y_final.Mz) / 2.0 / h;

      Zla.col( 1 ).topRows<3>() = (retP.y_final.Nt - retN.y_final.Nt) / 2.0 / h;
      Zla.col( 1 ).bottomRows<3>() = (retP.y_final.Mt - retN.y_final.Mt) / 2.0 / h;
   }
   {
      LoadedConfiguration3 qP = q;
      Eigen::Vector3d dir = Eigen::Vector3d::UnitZ();
      qP.PsiL += h * dir;
      LoadedConfiguration3 qN = q;
      qN.PsiL -= h * dir;

      auto retP = CTR::Kinematics( cannula, qP, OType(), 100, 10.0 / 0.001 );
      auto retN = CTR::Kinematics( cannula, qN, OType(), 100, 10.0 / 0.001 );

      Zta.col( 2 ).topRows<3>() = (retP.y_final.Psi - retN.y_final.Psi) / 2.0 / h;
      Zta.col( 2 ).bottomRows<3>() = (retP.y_final.Mz - retN.y_final.Mz) / 2.0 / h;

      Zla.col( 2 ).topRows<3>() = (retP.y_final.Nt - retN.y_final.Nt) / 2.0 / h;
      Zla.col( 2 ).bottomRows<3>() = (retP.y_final.Mt - retN.y_final.Mt) / 2.0 / h;
   }

   cout.precision( 4 );
   std::cout << "Zta (fd): " << std::endl << Zta.leftCols<3>() << std::endl << std::endl;
   std::cout << "Zta (prop): " << std::endl << ret2.y_final.Zta.leftCols<3>() << std::endl << std::endl;

   std::cout << "Zla (fd): " << std::endl << Zla.leftCols<3>() << std::endl << std::endl;
   std::cout << "Zla (prop): " << std::endl << ret2.y_final.Zla.leftCols<3>() << std::endl << std::endl;

   cout.precision( 10 );
   std::cout << "Moment sum: " << ret2.y_final.Mz.sum() << std::endl;
   std::cout << "r cross F: " << std::endl << (ret2.pTip - q.Beta( 0 )*Eigen::Vector3d::UnitZ()).cross( quaternion_to_rotation_matrix( ret2.qTip ) * q.Ftip ) << std::endl;
   std::cout << "qtip*Mt*inv(qtip): " << std::endl << quaternion_rotate( ret2.qTip, ret2.y_final.Mt ) << std::endl;
   {
      boost::io::ios_all_saver saver( std::cout );
      std::cout.precision( 6 );

      std::cout << "stability: " << CTR::GetStability( ret2.y_final ) << std::endl;
   }
   std::cout.precision( 3 );
   Eigen::Matrix<double, 6, 6> JJ = CTR::GetTipJacobianForTube1( ret.y_final );
   std::cout << "Zga transformed: " << std::endl << JJ << std::endl;

   Eigen::Matrix<double, 6, 6> J = CTR::GetLoadedFDJacobianForTube1( cannula, q, 1.0e-7, 1.0e-6 );
   std::cout << "Jb (FD) " << std::endl << J << std::endl;

   std::cout.precision( 3 );
   std::cout << "error : " << std::endl << (JJ - J) << std::endl;
   std::cout << "error norm of positional part (F-norm) : " << (JJ - J).topRows<3>().norm() << std::endl;

   std::cout << "Compliance Matrix: " << std::endl;
   std::cout << CTR::GetTipComplianceForTube1( ret2.y_final ) << std::endl;


}

void test_kinematics()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests                *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 59.886827696955997*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 18.703425486187498*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 44.894837830999300*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.376, 0.3622, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.2778, 0.2328, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1602, 0.1336, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   Configuration3 q;

   q.Beta << -0.348140191074130, -0.266418404773888, -0.154306582220192;
   q.PsiL << 0.351222965206277, 0.752910086847227, 2.540274845786610;

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options	OType;

   std::default_random_engine	eng;
   Utility::CodeTimer ct;
   ct.start();
   KinRet< State<3, OType> >	ret, ret2;

   for (int i = 0; i < 100000; ++i) {
      //q = get_random_config_for_cannula(cannula, eng);
      ret = CTR::Kinematics( cannula, q, OType(), 0, 0 );
   }
   ct.stop();
   ret2 = CTR::Kinematics( cannula, q, OType(), 200, 30. / 0.001 );

   std::cout << "position error: " << (ret.pTip - ret2.pTip).norm() << std::endl;
   std::cout << "quaternion mag error: " << ret.qTip.norm() - 1 << std::endl;
   std::cout << "quaternion angle err: " << Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( ret.qTip ), ret2.qTip ) ) << std::endl;

   Eigen::Matrix<double, 6, 6>	JJ = CTR::GetTipJacobianForTube1( ret.y_final );
   std::cout << "kinematics rate: " << 100000.0 / ct.elapsed() << std::endl;
   std::cout << "y_final = " << std::endl << ret.y_final << std::endl;

   using namespace Utility;
   std::cout << "position: " << ret.pTip.transpose() << std::endl;
   std::cout << "orientation: " << std::endl << quaternion_to_rotation_matrix( ret.qTip ) << std::endl;

   {
      boost::io::ios_all_saver saver( std::cout );
      std::cout.precision( 6 );

      std::cout << "stability: " << CTR::GetStability( ret.y_final ) << std::endl;
   }
   std::cout << "Zga transformed: " << std::endl << JJ << std::endl;

   Eigen::Matrix<double, 6, 6> J = CTR::GetFDJacobianForTube1( cannula, q, 1.0e-6, 1.0e-6 );
   std::cout << "Jb (FD) " << std::endl << J << std::endl;

   std::cout.precision( 3 );
   std::cout << "error : " << std::endl << (JJ - J) << std::endl;
   std::cout << "error norm of positional part (F-norm) : " << (JJ - J).topRows<3>().norm() << std::endl;

   std::cout << "New FD Jacobian: " << std::endl;
}

void test_dense_output()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests                *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k1_fun( 10.*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10.*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 10.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 0.25, 0.15, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.175, 0.1, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1, 0.05, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );

   Configuration3 q;

   q.Beta << -0.15, -0.1, -0.05;
   q.PsiL << M_PI / 4.0, -M_PI / 4.0, 0.0;

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options	OType;


   KinRetDense< State<3, OType> >	ret;

   ret = CTR::Kinematics_with_dense_output( cannula, q, OType(), 0, 0 );

}

template <typename MatrixType>
double singular_value_error( const Eigen::MatrixBase<MatrixType> &M )
{
   Eigen::MatrixXd A = M;
   Eigen::JacobiSVD< Eigen::MatrixXd > svd( A );

   return svd.singularValues().maxCoeff();
}

void test_kinematics_dense_points_convergence()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Convergence Test (Dense)        *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;

   constant_fun< Vector<2>::type > k_fun1( 17.2134*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k_fun2( 17.2258*Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k_fun3( 45.4668*Eigen::Vector2d::UnitX() );

   //Parameters for make_annular_tube are
   // Length (total), Length (transmission), OD, ID, curvature function, E (Young's Modulus), G (Shear Modulus)
   T1_type T1 = make_annular_tube( 0.2785, 0.2378, 1.0e-3, 0.8e-3, k_fun1, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 0.2778, 0.2419, 1.0e-3, 0.8e-3, k_fun2, 50e9, 50e9 / 2.0 / 1.33 );
   T3_type T3 = make_annular_tube( 0.1841, 0.1582, 1.0e-3, 0.8e-3, k_fun3, 50e9, 50e9 / 2.0 / 1.33 );

   struct Configuration
   {
      Vector<3>::type	PsiL;
      Vector<3>::type	Beta;
   } q;

   q.Beta << -0.0516, -0.0515, -0.0049;
   q.PsiL << 2.946, -2.3437, 2.1049;

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options	OType;

   auto cannula = std::make_tuple( T1, T2, T3 );

   auto kinRet = CTR::Kinematics( cannula, q, OType(), 500, 200 / 0.001 );

   std::cout << kinRet.y_final << std::endl;

   auto kinRet1 = CTR::Kinematics( cannula, q, OType(), 500, 0. / 0.001 );
   auto kinRet2 = CTR::Kinematics( cannula, q, OType(), 500, 3. / 0.001 );
   auto kinRet3 = CTR::Kinematics( cannula, q, OType(), 500, 5. / 0.001 );
   auto kinRet4 = CTR::Kinematics( cannula, q, OType(), 500, 10. / 0.001 );
   auto kinRet5 = CTR::Kinematics( cannula, q, OType(), 500, 20. / 0.001 );
   auto kinRet6 = CTR::Kinematics( cannula, q, OType(), 500, 30. / 0.001 );
   auto kinRet7 = CTR::Kinematics( cannula, q, OType(), 500, 50. / 0.001 );
   auto kinRet8 = CTR::Kinematics( cannula, q, OType(), 500, 100. / 0.001 );
   auto kinRet9 = CTR::Kinematics( cannula, q, OType(), 500, 150. / 0.001 );
   auto kinRet10 = CTR::Kinematics( cannula, q, OType(), 500, 175. / 0.001 );

   double perr[10];
   perr[0] = (kinRet1.pTip - kinRet.pTip).norm();
   perr[1] = (kinRet2.pTip - kinRet.pTip).norm();
   perr[2] = (kinRet3.pTip - kinRet.pTip).norm();
   perr[3] = (kinRet4.pTip - kinRet.pTip).norm();
   perr[4] = (kinRet5.pTip - kinRet.pTip).norm();
   perr[5] = (kinRet6.pTip - kinRet.pTip).norm();
   perr[6] = (kinRet7.pTip - kinRet.pTip).norm();
   perr[7] = (kinRet8.pTip - kinRet.pTip).norm();
   perr[8] = (kinRet9.pTip - kinRet.pTip).norm();
   perr[9] = (kinRet10.pTip - kinRet.pTip).norm();
   std::ofstream f_perr( "err_p_over_dense.txt" );
   std::copy( std::begin( perr ), std::end( perr ), std::ostream_iterator<double>( f_perr, "\n" ) );

   double qerr[10];
   qerr[0] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet1.qTip ) );
   qerr[1] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet2.qTip ) );
   qerr[2] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet3.qTip ) );
   qerr[3] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet4.qTip ) );
   qerr[4] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet5.qTip ) );
   qerr[5] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet6.qTip ) );
   qerr[6] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet7.qTip ) );
   qerr[7] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet8.qTip ) );
   qerr[8] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet9.qTip ) );
   qerr[9] = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( kinRet.qTip ), kinRet10.qTip ) );
   std::ofstream f_qerr( "err_q_over_dense.txt" );
   std::copy( std::begin( qerr ), std::end( qerr ), std::ostream_iterator<double>( f_qerr, "\n" ) );

   double Jerr[10];
   Jerr[0] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet1.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[1] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet2.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[2] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet3.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[3] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet4.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[4] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet5.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[5] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet6.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[6] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet7.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[7] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet8.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[8] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet9.y_final )).block<3, 3>( 0, 0 ) );
   Jerr[9] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet10.y_final )).block<3, 3>( 0, 0 ) );
   std::ofstream f_Jerr_pa( "err_J_pa_over_dense.txt" );
   std::copy( std::begin( Jerr ), std::end( Jerr ), std::ostream_iterator<double>( f_Jerr_pa, "\n" ) );

   Jerr[0] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet1.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[1] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet2.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[2] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet3.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[3] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet4.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[4] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet5.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[5] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet6.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[6] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet7.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[7] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet8.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[8] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet9.y_final )).block<3, 3>( 3, 0 ) );
   Jerr[9] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet10.y_final )).block<3, 3>( 3, 0 ) );
   std::ofstream f_Jerr_wa( "err_J_wa_over_dense.txt" );
   std::copy( std::begin( Jerr ), std::end( Jerr ), std::ostream_iterator<double>( f_Jerr_wa, "\n" ) );

   Jerr[0] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet1.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[1] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet2.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[2] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet3.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[3] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet4.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[4] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet5.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[5] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet6.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[6] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet7.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[7] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet8.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[8] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet9.y_final )).block<3, 3>( 0, 3 ) );
   Jerr[9] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet10.y_final )).block<3, 3>( 0, 3 ) );
   std::ofstream f_Jerr_pt( "err_J_pt_over_dense.txt" );
   std::copy( std::begin( Jerr ), std::end( Jerr ), std::ostream_iterator<double>( f_Jerr_pt, "\n" ) );

   Jerr[0] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet1.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[1] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet2.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[2] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet3.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[3] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet4.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[4] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet5.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[5] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet6.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[6] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet7.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[7] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet8.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[8] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet9.y_final )).block<3, 3>( 3, 3 ) );
   Jerr[9] = singular_value_error( (CTR::GetTipJacobianForTube1( kinRet.y_final ) - CTR::GetTipJacobianForTube1( kinRet10.y_final )).block<3, 3>( 3, 3 ) );
   std::ofstream f_Jerr_wt( "err_J_wt_over_dense.txt" );
   std::copy( std::begin( Jerr ), std::end( Jerr ), std::ostream_iterator<double>( f_Jerr_wt, "\n" ) );

   std::cout << "Jacobian error: " << std::endl << CTR::GetTipJacobianForTube1( kinRet10.y_final ) - CTR::GetTipJacobianForTube1( kinRet.y_final ) << std::endl;
}

void test_two_tube_kinematics()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Kinematics Tests (two tube)     *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;

   constant_fun< Vector<2>::type > k_fun( 68.*Eigen::Vector2d::UnitX() );
   T1_type T1 = make_annular_tube( 35e-3, 5e-3, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );
   T2_type T2 = make_annular_tube( 20e-3, 0e-3, 1.0e-3, 0.8e-3, k_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2 );

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian, Option::ComputeCompliance >::options	OType;

   struct Configuration
   {
      Eigen::Vector2d PsiL;
      Eigen::Vector2d Beta;
      Eigen::Vector3d Ftip;
      Eigen::Vector3d Ttip;
   } q;

   q.Beta << -5e-3, 0.0;
   q.PsiL << 0.0, 0.0;
   KinRet<State<2, OType> >	kinRet;

   Utility::CodeTimer ct;
   ct.start();

   for (int i = 0; i < 1; ++i)
      kinRet = CTR::Kinematics( cannula, q, OType(), 0, 0 );

   ct.stop();

   std::cout << kinRet.y_final << std::endl;

   cout.precision( 3 );
   std::cout << "Jb = " << std::endl << CTR::GetTipJacobianForTube1( kinRet.y_final ) << std::endl;
   std::cout << "pTip = " << kinRet.pTip.transpose() << std::endl;
   std::cout << "qTip = " << kinRet.qTip.transpose() << std::endl;
   std::cout << "Rtip = " << std::endl << Utility::quaternion_to_rotation_matrix( kinRet.qTip ) << std::endl;
   std::cout << "y.Zga = " << std::endl << kinRet.y_final.Zga << std::endl;

   {
      boost::io::ios_all_saver saver( std::cout );
      std::cout.precision( 6 );
      std::cout << "Stability = " << CTR::GetStability( kinRet.y_final ) << std::endl;
   }
   std::cout << "Rate = " << 10000.0 / ct.elapsed() << std::endl;
   std::cout << "Jacobian (FD) = " << std::endl << CTR::GetFDJacobianForTube1( cannula, q, 1.0e-5, 1.0e-8 ) << std::endl;
   std::cout << "Jacobian error = " << std::endl << CTR::GetFDJacobianForTube1( cannula, q, 1.0e-5, 1.0e-8 ) - CTR::GetTipJacobianForTube1( kinRet.y_final ) << std::endl;
}

void test_options()
{
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Options Tests                   *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   typedef CTR::DeclareOptions < CTR::Option::ComputeJacobian,
      CTR::Option::ExternalLoads > ::options	OType;

   std::cout << "Option number: ";
   std::cout << std::hex << std::showbase << std::setw( 10 ) << std::setfill( '0' ) << std::internal << OType::value << std::endl;

   CTR::Option::put_options< OType >( cout );

   CTR::State<3, OType >	S;
   cout.flush();
   std::cout << "size of state: " << std::dec << sizeof( CTR::State<3, OType> ) << std::endl;
}

void test_function_intervals()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Function Interval Tests         *\n"
      "*********************************************";
   std::cout << msg << std::endl;

   using namespace CTR;
   using namespace std;

   smooth_step_function<double> f( 0.0, 1.0, 0.1, 0.001 );

   TInterval::IntervalList il;
   add_intervals( f, il );
	std::cout << "function intervals for smooth_step_function" << std::endl;
   std::cout << il << std::endl;

   mollified_indicator_fun<double> fChi( 0.0, 1.0, 0.001 );
   TInterval::IntervalList ilChi;
   add_intervals( fChi, ilChi );
	std::cout << "function intervals for mollified_indicator_fun" << std::endl;
   std::cout << ilChi << std::endl;
}

void test_smooth_step_integration()
{
   using namespace Functions;
   boost::io::ios_all_saver iosaver( std::cout );
   const char *msg = "*********************************************\n"
      "*           Smooth Step Integration         *\n"
      "*********************************************";
   std::cout << msg << std::endl;
   using namespace CTR;
   smooth_step_function<double> step( 0, 1.0, 0.5, 1.0 );
   smooth_step_function<double> step2( 0, 2.0, 0.5, 1.0 );

   smooth_step_fun_der<double> step1_der = get_derivative( step );

   double d = Mathematics::rk4_final(
      [step, step2, step1_der]( double s, double y ) {
      return (step1_der( s )*(step2( s ) + 1.0)) / pow( (step( s ) + step2( s ) + 1.0), 2.0 );
   },
      0.0,
      Utility::linspace<double>( 0.0, 1.0, 1000 )
      );
   std::cout << d << std::endl;
}

double deg2rad( double angle )
{
   return angle*M_PI / 180.0;
}


Cannula3	get_random_3tube_cannula( std::default_random_engine &eng )
{
   uniform_real_distribution<double>	angle_dist( deg2rad( 0 ), deg2rad( 90 ) );;
   double theta1 = angle_dist( eng );
   double theta2 = angle_dist( eng );
   double theta3 = angle_dist( eng );

   uniform_real_distribution<double>	Lc_dist( 10.0e-3, 50.0e-3 );
   double Lc1 = Lc_dist( eng );
   double Lc2 = Lc_dist( eng );
   double Lc3 = Lc_dist( eng );

   double k1 = theta1 / Lc1;
   double k2 = theta2 / Lc2;
   double k3 = theta3 / Lc3;

   Vector2d k1_vec = { k1, 0.0 };
   constant_fun< Vector2d > k1_fun( k1_vec );

   Vector2d k2_vec = { k2, 0.0 };
   constant_fun< Vector2d > k2_fun( k2_vec );

   Vector2d k3_vec = { k3, 0.0 };
   constant_fun< Vector2d > k3_fun( k3_vec );

   typedef Tube< constant_fun< Vector2d > >	T1_type;
   typedef Tube< constant_fun< Vector2d > >	T2_type;
   typedef Tube< constant_fun< Vector2d > >	T3_type;

   double clear = 0.1e-3;
   double walls[3];
   uniform_real_distribution<double> wall_dist( 0.1e-3, 0.3e-3 );
   generate( walls, walls + 3, bind( wall_dist, eng ) );

   double ID1 = 0.86e-3;
   double OD1 = ID1 + walls[0];
   double ID2 = OD1 + clear;
   double OD2 = ID2 + walls[1];
   double ID3 = OD2 + clear;
   double OD3 = ID3 + walls[2];

   uniform_real_distribution<double> length_dist_3( Lc3, 3*Lc3 );
   double L3 = length_dist_3( eng );

   double minL2 = std::max( L3, Lc2 );
   double maxL2 = std::max( L3, 4 * Lc2 );
   uniform_real_distribution<double> length_dist_2( minL2, maxL2 );
   double L2 = length_dist_2( eng );

   double minL3 = std::max( L2, Lc1 );
   double maxL3 = std::max( L2, 5 * Lc1 );
   uniform_real_distribution<double>	length_dist_1( minL3, maxL3 );
   double L1 = length_dist_1( eng );

   T1_type	tube_1 = make_annular_tube( L1, L1 - Lc1, OD1, ID1, k1_fun, 50e9, 50e9 / 2. / 1.33 );
   T2_type tube_2 = make_annular_tube( L2, L2 - Lc2, OD2, ID2, k2_fun, 50e9, 50e9 / 2. / 1.33 );
   T2_type tube_3 = make_annular_tube( L3, L3 - Lc3, OD3, ID3, k3_fun, 50e9, 50e9 / 2. / 1.33 );

   return Cannula3( tube_1, tube_2, tube_3 );
}

struct min_max_mean
{
   min_max_mean() :
      min( std::numeric_limits<double>::max() ),
      max( std::numeric_limits<double>::min() ),
      mean( 0 ),
      N( 0 )
   {}

   void operator()( double item )
   {
      if (item > max)
         max = item;
      if (item < min)
         min = item;

      mean = 1.0 / static_cast<double>(N + 1) * (N * mean + item);
      N++;
   }

   min_max_mean operator+( min_max_mean const& other ) const
   {
      min_max_mean R;
      R.min = (other.min < this->min) ? other.min : this->min;
      R.max = (other.max > this->max) ? other.max : this->max;
      R.mean = 1.0 / (this->N + other.N) * (this->N * this->mean + other.N * other.mean);
      return R;
   }

   double min;
   double max;
   double mean;
   int N;
};

struct stats
{
   double J_va_err;
   double J_vt_err;
   double J_wa_err;
   double J_wt_err;
   double p_err;
   double q_err;
   double psi_err;
   double mz_err;
   int n_normal_points;
   double dense_density;

   double psi1;
   double psi2;
   double psi3;
   double beta1;
   double beta2;
   double beta3;
   double L1;
   double L2;
   double L3;
   double k1;
   double k2;
   double k3;
   double Lt1;
   double Lt2;
   double Lt3;
   double kb1;
   double kb2;
   double kb3;
   double ct1, ct2, ct3;
   double Fx, Fy, Fz;

   double stability;
   double q_mag_err;
   double q_angle_err;
   double psi_diff_mag;
   double mz_mag;
   double q_err_x;
   double q_err_y;
   double q_err_z;
   double normal_density;
   double normal_length;
   int normal_points_used;

   //loaded stats
   double tip_deflection;
};

struct stats_min_max_mean
{
   min_max_mean	p_err;
   min_max_mean	q_err;
   min_max_mean	psi_diff_mag;
   min_max_mean	mz_diff_mag;
};

stats_min_max_mean	operator+( stats_min_max_mean const& left, stats_min_max_mean const& right )
{
   stats_min_max_mean	s;
   s.p_err = left.p_err + right.p_err;
   s.q_err = left.q_err + right.q_err;
   s.mz_diff_mag = left.mz_diff_mag + right.mz_diff_mag;
   s.psi_diff_mag = left.psi_diff_mag + right.psi_diff_mag;
   return s;
}

ostream& operator<<( ostream& os, stats const& s )
{
   boost::io::ios_all_saver	saver( os );

   os.precision( 16 );
   os.width( 18 );

   //os << s.n_normal_points << " "; //1
   //os << s.dense_density << " ";   //2
   os << s.L1 << " ";              //3
   os << s.L2 << " ";            //4
   os << s.L3 << " ";             //5
   os << s.Lt1 << " ";             //6
   os << s.Lt2 << " ";             //7
   os << s.Lt3 << " ";             //8
   os << s.k1 << " ";             //9
   os << s.k2 << " ";             //10
   os << s.k3 << " ";             //11
   os << s.psi1 << " ";             //12
   os << s.psi2 << " ";             //13
   os << s.psi3 << " ";             //14
   os << s.beta1 << " ";             //15
   os << s.beta2 << " ";             //16
   os << s.beta3 << " ";             //17
   os << s.p_err << " ";             //18
   os << s.q_err << " ";             //19
   os << s.psi_err << " ";             //20
   os << s.mz_err << " ";             //21
   os << s.J_va_err << " ";             //22
   os << s.J_vt_err << " ";             //23
   os << s.J_wa_err << " ";             //24
   os << s.J_wt_err << " ";             //25
   os << s.kb1 << " ";             //26
   os << s.kb2 << " ";             //27
   os << s.kb3 << " ";             //28
   os << s.ct1 << " ";             //29
   os << s.ct2 << " ";             //30
   os << s.ct3 << " ";             //31
   os << s.stability << " ";       //32
   os << s.q_angle_err << " ";     //33
   os << s.q_mag_err << " ";       //34
   os << s.psi_diff_mag << " ";    //35
   os << s.mz_mag << " ";          //36
   //os << s.q_err_x << " ";         //37
   //os << s.q_err_y << " ";         //38
   //os << s.q_err_z << " ";         //39
   os << s.normal_density << " ";  //40
   os << s.normal_length << " ";   //41
   os << s.Fx << " ";              //42
   os << s.Fy << " ";              //43
   os << s.Fz << " ";              //44
   os << s.tip_deflection << " ";
   os << s.normal_points_used << " ";

   return os;
}



template <typename SType>
stats get_stats_from_rets( KinRet<SType> const& ret_test, KinRet<SType> const& ret2, int n_normal_points, double dense_density, Cannula3 const& C, Configuration3 const& q )
{
   stats s;

   s.n_normal_points = n_normal_points;
   s.dense_density = dense_density;
   s.p_err = (ret_test.pTip - ret2.pTip).norm();
   s.q_err =
      s.mz_err = (ret_test.y_final.Mz - ret2.y_final.Mz).norm();
   s.psi_err = (ret_test.y_final.Psi - ret2.y_final.Psi).norm();

   s.J_va_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 0, 0 ) );
   s.J_vt_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 0, 3 ) );
   s.J_wa_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 3, 0 ) );
   s.J_wt_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 3, 3 ) );

   s.beta1 = q.Beta( 0 );
   s.beta2 = q.Beta( 1 );
   s.beta3 = q.Beta( 2 );
   s.psi1 = q.PsiL( 0 );
   s.psi2 = q.PsiL( 1 );
   s.psi3 = q.PsiL( 2 );
   s.L1 = std::get<0>( C ).GetLength();
   s.L2 = std::get<1>( C ).GetLength();
   s.L3 = std::get<2>( C ).GetLength();
   s.Lt1 = std::get<0>( C ).GetTransmissionLength();
   s.Lt2 = std::get<1>( C ).GetTransmissionLength();
   s.Lt3 = std::get<2>( C ).GetTransmissionLength();

   using Mathematics::eval_pt;
   eval_pt pt1( (s.Lt1 + s.L1) / 2.0, s.Lt1, s.L1 );
   eval_pt pt2( (s.Lt2 + s.L2) / 2.0, s.Lt2, s.L2 );
   eval_pt pt3( (s.Lt3 + s.L3) / 2.0, s.Lt3, s.L3 );
   s.k1 = std::get<0>( C ).GetPrecurvature( pt1.s )(0);
   s.k2 = std::get<1>( C ).GetPrecurvature( pt2.s )(0);
   s.k3 = std::get<2>( C ).GetPrecurvature( pt3.s )(0);

   s.kb1 = std::get<0>( C ).GetBendingStiffness( pt1 );
   s.kb2 = std::get<1>( C ).GetBendingStiffness( pt2 );
   s.kb3 = std::get<2>( C ).GetBendingStiffness( pt3 );

   s.ct1 = std::get<0>( C ).GetTorsionalCompliance( pt1 );
   s.ct2 = std::get<1>( C ).GetTorsionalCompliance( pt2 );
   s.ct3 = std::get<2>( C ).GetTorsionalCompliance( pt3 );

   s.stability = CTR::GetStability( ret2.y_final );
   s.q_angle_err = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( ret_test.qTip ), ret2.qTip ) );
   s.q_mag_err = ret_test.qTip.norm() - 1.0;
   s.psi_diff_mag = (q.PsiL - ret2.y_final.Psi).norm();
   s.mz_mag = ret2.y_final.Mz.norm();

   Eigen::Vector3d	q_err_axis = (Utility::quaternion_mult( Utility::quaternion_conj( ret_test.qTip ), ret2.qTip )).template bottomRows<3>();
   s.q_err_x = 2.0*q_err_axis( 0 );
   s.q_err_y = 2.0*q_err_axis( 1 );
   s.q_err_z = 2.0*q_err_axis( 2 );

   DiscontinuityList	discs = internal::get_discontinuity_list( C, q );
   TInterval::IntervalList ivals = internal::get_integration_intervals( C, q, discs, typename internal::traits<SType>::Options() );
   double normal_length = internal::get_normal_length( ivals );

   s.normal_density = static_cast<double>(n_normal_points) / normal_length;
   s.normal_length = normal_length;

   return s;
}

template <typename SType, typename STypeUnloaded>
stats get_stats_from_rets_with_loads( KinRet<SType> const& ret_test, KinRet<SType> const& ret2, KinRet<STypeUnloaded> const &retUnloaded, int n_normal_points, double dense_density, Cannula3 const& C, LoadedConfiguration3 const& q )
{
   stats s;

   s.n_normal_points = n_normal_points;
   s.dense_density = dense_density;
   s.p_err = (ret_test.pTip - ret2.pTip).norm();
   s.q_err = (ret_test.qTip - ret2.qTip).norm();
   s.mz_err = (ret_test.y_final.Mz - ret2.y_final.Mz).norm();
   s.psi_err = (ret_test.y_final.Psi - ret2.y_final.Psi).norm();

   s.J_va_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 0, 0 ) );
   s.J_vt_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 0, 3 ) );
   s.J_wa_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 3, 0 ) );
   s.J_wt_err = singular_value_error( (CTR::GetTipJacobianForTube1( ret_test.y_final ) - CTR::GetTipJacobianForTube1( ret2.y_final )).block<3, 3>( 3, 3 ) );

   s.beta1 = q.Beta( 0 );
   s.beta2 = q.Beta( 1 );
   s.beta3 = q.Beta( 2 );
   s.psi1 = q.PsiL( 0 );
   s.psi2 = q.PsiL( 1 );
   s.psi3 = q.PsiL( 2 );
   s.L1 = std::get<0>( C ).GetLength();
   s.L2 = std::get<1>( C ).GetLength();
   s.L3 = std::get<2>( C ).GetLength();
   s.Lt1 = std::get<0>( C ).GetTransmissionLength();
   s.Lt2 = std::get<1>( C ).GetTransmissionLength();
   s.Lt3 = std::get<2>( C ).GetTransmissionLength();

   using Mathematics::eval_pt;
   eval_pt pt1( (s.Lt1 + s.L1) / 2.0, s.Lt1, s.L1 );
   eval_pt pt2( (s.Lt2 + s.L2) / 2.0, s.Lt2, s.L2 );
   eval_pt pt3( (s.Lt3 + s.L3) / 2.0, s.Lt3, s.L3 );
   s.k1 = std::get<0>( C ).GetPrecurvature( pt1.s )(0);
   s.k2 = std::get<1>( C ).GetPrecurvature( pt2.s )(0);
   s.k3 = std::get<2>( C ).GetPrecurvature( pt3.s )(0);

   s.kb1 = std::get<0>( C ).GetBendingStiffness( pt1 );
   s.kb2 = std::get<1>( C ).GetBendingStiffness( pt2 );
   s.kb3 = std::get<2>( C ).GetBendingStiffness( pt3 );

   s.ct1 = std::get<0>( C ).GetTorsionalCompliance( pt1 );
   s.ct2 = std::get<1>( C ).GetTorsionalCompliance( pt2 );
   s.ct3 = std::get<2>( C ).GetTorsionalCompliance( pt3 );

   s.stability = CTR::GetStability( ret2.y_final );
   s.q_angle_err = Utility::quaternion_angle( Utility::quaternion_mult( Utility::quaternion_conj( ret_test.qTip ), ret2.qTip ) );
   s.q_mag_err = ret_test.qTip.norm() - 1.0;
   s.psi_diff_mag = (q.PsiL - ret2.y_final.Psi).norm();
   s.mz_mag = ret2.y_final.Mz.norm();

   Eigen::Vector3d	q_err_axis = (Utility::quaternion_mult( Utility::quaternion_conj( ret_test.qTip ), ret2.qTip )).template bottomRows<3>();
   s.q_err_x = 2.0*q_err_axis( 0 );
   s.q_err_y = 2.0*q_err_axis( 1 );
   s.q_err_z = 2.0*q_err_axis( 2 );

   DiscontinuityList	discs = internal::get_discontinuity_list( C, q );
   TInterval::IntervalList ivals = internal::get_integration_intervals( C, q, discs, typename internal::traits<SType>::Options() );
   double normal_length = internal::get_normal_length( ivals );

   s.normal_density = static_cast<double>(n_normal_points) / normal_length;
   s.normal_length = normal_length;

   s.Fx = q.Ftip( 0 );
   s.Fy = q.Ftip( 1 );
   s.Fz = q.Ftip( 2 );

   s.tip_deflection = (retUnloaded.pTip - ret2.pTip).norm();

   s.normal_points_used = ret_test.normal_points_used;
   return s;
}

template <typename SType>
void get_stats_from_rets_2( KinRet<SType> const& ret_test, KinRet<SType> const& ret2, int n_normal_points, double dense_density, Cannula3 const& C, Configuration3 const& q, stats_min_max_mean& stats )
{

   //s.n_normal_points = n_normal_points;
   //s.dense_density = dense_density;
   stats.p_err( (ret_test.pTip - ret2.pTip).norm() );
   stats.q_err( (ret_test.qTip - ret2.qTip).norm() );
   stats.psi_diff_mag( (ret_test.y_final.Psi - ret2.y_final.Psi).norm() );
   stats.mz_diff_mag( (ret_test.y_final.Mz - ret2.y_final.Mz).norm() );
   //s.q_err =
   //	s.mz_err = (ret_test.y_final.Mz - ret2.y_final.Mz).norm();
   //s.psi_err = (ret_test.y_final.Psi - ret2.y_final.Psi).norm();

   //s.J_va_err = singular_value_error((CTR::GetTipJacobianForTube1(ret_test.y_final) - CTR::GetTipJacobianForTube1(ret2.y_final)).block<3, 3>(0, 0));
   //s.J_vt_err = singular_value_error((CTR::GetTipJacobianForTube1(ret_test.y_final) - CTR::GetTipJacobianForTube1(ret2.y_final)).block<3, 3>(0, 3));
   //s.J_wa_err = singular_value_error((CTR::GetTipJacobianForTube1(ret_test.y_final) - CTR::GetTipJacobianForTube1(ret2.y_final)).block<3, 3>(3, 0));
   //s.J_wt_err = singular_value_error((CTR::GetTipJacobianForTube1(ret_test.y_final) - CTR::GetTipJacobianForTube1(ret2.y_final)).block<3, 3>(3, 3));

   //s.beta1 = q.Beta(0);
   //s.beta2 = q.Beta(1);
   //s.beta3 = q.Beta(2);
   //s.psi1 = q.PsiL(0);
   //s.psi2 = q.PsiL(1);
   //s.psi3 = q.PsiL(2);
   //s.L1 = std::get<0>(C).GetLength();
   //s.L2 = std::get<1>(C).GetLength();
   //s.L3 = std::get<2>(C).GetLength();
   //s.Lt1 = std::get<0>(C).GetTransmissionLength();
   //s.Lt2 = std::get<1>(C).GetTransmissionLength();
   //s.Lt3 = std::get<2>(C).GetTransmissionLength();

   //using Mathematics::eval_pt;
   //eval_pt pt1((s.Lt1 + s.L1) / 2.0, s.Lt1, s.L1);
   //eval_pt pt2((s.Lt2 + s.L2) / 2.0, s.Lt2, s.L2);
   //eval_pt pt3((s.Lt3 + s.L3) / 2.0, s.Lt3, s.L3);
   //s.k1 = std::get<0>(C).GetPrecurvature(pt1)(0);
   //s.k2 = std::get<1>(C).GetPrecurvature(pt2)(0);
   //s.k3 = std::get<2>(C).GetPrecurvature(pt3)(0);

   //s.kb1 = std::get<0>(C).GetBendingStiffness(pt1);
   //s.kb2 = std::get<1>(C).GetBendingStiffness(pt2);
   //s.kb3 = std::get<2>(C).GetBendingStiffness(pt3);

   //s.ct1 = std::get<0>(C).GetTorsionalCompliance(pt1);
   //s.ct2 = std::get<1>(C).GetTorsionalCompliance(pt2);
   //s.ct3 = std::get<2>(C).GetTorsionalCompliance(pt3);

   //s.stability = CTR::GetStability(ret2.y_final);
   //s.q_angle_err = Utility::quaternion_angle(Utility::quaternion_mult(Utility::quaternion_conj(ret_test.qTip), ret2.qTip));
   //s.q_mag_err = ret_test.qTip.norm() - 1.0;
   //s.psi_diff_mag = (q.PsiL - ret2.y_final.Psi).norm();
   //s.mz_mag = ret2.y_final.Mz.norm();

   //Eigen::Vector3d	q_err_axis = (Utility::quaternion_mult(Utility::quaternion_conj(ret_test.qTip), ret2.qTip)).bottomRows<3>();
   //s.q_err_x = 2.0*q_err_axis(0);
   //s.q_err_y = 2.0*q_err_axis(1);
   //s.q_err_z = 2.0*q_err_axis(2);

   //DiscontinuityList	discs = internal::get_discontinuity_list(C, q);
   //TInterval::IntervalList ivals = internal::get_integration_intervals(C, q, discs);
   //double normal_length = internal::get_normal_length(ivals);

   //s.normal_density = static_cast<double>(n_normal_points) / normal_length;
   //s.normal_length = normal_length;
}

//vector<stats>	get_stats_for_cannula( Cannula3 const& C, std::default_random_engine eng )
//{
//   int NConfigs = 1000;
//   vector<stats> vstats;
//
//   uniform_int_distribution<int> n_normal_points_dist( 0, 100 );
//   uniform_int_distribution<int> n_dense_dist( 3, 20 );
//
//   typedef DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options OType;
//   typedef State<3, OType> SType;
//
//   for ( int i = 0; i < NConfigs; ++i )
//   {
//      Configuration3 q = get_random_config_for_cannula( C, eng );
//      int n_normal_points = 0;// n_normal_points_dist(eng);
//      int dense_points = 0;// n_dense_dist(eng);
//
//      KinRet<SType>	ret_test = Kinematics( C, q, OType(), n_normal_points, static_cast<double>(dense_points) / 0.001 );
//      KinRet<SType>	ret_ref = Kinematics( C, q, OType(), 250, 30.0 / 0.001 );
//
//      stats ss = get_stats_from_rets( ret_test, ret_ref, n_normal_points, static_cast<double>(dense_points) / 0.001, C, q );
//      vstats.push_back( ss );
//   }
//
//
//   return vstats;
//}

//vector<stats>	get_stats_for_cannula_with_loads( Cannula3 const& C, std::default_random_engine eng )
//{
//   int NConfigs = 1000;
//   vector<stats> vstats;
//
//   uniform_int_distribution<int> n_normal_points_dist( 0, 100 );
//   uniform_int_distribution<int> n_dense_dist( 3, 20 );
//
//   typedef DeclareOptions< Option::ComputeGeometry, Option::ExternalLoads, Option::ComputeJacobian >::options OType;
//   typedef DeclareOptions< Option::ComputeGeometry >::options OTypeUnloaded;
//   typedef State<3, OType> SType;
//
//   for ( int i = 0; i < NConfigs; ++i )
//   {
//      LoadedConfiguration3 q = get_random_config_for_cannula_with_loads( C, eng );
//      int n_normal_points = 0;// n_normal_points_dist(eng);
//      int dense_points = 0;// n_dense_dist(eng);
//
//      KinRet<SType>	ret_test = Kinematics( C, q, OType(), n_normal_points, static_cast<double>(dense_points) / 0.001 );
//      KinRet<SType>	ret_ref = Kinematics( C, q, OType(), 100, 5.0 / 0.001 );
//      KinRet< State<3, OTypeUnloaded> > ret_unloaded = Kinematics( C, q, OTypeUnloaded(), 50, 5.0 / 0.001 );
//
//      stats ss = get_stats_from_rets_with_loads( ret_test, ret_ref, ret_unloaded, n_normal_points, static_cast<double>(dense_points) / 0.001, C, q );
//      vstats.push_back( ss );
//   }
//
//   return vstats;
//}

//void	get_stats_for_cannula_2( Cannula3 const& C, std::default_random_engine eng, stats_min_max_mean& stats )
//{
//   int NConfigs = 1000;
//
//   uniform_int_distribution<int> n_normal_points_dist( 0, 100 );
//   uniform_int_distribution<int> n_dense_dist( 3, 20 );
//
//   typedef DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options OType;
//   typedef State<3, OType> SType;
//
//   for ( int i = 0; i < NConfigs; ++i )
//   {
//      Configuration3 q = get_random_config_for_cannula( C, eng );
//      int n_normal_points = 0;// n_normal_points_dist(eng);
//      int dense_points = 0;// n_dense_dist(eng);
//
//      KinRet<SType>	ret_test = Kinematics( C, q, OType(), n_normal_points, static_cast<double>(dense_points) / 0.001 );
//      KinRet<SType>	ret_ref = Kinematics( C, q, OType(), 250, 30.0 / 0.001 );
//
//      get_stats_from_rets_2( ret_test, ret_ref, n_normal_points, static_cast<double>(dense_points) / 0.001, C, q, stats );
//   }
//}

//void test_random_convergence_2()
//{
//   const int NTrials = 1000000;
//
//   int count = 0;
//   int tid;
//   int numthreads;
//
////   int desiredthreads = omp_get_max_threads();
// //  if ( desiredthreads > 1 )
// //     --desiredthreads;
// //  omp_set_num_threads( desiredthreads ); //reserve one core
//
//   std::random_device dev;
//   stats_min_max_mean stats;
//
//   std::mt19937	eng;
//
//   const int DECIMATION = 10;
//   vector<double> max_stats_p;
//
//#pragma omp parallel for private(tid, eng) shared(numthreads,count,stats)
//   for ( int i = 0; i < NTrials; ++i )
//   {
//#pragma omp critical
//      {
//         eng.seed( dev() );
//      }
//      stats_min_max_mean private_stats;
//      Cannula3 CC = get_random_3tube_cannula( eng );
//      Configuration3 q = get_random_config_for_cannula( CC, eng );
//      int n_normal_points = 0;// n_normal_points_dist(eng);
//      int dense_points = 0;// n_dense_dist(eng);
//      typedef DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options OType;
//      typedef State<3, OType> SType;
//      KinRet<SType>	ret_test = Kinematics( CC, q, OType(), n_normal_points, static_cast<double>(dense_points) / 0.001 );
//      KinRet<SType>	ret_ref = Kinematics( CC, q, OType(), 250, 30.0 / 0.001 );
//      get_stats_from_rets_2( ret_test, ret_ref, n_normal_points, static_cast<double>(dense_points) / 0.001, CC, q, private_stats );
//
//#pragma omp critical
//      {
//         stats = stats + private_stats;
//      }
//
//#pragma omp atomic
//      count += 1;
////      tid = omp_get_thread_num();
//      if ( tid == 0 )
//      {
//         boost::io::ios_all_saver saver( std::cout );
////         numthreads = omp_get_num_threads();
//         std::cout << setw( 6 ) << setprecision( 2 ) << std::fixed << std::right << static_cast<double>(count) / static_cast<double>(NTrials)* 100 << "% complete" << std::endl;
//      }
//
//      if ( tid == 0 && count % DECIMATION == 0 )
//      {
//         max_stats_p.push_back( stats.p_err.max );
//      }
//   }
//
//   std::cout << "Number of threads: " << numthreads << std::endl;
//   std::cout << "count: " << count << std::endl;
//
//   std::cout << "Position statistics: " << std::endl;
//   std::cout << "\t mean: " << stats.p_err.mean << std::endl;
//   std::cout << "\t min:  " << stats.p_err.min << std::endl;
//   std::cout << "\t max:  " << stats.p_err.max << std::endl;
//
//   std::cout << "Quaternion statistics: " << std::endl;
//   std::cout << "\t mean: " << stats.q_err.mean << std::endl;
//   std::cout << "\t min:  " << stats.q_err.min << std::endl;
//   std::cout << "\t max:  " << stats.q_err.max << std::endl;
//
//   std::cout << "Psi statistics: " << std::endl;
//   std::cout << "\t mean: " << stats.psi_diff_mag.mean << std::endl;
//   std::cout << "\t min:  " << stats.psi_diff_mag.min << std::endl;
//   std::cout << "\t max:  " << stats.psi_diff_mag.max << std::endl;
//
//   std::cout << "Mz statistics: " << std::endl;
//   std::cout << "\t mean: " << stats.mz_diff_mag.mean << std::endl;
//   std::cout << "\t min:  " << stats.mz_diff_mag.min << std::endl;
//   std::cout << "\t max:  " << stats.mz_diff_mag.max << std::endl;
//
//   std::ofstream max_stats_p_file( "max_p_err.txt" );
//   std::copy( max_stats_p.begin(), max_stats_p.end(), std::ostream_iterator<double>( max_stats_p_file, "\n" ) );
//}

//void test_random_convergence()
//{
//   const int NDesigns = 10000;
//   vector< Cannula3, Eigen::aligned_allocator<Cannula3> >	Cannulae( NDesigns );
//   std::default_random_engine	eng( std::random_device{}() );
//   generate( Cannulae.begin(), Cannulae.end(), bind( get_random_3tube_cannula, std::ref( eng ) ) );
//
//   int count = 0;
//   int tid;
//   int numthreads;
//
////   int desiredthreads = omp_get_max_threads();
//   //if ( desiredthreads > 1 )
//  //    --desiredthreads;
// //  omp_set_num_threads( desiredthreads ); //reserve one core
//
//   vector<stats>	vstats;
//
//   std::random_device dev;
//
//#pragma omp parallel for private(tid, eng) shared(numthreads,count,vstats,Cannulae)
//   for ( int i = 0; i < NDesigns; ++i )
//   {
//#pragma omp critical
//      {
//         eng.seed( dev() );
//      }
//
//      vector<stats>	stats_for_design = get_stats_for_cannula( Cannulae[i], eng );
//
//#pragma omp critical
//      {
//         copy( stats_for_design.begin(), stats_for_design.end(), back_inserter( vstats ) );
//      }
//
//#pragma omp atomic
//      count += 1;
////      tid = omp_get_thread_num();
//      if ( tid == 0 )
//      {
//   //      numthreads = omp_get_num_threads();
//         std::cout << setw( 6 ) << setprecision( 2 ) << std::fixed << std::right << static_cast<double>(count) / static_cast<double>(NDesigns)* 100 << "% complete" << std::endl;
//      }
//   }
//
//   std::cout << "Number of threads: " << numthreads << std::endl;
//   std::cout << "count: " << count << std::endl;
//
//   ofstream data( "random_convergence_data.txt" );
//   std::copy( vstats.begin(), vstats.end(), ostream_iterator<stats>( data, "\n" ) );
//}

//void test_random_convergence_with_loads()
//{
//   const int NDesigns = 1000;
//   vector< Cannula3, Eigen::aligned_allocator<Cannula3> >	Cannulae( NDesigns );
//   std::default_random_engine	eng( std::random_device{}() );
//   generate( Cannulae.begin(), Cannulae.end(), bind( get_random_3tube_cannula, std::ref( eng ) ) );
//
//   int count = 0;
//   int tid;
//   int numthreads;
//
// //  int desiredthreads = omp_get_max_threads();
// //  if ( desiredthreads > 1 )
// //     --desiredthreads;
////   omp_set_num_threads( desiredthreads ); //reserve one core
//
//   vector<stats>	vstats;
//
//   std::random_device dev;
//
//#pragma omp parallel for private(tid, eng) shared(numthreads,count,vstats,Cannulae)
//   for ( int i = 0; i < NDesigns; ++i )
//   {
//#pragma omp critical
//      {
//         eng.seed( dev() );
//      }
//
//      vector<stats>	stats_for_design = get_stats_for_cannula_with_loads( Cannulae[i], eng );
//
//#pragma omp critical
//      {
//         copy( stats_for_design.begin(), stats_for_design.end(), back_inserter( vstats ) );
//      }
//
//#pragma omp atomic
//      count += 1;
//  //    tid = omp_get_thread_num();
//      if ( tid == 0 )
//      {
//  //       numthreads = omp_get_num_threads();
//         std::cout << setw( 6 ) << setprecision( 2 ) << std::fixed << std::right << static_cast<double>(count) / static_cast<double>(NDesigns)* 100 << "% complete" << std::endl;
//      }
//   }
//
//   std::cout << "Number of threads: " << numthreads << std::endl;
//   std::cout << "count: " << count << std::endl;
//
//   ofstream data( "random_convergence_data_with_loads.txt" );
//   std::copy( vstats.begin(), vstats.end(), ostream_iterator<stats>( data, "\n" ) );
//}

//void test_tmech_optimized()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           Kinematics Tests                *\n"
//      "*********************************************";
//   std::cout << msg << std::endl;
//   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;
//
//   constant_fun< Vector<2>::type > k1_fun( 17.4*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k2_fun( 6.8*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k3_fun( 3.5*Eigen::Vector2d::UnitX() );
//   T1_type T1 = make_annular_tube( 0.475, 0.4356, 1.17e-3, 0.76e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 0.296, 0.2414, 1.68e-3, 1.35e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 0.148, 0.0946, 2.32e-3, 1.87e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto cannula = std::make_tuple( T1, T2, T3 );
//   std::default_random_engine	eng( std::random_device{}() );
//
//   std::random_device  dev;
//   eng.seed( dev() ); //seed the rng with a real random number
//
//   int count = 0;
//   int tid;
//   int numthreads;
//
//   int NConfigs = 100000;
//   typedef CTR::DeclareOptions < 
//      CTR::Option::ComputeGeometry,
//      CTR::Option::ExternalLoads,
//      CTR::Option::ComputeJacobian 
//   >::options OType;
//
//   typedef CTR::DeclareOptions <
//      CTR::Option::ComputeGeometry
//   >::options OTypeUnloaded;
//
//   std::vector<stats>   vstats;
//
//   HWND console = GetConsoleWindow();
//   COORD pos;
//   pos.X = 0;
//   pos.Y = 4;
//
//#pragma omp parallel private(tid, eng) shared(numthreads,count,dev,vstats)
//   {
//#pragma omp critical
//      {
//         eng.seed( dev() );
//      }
// //     tid = omp_get_thread_num();
//      if ( tid == 0 )
// //        numthreads = omp_get_num_threads();
//#pragma omp for
//      for ( int i = 0; i < NConfigs; ++i )
//      {
//         LoadedConfiguration3 q = get_random_config_for_cannula_with_loads( cannula, eng );
//         auto ret = CTR::Kinematics( cannula, q, OType(), 10, 0.0 );
//         auto ret2 = CTR::Kinematics( cannula, q, OType(), 50, 5.0 / 0.001 );
//         auto ret3 = CTR::Kinematics( cannula, q, OTypeUnloaded(), 50, 5.0 / 0.001 );
//
//         stats s = get_stats_from_rets_with_loads( ret, ret2, ret3, 0, 0, cannula, q );
//
//#pragma omp critical
//         {
//            vstats.push_back( s );
//         }
//
//#pragma omp atomic
//         count += 1;
//
//         if ( tid == 0 )
//         {
//            SetConsoleCursorPosition( GetStdHandle(STD_OUTPUT_HANDLE), pos );
//            std::cout << setw( 6 ) << setprecision( 2 ) << std::fixed << std::right << static_cast<double>(count) / static_cast<double>(NConfigs)* 100 << "% complete" << std::endl;
//         }
//      }
//   }
//
//   std::cout << "Number of threads: " << numthreads << std::endl;
//   std::cout << "count: " << count << std::endl;
//
//   ofstream data( "tmech_design_convergence_data.txt" );
//   std::copy( vstats.begin(), vstats.end(), ostream_iterator<stats>( data, "\n" ) );
//
//}

//void test_difficult_case()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           Kinematics Tests                *\n"
//      "*********************************************";
//   std::cout << msg << std::endl;
//   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;
//
//   constant_fun< Vector<2>::type > k1_fun( 43.78965813966766*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k2_fun( 11.09403286485417*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k3_fun( 42.48967114625794*Eigen::Vector2d::UnitX() );
//   T1_type T1 = make_annular_tube( 0.348662246433066, 0.324559837621616, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 0.08062868572705952, 0.04577483549648491, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 0.06973303994407346, 0.04080058371891667, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto cannula = std::make_tuple( T1, T2, T3 );
//
//   LoadedConfiguration3 q;
//
//   q.Beta << -0.1168518631846467, -0.05175444135400894, -0.051420154317757;
//   q.PsiL << -0.483470616206056, -1.77829075500498, 1.336457625395942;
//   q.Ftip << 1.258503366261721, 2.919838404050097, 1.66758547257632;
//   //q.Ftip << 0.0, 0.0, 0.0;
//   q.Ttip << 0.0, 0.0, 0.0;
//
//   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian, Option::ExternalLoads >::options	OType;
//
//   KinRet<State<3, OType> > ret = CTR::Kinematics( cannula, q, OType(), 10, 0 );
//
//   std::cout << ret.y_final << std::endl;
//}


//void test_difficult_case_2()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           Kinematics Tests                *\n"
//      "*********************************************";
//   std::cout << msg << std::endl;
//   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;
//
//   constant_fun< Vector<2>::type > k1_fun( 17.4*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k2_fun( 6.8*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k3_fun( 3.5*Eigen::Vector2d::UnitX() );
//   T1_type T1 = make_annular_tube( 0.475, 0.4356, 1.17e-3, 0.76e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 0.296, 0.2414, 1.68e-3, 1.35e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 0.148, 0.0946, 2.32e-3, 1.87e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto cannula = std::make_tuple( T1, T2, T3 );
//   std::default_random_engine	eng( std::random_device{}() );
//
//   std::random_device  dev;
//   eng.seed( dev() ); //seed the rng with a real random number
//
//   int count = 0;
//   int tid;
//   int numthreads;
//
//   int NConfigs = 1000000;
//   typedef CTR::DeclareOptions <
//      CTR::Option::ComputeGeometry,
//      CTR::Option::ExternalLoads,
//      CTR::Option::ComputeJacobian
//   >::options OType;
//
//   LoadedConfiguration3 q;
//   q.Beta << -0.311350601073354, -0.263619970780423, -0.123007991162321;
//   q.PsiL << 2.488903552868916, 1.409117439717429, 1.353368654105309;
//   q.Ttip << 0.0, 0.0, 0.0;
//   //q.Ftip << 0.022573943249881, 0.021301753353328, 0.002035699319094;
//   q.Ftip << 0.0, 0.0, 0.0;
//
//   auto ret1 = CTR::Kinematics( cannula, q, OType() );
//   auto ret2 = CTR::Kinematics( cannula, q, OType(), 10, 0.0 );
//   
//   std::cout << "quaternion norm: " << ret1.qTip.norm() << std::endl;
//
//   std::cout << ret1.y_final << std::endl;
//   std::cout << ret2.y_final << std::endl;
//
//   std::cout << "Position error: " << (ret1.pTip - ret2.pTip).norm() << std::endl;
//}

//template <typename Configuration, typename Ret>
//void writeDenseDataToFile( std::ofstream& file, Configuration const &q, Ret const &ret )
//{
//   int N = static_cast<int>(ret.dense_state_output.size());
//   for ( int i = 0; i < N; ++i )
//   {
//      auto S = ret.dense_state_output[i];
//      file << ret.arc_length_points[i] << " ";
//      file << S.Psi.transpose() << " ";
//      file << S.Mz.transpose() << " ";
//      file << S.p.transpose() << " ";
//
//      Eigen::Matrix3d R = Utility::quaternion_to_rotation_matrix( S.q );
//      Eigen::Map<Eigen::Matrix<double, 9, 1> > m( R.data() );
//      file << m.transpose() << " ";
//      file << std::endl;
//   }
//}

//void make_bifurcation_movie_points()
//{
//   using namespace Functions;
//   boost::io::ios_all_saver iosaver( std::cout );
//   const char *msg = "*********************************************\n"
//      "*           Making Movie Points             *\n"
//      "*********************************************";
//   std::cout << msg << std::endl;
//   typedef Tube< constant_fun< Vector<2>::type> >	T1_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T2_type;
//   typedef Tube< constant_fun< Vector<2>::type> >	T3_type;
//
//   constant_fun< Vector<2>::type > k1_fun( 27.2969*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k2_fun( 27.2969*Eigen::Vector2d::UnitX() );
//   constant_fun< Vector<2>::type > k3_fun( 27.2969*Eigen::Vector2d::UnitX() );
//   T1_type T1 = make_annular_tube( 0.051, 0.0, 1.0e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T2_type T2 = make_annular_tube( 0.050, 0.0, 1.0e-3, 0.86e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
//   T3_type T3 = make_annular_tube( 0.050, 0.0, 1.0e-3, 0.86e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );
//
//   auto cannula = std::make_tuple( T1, T2, T3 );
//
//   typedef CTR::DeclareOptions< Option::ComputeGeometry >::options   OType;
//
//   std::vector<double> data( { 0.3142,    0.3142,
//                             0.3321,    0.3321,
//                             0.3501,    0.3501,
//                             0.3680,    0.3680,
//                             0.3860,    0.3860,
//                             0.4039,    0.4039,
//                             0.4219,    0.4219,
//                             0.4399,    0.4399,
//                             0.4578,    0.4578,
//                             0.4758,    0.4758,
//                             0.4937,    0.4937,
//                             0.5117,    0.5117,
//                             0.5296,    0.5296,
//                             0.5476,    0.5476,
//                             0.5656,    0.5656,
//                             0.5835,    0.5835,
//                             0.6015,    0.6015,
//                             0.6194,    0.6194,
//                             0.6374,    0.6374,
//                             0.6553,    0.6553,
//                             0.6733,    0.6733,
//                             0.6913,    0.6913,
//                             0.7092,    0.7092,
//                             0.7272,    0.7272,
//                             0.7451,    0.7451,
//                             0.7631,    0.7631,
//                             0.7810,    0.7810,
//                             0.7990,    0.7990,
//                             0.8170,    0.8170,
//                             0.8349,    0.8349,
//                             0.8529,    0.8529,
//                             0.8708,    0.8708,
//                             0.8888,    0.8888,
//                             0.9067,    0.9067,
//                             0.9247,    0.9247,
//                             0.9427,    0.9427,
//                             0.9606,    0.9606,
//                             0.9786,    0.9786,
//                             0.9965,    0.9965,
//                             1.0145,    1.0145,
//                             1.0324,    1.0324,
//                             1.0504,    1.0504,
//                             1.0683,    1.0683,
//                             1.0863,    1.0863,
//                             1.1043,    1.1043,
//                             1.1222,    1.1222,
//                             1.1402,    1.1402,
//                             1.1581,    1.1581,
//                             1.1761,    1.1761,
//                             1.1940,    1.1940,
//                             1.2120,    1.2120,
//                             1.2300,    1.2300,
//                             1.2479,    1.2479,
//                             1.2659,    1.2659,
//                             1.2838,    1.2838,
//                             1.3018,    1.3018,
//                             1.3197,    1.3197,
//                             1.3377,    1.3377,
//                             1.3557,    1.3557,
//                             1.3736,    1.3736,
//                             1.3916,    1.3916,
//                             1.4095,    1.4095,
//                             1.4275,    1.4275,
//                             1.4454,    1.4454,
//                             1.4634,    1.4634,
//                             1.4814,    1.4814,
//                             1.4993,    1.4993,
//                             1.5173,    1.5173,
//                             1.5352,    1.5352,
//                             1.5532,    1.5532,
//                             1.5711,    1.5711,
//                             1.5891,    1.5891,
//                             1.6071,    1.6071,
//                             1.6250,    1.6250,
//                             1.6430,    1.6430,
//                             1.6609,    1.6609,
//                             1.6789,    1.6789,
//                             1.6968,    1.6968,
//                             1.7148,    1.7148,
//                             1.7328,    1.7328,
//                             1.7507,    1.7507,
//                             1.7687,    1.7687,
//                             1.7866,    1.7866,
//                             1.8046,    1.8046,
//                             1.8225,    1.8225,
//                             1.8405,    1.8405,
//                             1.8585,    1.8585,
//                             1.8764,    1.8764,
//                             1.8944,    1.8944,
//                             1.9123,    1.9123,
//                             1.9303,    1.9303,
//                             1.9482,    1.9482,
//                             1.9662,    1.9662,
//                             1.9842,    1.9842,
//                             2.0021,    2.0021,
//                             2.0201,    2.0201,
//                             2.0380,    2.0380,
//                             2.0560,    2.0560,
//                             2.0739,    2.0739,
//                             2.0919,    2.0919,
//                             2.1099,    2.1099,
//                             2.1278,    2.1278,
//                             2.1458,    2.1458,
//                             2.1637,    2.1637,
//                             2.1817,    2.1817,
//                             2.1996,    2.1996,
//                             2.2176,    2.2176,
//                             2.2355,    2.2355,
//                             2.2535,    2.2535,
//                             2.2715,    2.2715,
//                             2.2894,    2.2894,
//                             2.3074,    2.3074,
//                             2.3253,    2.3253,
//                             2.3433,    2.3433,
//                             2.3612,    2.3612,
//                             2.3792,    2.3792,
//                             2.3972,    2.3972,
//                             2.4151,    2.4151,
//                             2.4331,    2.4331,
//                             2.4510,    2.4510,
//                             2.4690,    2.4690,
//                             2.4869,    2.4869,
//                             2.5049,    2.5049,
//                             2.5229,    2.5229,
//                             2.5408,    2.5408,
//                             2.5588,    2.5588,
//                             2.5767,    2.5767,
//                             2.5947,    2.5947,
//                             2.6126,    2.6126,
//                             2.6306,    2.6306,
//                             2.6486,    2.6486,
//                             2.6665,    2.6665,
//                             2.6845,    2.6845,
//                             2.7024,    2.7024,
//                             2.7204,    2.7204,
//                             2.7383,    2.7383,
//                             2.7563,    2.7563,
//                             2.7743,    2.7743,
//                             2.7922,    2.7922,
//                             2.8102,    2.8102,
//                             2.8281,    2.8281,
//                             2.8461,    2.8461,
//                             2.8640,    2.8640,
//                             2.8820,    2.8820,
//                             2.9000,    2.9000,
//                             2.9179,    2.9179,
//                             3.6757,    3.6757,
//                             3.6945,    3.6945,
//                             3.7134,    3.7134,
//                             3.7322,    3.7322,
//                             3.7511,    3.7511,
//                             3.7699,    3.7699,
//                             3.7888,    3.7888,
//                             3.8076,    3.8076,
//                             3.8265,    3.8265,
//                             3.8453,    3.8453,
//                             3.8642,    3.8642,
//                             3.8830,    3.8830,
//                             3.9019,    3.9019,
//                             3.9207,    3.9207,
//                             3.9396,    3.9396,
//                             3.9584,    3.9584,
//                             3.9773,    3.9773,
//                             3.9961,    3.9961,
//                             4.0150,    4.0150,
//                             4.0338,    4.0338,
//                             4.0527,    4.0527,
//                             4.0715,    4.0715,
//                             4.0904,    4.0904,
//                             4.1092,    4.1092,
//                             4.1281,    4.1281,
//                             4.1469,    4.1469,
//                             4.1658,    4.1658,
//                             4.1846,    4.1846,
//                             4.2035,    4.2035,
//                             4.2223,    4.2223,
//                             4.2412,    4.2412,
//                             4.2600,    4.2600,
//                             4.2788,    4.2788,
//                             4.2977,    4.2977,
//                             4.3165,    4.3165,
//                             4.3354,    4.3354,
//                             4.3542,    4.3542,
//                             4.3731,    4.3731,
//                             4.3919,    4.3919,
//                             4.4108,    4.4108,
//                             4.4296,    4.4296,
//                             4.4485,    4.4485,
//                             4.4673,    4.4673,
//                             4.4862,    4.4862,
//                             4.5050,    4.5050,
//                             4.5239,    4.5239,
//                             4.5427,    4.5427,
//                             4.5616,    4.5616,
//                             4.5804,    4.5804,
//                             4.5993,    4.5993,
//                             4.6181,    4.6181,
//                             4.6370,    4.6370,
//                             4.6558,    4.6558,
//                             4.6747,    4.6747,
//                             4.6935,    4.6935,
//                             4.7124,    4.7124,
//                             4.7312,    4.7312,
//                             4.7501,    4.7501,
//                             4.7689,    4.7689,
//                             4.7878,    4.7878,
//                             4.8066,    4.8066,
//                             4.8255,    4.8255,
//                             4.8443,    4.8443,
//                             4.8632,   4.8632,
//                             4.8820,    4.8820,
//                             4.9009,    4.9009,
//                             4.9197,    4.9197,
//                             4.9386,    4.9386,
//                             4.9574,    4.9574,
//                             4.9763,    4.9763,
//                             4.9951,    4.9951,
//                             5.0140,    5.0140,
//                             5.0328,    5.0328,
//                             5.0517,    5.0517,
//                             5.0705,    5.0705,
//                             5.0894,    5.0894,
//                             5.1082,    5.1082,
//                             5.1271,    5.1271,
//                             5.1459,    5.1459,
//                             5.1648,    5.1648,
//                             5.1836,    5.1836,
//                             5.2025,    5.2025,
//                             5.2213,    5.2213,
//                             5.2402,    5.2402,
//                             5.2590,    5.2590,
//                             5.2779,    5.2779,
//                             5.2967,    5.2967,
//                             5.3156,    5.3156,
//                             5.3344,    5.3344,
//                             5.3533,    5.3533,
//                             5.3721,    5.3721,
//                             5.3910,    5.3910,
//                             5.4098,    5.4098,
//                             5.4287,    5.4287,
//                             5.4475,    5.4475,
//                             5.4664,    5.4664,
//                             5.4852,    5.4852,
//                             5.5041,    5.5041,
//                             5.5229,    5.5229,
//                             5.5418,    5.5418,
//                             5.5606,    5.5606,
//                             5.5795,    5.5795,
//                             5.5983,    5.5983,
//                             5.6172,    5.6172,
//                             5.6360,    5.6360,
//                             5.6549,    5.6549,
//                             5.6737,    5.6737,
//                             5.6926,    5.6926,
//                             5.7114,    5.7114,
//                             5.7303,    5.7303,
//                             5.7491,    5.7491,
//                             5.7680,    5.7680,
//                             5.7868,    5.7868,
//                             5.8057,    5.8057,
//                             5.8245,    5.8245,
//                             5.8434,    5.8434,
//                             5.8622,    5.8622,
//                             5.8811,    5.8811,
//                             5.8999,    5.8999,
//                             5.9188,    5.9188,
//                             5.9376,    5.9376,
//                             5.9565,    5.9565,
//                             5.9753,    5.9753,
//                             5.9942,    5.9942,
//                             6.0130,    6.0130,
//                             6.0319,    6.0319,
//                             6.0507,    6.0507,
//                             6.0696,    6.0696,
//                             6.0884,    6.0884,
//                             6.1073,    6.1073,
//                             6.1261,    6.1261 } );
//
//   for ( int i = 0; i < data.size()/2; ++i )
//   {
//      int j = 2 * i;
//      double psiL1 = 0.0;
//      double psiL2 = data[j];
//      double psiL3 = data[j + 1];
//      Configuration3 qfirst;
//      qfirst.Beta << 0.0, 0.0, 0.0;
//      qfirst.PsiL << psiL1, psiL2, psiL3;
//
//      //Configuration3 qsecond;
//      //qsecond.Beta << 0.0, 0.0, 0.0;
//      //j = 2 * (i + 1);
//      //double psiL12 = 0.0;
//      //double psiL22 = data[j];
//      //double psiL32 = data[j + 1];
//      //qsecond.PsiL << psiL12, psiL22, psiL32;
////
////      auto ret = CTR::Kinematics_with_dense_output( cannula, qfirst, OType(), 50, 5.0 / 0.001 );
////      std::string number = std::to_string( i );
////      std::ofstream file3( "dense_output_data_" + number + ".txt" );
////      writeDenseDataToFile( file3, qfirst, ret );
////   }
////}

void test_simulink_crash()
{
   using Eigen::Matrix;
   using namespace Functions;
   constant_fun< Vector<2>::type > k1_fun( 20 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k2_fun( 10 * Eigen::Vector2d::UnitX() );
   constant_fun< Vector<2>::type > k3_fun( 5 * Eigen::Vector2d::UnitX() );
   auto T1 = make_annular_tube( 0.2, 0.15, 1.111e-3, 0.86e-3, k1_fun, 50e9, 50e9 / 2.0 / 1.33 );
   auto T2 = make_annular_tube( 0.1, 0.05, 1.4789e-3, 1.21e-3, k2_fun, 50e9, 50e9 / 2.0 / 1.33 );
   auto T3 = make_annular_tube( 0.05, 0.025, 1.8398e-3, 1.5789e-3, k3_fun, 50e9, 50e9 / 2.0 / 1.33 );

   auto cannula = std::make_tuple( T1, T2, T3 );
   Configuration3 q;

   q.PsiL << 1.0154, 1.6505, 2.6347;
   q.Beta << -0.1950, 0.0240, 0.0051;

   typedef CTR::DeclareOptions< Option::ComputeGeometry, Option::ComputeJacobian >::options	OType;
   KinRet< State<3, OType> >	ret;
   ret = CTR::Kinematics( cannula, q, OType(), 20, 0 );

   cout << ret.pTip << endl;
   cout << ret.qTip << endl;
}

class histogram_binner
{
public:
   histogram_binner( double min, double max, int nBins ) :
      m_bins(nBins, 0),
      m_min(min),
      m_max(max)
   {

   }

   void add_measurement( double val )
   {
      m_bins[bin_number( val )]++;
   }

protected:
   unsigned int bin_number( double val )
   {
      double chi = (val - m_min) / (m_max - m_min);
      int nBins = m_bins.size();
      int number = (int)floor( chi*nBins );
      if ( number < m_bins.size() && number >= 0 )
      {
         return (unsigned int)number;
      }
      else if ( number < 0 )
      {
         return 0;
      }
      else
      {
         return m_bins.size() - 1;
      }
   }


private:
   std::vector<unsigned int> m_bins;
   double m_min;
   double m_max;

   friend ostream& operator<<(ostream &os, histogram_binner const &binner);
};

ostream& operator<<(ostream & os, histogram_binner const &binner)
{
   std::copy( binner.m_bins.begin(), binner.m_bins.end(), ostream_iterator<unsigned int>( os, "\n" ) );
   return os;
}

void test_random_designs_and_configs()
{
   default_random_engine eng;

   const long int N = 100000;
   histogram_binner binner( -20, 0, 40 );

   typedef CTR::DeclareOptions< Option::ComputeGeometry >::options OType;

   vector<double> errs( N );

   for ( int i = 0; i < N; ++i )
   {
      auto cannula = get_random_3tube_cannula( eng );
      auto configuration = get_random_config_for_cannula( cannula, eng );

      auto ret1 = Kinematics( cannula, configuration, OType(), 0, 0 );
      auto ret_good = Kinematics( cannula, configuration, OType(), 50, 5.0 / 0.001 );

      errs[i] = log10( (ret1.pTip - ret_good.pTip).norm() );
      if ( errs[i] == -std::numeric_limits<double>::infinity() )
      {
         errs[i] = -std::numeric_limits<double>::max();
      }
      binner.add_measurement( errs[i] );
   }

   cout << binner;
   ofstream p_errs( "p_tip_errs.txt" );
   std::copy( errs.begin(), errs.end(), ostream_iterator<double>( p_errs, "\n" ) );
}

int main( int argc, char* argv[] )
{
#ifdef WIN32
   BOOL ok = SetPriorityClass( GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
   if (!ok) {
      return 0;
   }
   ok = SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
   if (!ok) {
      return 0;
   }
#endif
   {
      //do_test(test_smooth_step_function);
      //do_test(test_identity_zero_functions);
      //do_test(test_tube_class);
      //do_test(test_get_precurvature);
      //do_test(test_tagged_intervals);
      //do_test(test_rk4_integration);
      //do_test(test_tuple_transform);
      //do_test(test_kinematics);
      //do_test(test_kinematics_with_loads);

      //do_test(test_two_tube_kinematics);
      //do_test(test_options);
      //do_test(test_mollified_indicator_function);
      //do_test(test_function_intervals);
      //do_test(test_smooth_step_integration);
      //do_test(test_kinematics_dense_points_convergence);
      //do_test(test_random_convergence);
      //do_test(test_random_convergence_2);
      //do_test( test_random_convergence_with_loads );
      //do_test( test_difficult_case );
      //do_test( test_tmech_optimized );
      //do_test( test_difficult_case_2 );
      do_test( test_dense_output );

      //make_bifurcation_movie_points();

      //Tests 12/30/2015
      //do_test( test_Zta_with_loads );
      //do_test( test_Zlf_with_loads );
      //do_test( test_Zgf_with_loads );
      //do_test( test_Zga_with_loads );
      //do_test( test_Zta_with_loads );
      //do_test( test_Ztf_with_loads );
      //do_test( test_kinematics_geom_only );
      //do_test( test_kinematics_torsion_only );
      //do_test( test_kinematics_Zga );
      //do_test( test_kinematics_loads );
      //do_test( test_kinematics_jacobians_loads );

      //Tests 12/30/2015
      //do_test( test_kinematics_geom_only );

      //do_test( test_simulink_crash );

      //Tests 12/31/2015
      //do_test( test_random_designs_and_configs );
   }
   return 0;
}

