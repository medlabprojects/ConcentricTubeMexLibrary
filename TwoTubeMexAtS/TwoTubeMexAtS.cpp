/***************************************
*  Copyright 2016 Vanderbilt University
*	Author: Hunter Gilbert
****************************************/
#include "BasicFunctions.h"
#include "Tube.h"
#include "Kinematics.h"
#include "SE3.h"
#include <mex.h>
#include "matrix.h"

using namespace CTR;
using namespace CTR::Functions;

#define PROGRAM_NAME "TwoTubeMexAtS"

const char * ERR_MSG_NOT_ENOUGH_ARGS = PROGRAM_NAME" did not receive the correct number of arguments";
const char * ERR_MSG_WRONG_NUMBER_OF_POINTS = "The number of extra points must be greater than 0";
const char * ERR_MSG_WRONG_S = "The evaluation arclength must be within the allowwable domain";

struct Configuration
{
   Eigen::Vector2d PsiL;
   Eigen::Vector2d Beta;
} q;

typedef constant_fun<Eigen::Vector2d>   CurvFun;
typedef std::tuple< Tube<CurvFun>, Tube<CurvFun> > CannulaT;
typedef DeclareOptions< Option::ComputeJacobian, Option::ComputeGeometry >::options OType;

void printHelp()
{
   const char * helpMsg = "THIS HELP MESSAGE HAS BEEN INCORRECTLY LEFT UNREVISED."

						  "function [] = ThreeTubeMex(T1,T2,psiL,Beta)   \n"
                          "   Inputs:\n"
                          "      T1:   Tube structure for tube 1 (see definition below)\n"
                          "      T2:   Tube structure for tube 2 (see definition below)\n"
                          "      psiL: The distal tip angles\n"
                          "      Beta: The tube translations\n"
                          "   Outputs: \n"
                          "      p_tip: The tip position\n"
                          "      q_tip: The tip orientation\n"
                          "      J_tip: The jacobian matrix of p_tip and q_tip with respect to psiL and Beta\n"
                          "\n"
                          "   Tube structure definition: \n"
                          "      OD: Outer diameter\n"
                          "      ID: Inner diameter\n"
                          "      k:  Tube precurvature (scalar)\n"
                          "      L:  Total tube length\n"
                          "      Lt: Transmission length\n"
                          "      E:  Young's Modulus\n"
                          "      G:  Shear Modulus\n";
   mexPrintf( helpMsg );
}

bool checkTube( const mxArray *tube )
{
   if ( mxGetClassID( tube ) != mxSTRUCT_CLASS ) return false;
   if ( mxGetNumberOfFields( tube ) != 7 ) return false; //tubes need 7 pieces of information (see help info above)

   if ( mxGetField( tube, 0, "OD" ) == 0 ) return false;
   if ( mxGetField( tube, 0, "ID" ) == 0 ) return false;
   if ( mxGetField( tube, 0, "k" ) == 0 )  return false;
   if ( mxGetField( tube, 0, "L" ) == 0 )  return false;
   if ( mxGetField( tube, 0, "Lt" ) == 0 ) return false;
   if ( mxGetField( tube, 0, "E" ) == 0 )  return false;
   if ( mxGetField( tube, 0, "G" ) == 0 )  return false;

   return true;
}

bool checkConfig( const mxArray *config )
{
   if ( mxGetClassID( config ) != mxDOUBLE_CLASS ) return false;
   if ( mxGetNumberOfElements( config ) != 2 ) return false;
   return true;
}
bool checkS( const mxArray *S0 )
{
	return true; // implement input validation later
}

bool checkArgs( int nrhs, const mxArray *prhs[] )
{
   //tube1, tube2, psi, beta, extra, S0
   if ( nrhs < 5 || nrhs > 6 )
   {
      mexWarnMsgTxt( ERR_MSG_NOT_ENOUGH_ARGS );
      return false;
   }
   
   //if psi is wrong
   if ( !checkConfig( prhs[2] ) ) return false;
   //if beta is wrong
   if ( !checkConfig( prhs[3] ) ) return false;
   //if tube 1 is bad
   if ( !checkTube( prhs[0] ) ) return false;
   //if tube 2 is bad
   if ( !checkTube( prhs[1] ) ) return false;
   //if s is outside domain
   if ( !checkS( prhs[4] ) ) return false;
   return true;
}

Tube<CurvFun> GetTubeFromInput(const mxArray *tube)
{
   double ODval, IDval, kval, Lval, Ltval, Eval, Gval;
   mxArray *OD = mxGetField( tube, 0, "OD" );
   ODval = *mxGetPr( OD );
   mxArray *ID = mxGetField( tube, 0, "ID" );
   IDval = *mxGetPr( ID );
   mxArray *L = mxGetField( tube, 0, "L" );
   Lval = *mxGetPr( L );
   mxArray *Lt = mxGetField( tube, 0, "Lt" );
   Ltval = *mxGetPr( Lt );
   mxArray *E = mxGetField( tube, 0, "E" );
   Eval = *mxGetPr( E );
   mxArray *G = mxGetField( tube, 0, "G" );
   Gval = *mxGetPr( G );
   mxArray *k = mxGetField( tube, 0, "k" );
   kval = *mxGetPr( k );

   CurvFun kfun( kval*Eigen::Vector2d::UnitX() );

   return make_annular_tube( Lval, Ltval, ODval, IDval, kfun, Eval, Gval );
}

Configuration GetConfigFromInput( const mxArray *psiL, const mxArray *Beta )
{
   const double *psiL_data = mxGetPr( psiL );
   const double *Beta_data = mxGetPr( Beta );

   Configuration q;
   q.PsiL << psiL_data[0], psiL_data[1];
   q.Beta << Beta_data[0], Beta_data[1];
   return q;
}

void copy_outputs_to_matlab( KinRetDense<State<2,OType> > const &ret1, mxArray *plhs[] )
{
	using SE3::transformation;
	using SE3::transformation_tie;
	using SE3::transformation_inverse;

	mxArray *p_tip = mxCreateDoubleMatrix( 3, 1, mxREAL );
	mxArray *q_tip = mxCreateDoubleMatrix( 4, 1, mxREAL );
	mxArray *J_tip = mxCreateDoubleMatrix( 6, 4, mxREAL );

    int Npts = static_cast<int>(ret1.arc_length_points.size());
	mxArray *s =     mxCreateDoubleMatrix( Npts, 1, mxREAL );
	mxArray *p =     mxCreateDoubleMatrix( Npts, 3, mxREAL );
	mxArray *qq =    mxCreateDoubleMatrix( Npts, 4, mxREAL );

	double *pdata = mxGetPr( p_tip );
	double *qdata = mxGetPr( q_tip );
	double *Jdata = mxGetPr( J_tip );
	double *sdata = mxGetPr( s );
	double *pdensedata = mxGetPr( p );
	double *qdensedata = mxGetPr( qq );

	std::copy( ret1.pTip.data(), ret1.pTip.data() + 3, pdata );
	std::copy( ret1.qTip.data(), ret1.qTip.data() + 4, qdata );   
	
	Eigen::Matrix<double, 6, 4> J = GetTipJacobianForTube1( ret1.y_final );
	std::copy( J.data(), J.data() + 24, Jdata );

	std::copy( ret1.arc_length_points.begin(), ret1.arc_length_points.end(), sdata );
	transformation g_base = transformation(Eigen::Vector3d::UnitZ()*ret1.arc_length_points.back(), 
														Eigen::Vector4d::UnitX()) 
									* transformation_inverse( ret1.y_final.p, ret1.y_final.q );
	
	for ( int i = 0; i < Npts; ++i )
	{
		Eigen::Vector3d const &pi = ret1.dense_state_output.at( i ).p;
		Eigen::Vector4d const &qi = ret1.dense_state_output.at( i ).q;
		
		transformation g_star = g_base * transformation(pi,qi);
	
		*(pdensedata + i) = g_star.p[0];
		*(pdensedata + i + Npts) = g_star.p[1];
		*(pdensedata + i + 2 * Npts) = g_star.p[2];

		*(qdensedata + i) =              g_star.q[0];
		*(qdensedata + i + Npts) =       g_star.q[1];
		*(qdensedata + i + 2 * Npts) =   g_star.q[2];
		*(qdensedata + i + 3 * Npts) =   g_star.q[3];
	}

	mxSetField( plhs[0], 0, "p_tip", p_tip );
	mxSetField( plhs[0], 0, "q_tip", q_tip );
	mxSetField( plhs[0], 0, "J_tip", J_tip );
	mxSetField( plhs[0], 0, "s", s );
	mxSetField( plhs[0], 0, "p", p );
	mxSetField( plhs[0], 0, "q", qq );

	//handle the Jacobian states a little differently, as 
	// a struct array
	const size_t dimsJJ[1] = { Npts };
	mxArray *JJ = mxCreateCellArray(1, dimsJJ);
	Eigen::Matrix<double,6,6> Ad_beta = Utility::Adjoint_p_q(ret1.y_final.p, ret1.y_final.q);
	Eigen::Matrix<double,6,4>	const &Zga_beta = ret1.y_final.Zga;
	for (int i = 0; i < Npts; ++i) 
	{
		mxArray *Ji = mxCreateDoubleMatrix(6, 4, mxREAL);
		double *JiData = mxGetPr(Ji);	
		Eigen::Matrix<double, 6, 4> const &Zga_i = ret1.dense_state_output.at(i).Zga;
		Eigen::Vector3d p_inv;
		Eigen::Vector4d q_inv;
		transformation_tie(p_inv, q_inv) = transformation_inverse( ret1.dense_state_output.at(i).p, ret1.dense_state_output.at(i).q );
		Eigen::Matrix<double,6,6>	Ad_s = Utility::Adjoint_p_q( p_inv, q_inv );

		Eigen::Matrix<double,6,4> Jgi = Zga_i - Ad_s * Ad_beta * Zga_beta;
		
		std::copy( Jgi.data(), Jgi.data()+24, JiData);
		mxSetCell(JJ, i, Ji);	
	}
	mxSetField( plhs[0], 0, "J", JJ );
}

int GetNNormalPoints( const mxArray *pNum )
{
    double *d = mxGetPr(pNum);
    int n = static_cast<int>(*d);
    if (n < 0) {
        mexWarnMsgTxt( ERR_MSG_WRONG_NUMBER_OF_POINTS );
        n = 0;
    }
    return n;
}

double * GetSEval( const mxArray *pNum )
{
	double *d = mxGetPr(pNum);
	return d;
}

extern "C"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
   if ( !checkArgs(nrhs,prhs) )
   {
      printHelp();
      return;
   }
 
   if ( nlhs == 1 )
   {
      const char *fields[] = { "p_tip", "q_tip", "J_tip" ,"s", "p", "q", "J" };
	  //const char *fields[] = { "p_sEval", "q_sEval" };
      const size_t dims[1] = { 1 };
      plhs[0] = mxCreateStructArray( 1, dims, 7, fields );

      using namespace CTR;
      Tube<CurvFun> t1 = GetTubeFromInput( prhs[0] );
      Tube<CurvFun> t2 = GetTubeFromInput( prhs[1] );

      CannulaT cannula = std::make_tuple( t1, t2 );
      Configuration q = GetConfigFromInput( prhs[2], prhs[3] );
      
      int n_normal_points = 10;
	  double * sEval = nullptr;
      if (nrhs == 6) {
        n_normal_points = GetNNormalPoints( prhs[4] );
		sEval = GetSEval( prhs[5] );
      }
	  
	  //mexPrintf("sEval = %f\n", *sEval);
	  
      KinRetDense<State<2,OType> > ret1 = Kinematics_with_dense_output( cannula, q, OType(), n_normal_points, 0, *sEval );
	  //KinRetDense<State<2,OType> > ret1 = Kinematics_with_dense_output( cannula, q, OType(), 0, 0, *sEval );
	  
      copy_outputs_to_matlab(ret1, plhs);
   }
   else
   {
      printHelp();
   }
}

