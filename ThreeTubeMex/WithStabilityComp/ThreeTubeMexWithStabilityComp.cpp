/***************************************
*  Copyright 2016 Vanderbilt University
*	Author: Hunter Gilbert
****************************************/
#include "BasicFunctions.h"
#include "Tube.h"
#include "Kinematics.h"
#include "SE3.h"
#include <mex.h>

using namespace CTR;
using namespace CTR::Functions;

#define PROGRAM_NAME "ThreeTubeMexWithStabilityComp"

const char * ERR_MSG_NOT_ENOUGH_ARGS = PROGRAM_NAME" did not receive the correct number of arguments";
const char * ERR_MSG_WRONG_NUMBER_OF_POINTS = "The number of extra points must be greater than 0";

struct Configuration
{
   Eigen::Vector3d PsiL;
   Eigen::Vector3d Beta;
   Eigen::Vector3d Ftip;
   Eigen::Vector3d Ttip;
} q;

typedef constant_fun<Eigen::Vector2d>   CurvFun;
typedef std::tuple< Tube<CurvFun>, Tube<CurvFun>, Tube<CurvFun> > CannulaT;
typedef DeclareOptions< Option::ComputeJacobian, Option::ComputeGeometry, Option::ComputeStability, Option::ComputeCompliance>::options OType;

void printHelp()
{
   const char * helpMsg = "function [] = ThreeTubeMexWithPsi(T1,T2,T3,psiL,Beta)   \n"
                          "   Inputs:\n"
                          "      T1:   Tube structure for tube 1 (see definition below)\n"
                          "      T2:   Tube structure for tube 2 (see definition below)\n"
                          "      T3:   Tube structure for tube 2 (see definition below)\n"
                          "      psiL: The distal tip angles\n"
                          "      Beta: The tube translations\n"
                          "      n:    (optional) number of sampled points\n"
                          "   Outputs: \n"
                          "      p_tip: [3x1] tip position\n"
                          "      q_tip: [4x1] quaternion of tip orientation\n"
                          "      J_tip: [6x6] Jacobian matrix of p_tip and q_tip with respect to psiL and Beta\n"
                          "      s:     [nx1] arclengths of sampled points\n"
                          "      p:     [nx3] positions of sampled points\n"
                          "      q:     [nx4] quaternion orientation of sampled points\n"
                          "      J:     [nx1 cells] jacobian of sampled points\n"
                          "      psi:   [nx3] psi (rotations) of sampled points\n"
                          "      stability:  [double] stability at tip\n"
                          "      C_tip: [6x6] compliance matrix at tip\n"
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

// make sure that the passed in struct contains all of the necessary fields
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
   if ( mxGetNumberOfElements( config ) != 3 ) return false;
   return true;
}

bool checkArgs( int nrhs, const mxArray *prhs[] )
{
   if ( nrhs < 5 || nrhs > 6 )
   {
      mexWarnMsgTxt( ERR_MSG_NOT_ENOUGH_ARGS );
      return false;
   }
   
   if ( !checkConfig( prhs[3] ) ) return false;
   if ( !checkConfig( prhs[4] ) ) return false;

   if ( !checkTube( prhs[0] ) ) return false;
   if ( !checkTube( prhs[1] ) ) return false;
   if ( !checkTube( prhs[2] ) ) return false;

   return true;
}

// parse the tube info
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
   q.PsiL << psiL_data[0], psiL_data[1], psiL_data[2];
   q.Beta << Beta_data[0], Beta_data[1], Beta_data[2];
   q.Ftip << 0.0, 0.0, 0.0; // assume wrench at tip is all 0
   q.Ttip << 0.0, 0.0, 0.0; 
   return q;
}

void copy_outputs_to_matlab( KinRetDense<State<3,OType> > const &ret1, mxArray *plhs[] )
{
	using SE3::transformation;
	using SE3::transformation_tie;
	using SE3::transformation_inverse;
    
    // create empty pointers for tip matrices
	mxArray *p_tip = mxCreateDoubleMatrix( 3, 1, mxREAL );
	mxArray *q_tip = mxCreateDoubleMatrix( 4, 1, mxREAL );
	mxArray *J_tip = mxCreateDoubleMatrix( 6, 6, mxREAL );
	mxArray *C_tip = mxCreateDoubleMatrix( 6, 6, mxREAL );
    mxArray *stability = mxCreateDoubleMatrix( 1, 1, mxREAL );
    
    // create empty pointers for data at each point on backbone
	int Npts = static_cast<int>(ret1.arc_length_points.size());
	mxArray *s =     mxCreateDoubleMatrix( Npts, 1, mxREAL );
	mxArray *p =     mxCreateDoubleMatrix( Npts, 3, mxREAL );
	mxArray *qq =    mxCreateDoubleMatrix( Npts, 4, mxREAL );
	mxArray *psi =   mxCreateDoubleMatrix( Npts, 3, mxREAL );

	
	double *pdata = mxGetPr( p_tip );
	double *qdata = mxGetPr( q_tip );
	double *Jdata = mxGetPr( J_tip );
    double *Cdata = mxGetPr( C_tip );
	double *sdata = mxGetPr( s );
	double *pdensedata = mxGetPr( p );
	double *qdensedata = mxGetPr( qq );
	double *psidensedata=mxGetPr( psi );
    double *stabilityData = mxGetPr( stability );
	
    // copy data into results
	std::copy( ret1.pTip.data(), ret1.pTip.data() + 3, pdata );
	std::copy( ret1.qTip.data(), ret1.qTip.data() + 4, qdata );   
	
	Eigen::Matrix<double, 6, 6> J = GetTipJacobianForTube1( ret1.y_final );
	std::copy( J.data(), J.data() + 36, Jdata );

    Eigen::Matrix<double, 6, 6> C = GetTipComplianceForTube1( ret1.y_final );
    std::copy( C.data(), C.data() + 36, Cdata );

    double S = GetStability(ret1.y_final);
    *stabilityData = S;

	std::copy( ret1.arc_length_points.begin(), ret1.arc_length_points.end(), sdata );
	transformation g_base = transformation(Eigen::Vector3d::UnitZ()*ret1.arc_length_points.back(), 
														Eigen::Vector4d::UnitX()) 
									* transformation_inverse( ret1.y_final.p, ret1.y_final.q );
	
	for ( int i = 0; i < Npts; ++i )
	{
		Eigen::Vector3d const &pi = ret1.dense_state_output.at( i ).p;
		Eigen::Vector4d const &qi = ret1.dense_state_output.at( i ).q;
		Eigen::VectorXd const &psI = ret1.dense_state_output.at( i ).Psi;
		
		*(psidensedata + i)			   = psI[0];
		*(psidensedata + i + Npts)	   = psI[1];
		*(psidensedata + i + 2 * Npts) = psI[2];
		
		transformation g_star = g_base * transformation(pi,qi);
	
		*(pdensedata + i) = g_star.p[0];
		*(pdensedata + i + Npts) = g_star.p[1];
		*(pdensedata + i + 2 * Npts) = g_star.p[2];

		*(qdensedata + i) =              g_star.q[0];
		*(qdensedata + i + Npts) =       g_star.q[1];
		*(qdensedata + i + 2 * Npts) =   g_star.q[2];
		*(qdensedata + i + 3 * Npts) =   g_star.q[3];
	}
    
    // organize struct output
	mxSetField( plhs[0], 0, "p_tip", p_tip );
	mxSetField( plhs[0], 0, "q_tip", q_tip );
	mxSetField( plhs[0], 0, "J_tip", J_tip );
	mxSetField( plhs[0], 0, "s", s );
	mxSetField( plhs[0], 0, "p", p );
	mxSetField( plhs[0], 0, "q", qq );
	mxSetField( plhs[0], 0, "psi", psi);
    mxSetField( plhs[0], 0, "C_tip", C_tip );
    mxSetField( plhs[0], 0, "stability", stability);
	
	//handle the Jacobian at each point states a little differently, as 
	// a struct array
	size_t dimsJJ[1] = {static_cast<size_t>(Npts)};
	mxArray *JJ = mxCreateCellArray(1, dimsJJ);
	Eigen::Matrix<double,6,6> Ad_beta = Utility::Adjoint_p_q(ret1.y_final.p, ret1.y_final.q);
	Eigen::Matrix<double,6,6>	const &Zga_beta = ret1.y_final.Zga;
	for (int i = 0; i < Npts; ++i) 
	{
		mxArray *Ji = mxCreateDoubleMatrix(6, 6, mxREAL);
		double *JiData = mxGetPr(Ji);	
		Eigen::Matrix<double, 6, 6> const &Zga_i = ret1.dense_state_output.at(i).Zga;
		Eigen::Vector3d p_inv;
		Eigen::Vector4d q_inv;
		transformation_tie(p_inv, q_inv) = transformation_inverse( ret1.dense_state_output.at(i).p, ret1.dense_state_output.at(i).q );
		Eigen::Matrix<double,6,6>	Ad_s = Utility::Adjoint_p_q( p_inv, q_inv );

		Eigen::Matrix<double,6,6> Jgi = Zga_i - Ad_s * Ad_beta * Zga_beta;
		
		std::copy( Jgi.data(), Jgi.data()+36, JiData);
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
      const char *fields[] = { "p_tip", "q_tip", "J_tip" ,"s", "p", "q", "J", "psi", "C_tip", "stability"};
      size_t dims[1] = { 1 };
      plhs[0] = mxCreateStructArray( 1, dims, 10, fields );

      using namespace CTR;
      Tube<CurvFun> t1 = GetTubeFromInput( prhs[0] );
      Tube<CurvFun> t2 = GetTubeFromInput( prhs[1] );
	  Tube<CurvFun> t3 = GetTubeFromInput( prhs[2] );

      CannulaT cannula = std::make_tuple( t1, t2, t3 );
      Configuration q = GetConfigFromInput( prhs[3], prhs[4] );
      
      // optional input for variable number of points to sample
      int n_normal_points = 20;
      if (nrhs == 6) {
        n_normal_points = GetNNormalPoints( prhs[5] );
      }
//       int n_normal_points = GetNNormalPoints( prhs[5] );

      KinRetDense<State<3,OType> > ret1 = Kinematics_with_dense_output( cannula, q, OType(), n_normal_points );

      copy_outputs_to_matlab(ret1, plhs);
   }
   else
   {
      printHelp();
   }
}

