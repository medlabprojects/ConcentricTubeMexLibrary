#include "BasicFunctions.h"
#include "Tube.h"
#include "Kinematics.h"
#include <mex.h>

using namespace CTR;
using namespace CTR::Functions;

#define PROGRAM_NAME "ThreeTubeMex"

const char * ERR_MSG_NOT_ENOUGH_ARGS = PROGRAM_NAME" did not receive the correct number of arguments";

struct Configuration
{
   Eigen::Vector3d PsiL;
   Eigen::Vector3d Beta;
} q;

typedef constant_fun<Eigen::Vector2d>   CurvFun;
typedef std::tuple< Tube<CurvFun>, Tube<CurvFun>, Tube<CurvFun> > CannulaT;
typedef DeclareOptions< Option::ComputeJacobian, Option::ComputeGeometry >::options OType;

void printHelp()
{
   const char * helpMsg = "function [] = ThreeTubeMex(T1,T2,T3,psiL,Beta)   \n"
                          "   Inputs:\n"
                          "      T1:   Tube structure for tube 1 (see definition below)\n"
                          "      T2:   Tube structure for tube 2 (see definition below)\n"
                          "      T3:   Tube structure for tube 3 (see definition below)\n"
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
                          "      k:  Tube precurvature (const)\n"
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
   if ( mxGetNumberOfElements( config ) != 3 ) return false;
   return true;
}

bool checkArgs( int nrhs, const mxArray *prhs[] )
{
   if ( nrhs != 6 )
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
   return q;
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
      const char *fields[] = { "p_tip", "q_tip", "J_tip" ,"s", "p", "q" };
      const size_t dims[1] = { 1 };
      plhs[0] = mxCreateStructArray( 1, dims, 6, fields );
	  
      using namespace CTR;
      Tube<CurvFun> t1 = GetTubeFromInput( prhs[0] );
      Tube<CurvFun> t2 = GetTubeFromInput( prhs[1] );
      Tube<CurvFun> t3 = GetTubeFromInput( prhs[2] );

      CannulaT cannula = std::make_tuple( t1, t2, t3 );
      Configuration q = GetConfigFromInput( prhs[3], prhs[4] );
	  
	  if (nrhs == 6) {
        n_normal_points = GetNNormalPoints( prhs[4] );
		sEval = GetSEval( prhs[5] );
      }
	  
      auto ret1 = Kinematics_with_dense_output( cannula, q, OType() );
      mxArray *p_tip = mxCreateDoubleMatrix( 3, 1, mxREAL );
      mxArray *q_tip = mxCreateDoubleMatrix( 4, 1, mxREAL );
      mxArray *J_tip = mxCreateDoubleMatrix( 6, 6, mxREAL );
      int Npts = ret1.arc_length_points.size();
      mxArray *s =     mxCreateDoubleMatrix( Npts, 1, mxREAL );
      mxArray *p =     mxCreateDoubleMatrix( Npts, 3, mxREAL );
      mxArray *qq =     mxCreateDoubleMatrix( Npts, 4, mxREAL );

      double *pdata = mxGetPr( p_tip );
      double *qdata = mxGetPr( q_tip );
      double *Jdata = mxGetPr( J_tip );
      double *sdata = mxGetPr( s );
      double *pdensedata = mxGetPr( p );
      double *qdensedata = mxGetPr( qq );

      std::copy( ret1.pTip.data(), ret1.pTip.data() + 3, pdata );
      std::copy( ret1.qTip.data(), ret1.qTip.data() + 4, qdata );   
      
      Eigen::Matrix<double, 6, 6> J = GetTipJacobianForTube1( ret1.y_final );
      std::copy( J.data(), J.data() + 36, Jdata );

      std::copy( ret1.arc_length_points.begin(), ret1.arc_length_points.end(), sdata );
      for ( int i = 0; i < Npts; ++i )
      {
         Eigen::Vector3d pi = ret1.dense_state_output.at( i ).p;
         *(pdensedata + i) = pi[0];
         *(pdensedata + i + Npts) = pi[1];
         *(pdensedata + i + 2 * Npts) = pi[2];
      }

      for ( int i = 0; i < Npts; ++i )
      {
         Eigen::Vector4d qi = ret1.dense_state_output.at( i ).q;
         *(qdensedata + i) =              qi[0];
         *(qdensedata + i + Npts) =       qi[1];
         *(qdensedata + i + 2 * Npts) =   qi[2];
         *(qdensedata + i + 3 * Npts) =   qi[3];
      }

      mxSetField( plhs[0], 0, "p_tip", p_tip );
      mxSetField( plhs[0], 0, "q_tip", q_tip );
      mxSetField( plhs[0], 0, "J_tip", J_tip );
      mxSetField( plhs[0], 0, "s", s );
      mxSetField( plhs[0], 0, "p", p );
      mxSetField( plhs[0], 0, "q", qq );

   }
   else
   {
      printHelp();
   }
}
