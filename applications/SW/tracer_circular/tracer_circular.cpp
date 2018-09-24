/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"

#include "TransientSystem.hpp"
#include "LinearImplicitSystem.hpp"

#include "slepceps.h"
#include <slepcmfn.h>

using namespace femus;

bool assembly = true;

double dt = 1.; //= dx / maxWaveSpeed * 0.85;

double k_v = 0.0001;

double pi = acos ( -1. );
double k_h = 1 / ( 10 * pi );

const unsigned NumberOfLayers = 20;

bool twostage = false;

//const double hRest[10]={1,1,1,1,1,1,1,1,1,1};
//const double hRest[10]={2,2,2,2,2,2,2,2,2,2};
//const double hRest[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const double hRest[20] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

double InitalValueV ( const std::vector < double >& x ) {
  return 1. / 10.;
}

// double InitalValueV0 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 +hRest[0]*(NumberOfLayers-1);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV1 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-2);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV2 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-3);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV3 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-4);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV4 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-5);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV5 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-6);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV6 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);;
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-7);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV7 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-8);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV8 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-9);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV9 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-10);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV10 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-11);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV11 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-12);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV12 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-13);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV13 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-14);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV14 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-15);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV15 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-16);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV16 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-17);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV17 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-18);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV18 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-19);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }
// double InitalValueV19 ( const std::vector < double >& x ) {
//   double psi1 = (10.-x[0])*x[0]/(5*5);
//   double z = -10 + hRest[0]/2 + hRest[0]*(NumberOfLayers-20);
//   double d_psi2 = (2*z*(z+10)*(2*z+10))/(5*5*5*5);
//   //double d_psi2 = (10 - 2*x[0])/(5*5);
//   return psi1*d_psi2;
// }


double InitalValueH ( const std::vector < double >& x ) {
  return hRest[0];
}

double InitalValueT ( const std::vector < double >& x ) {
  double pi = acos ( -1. );
//   return 17.5 + 25/pi * atan(x[0]/100.);
  if ( x[0] < 5 ) return 5;
  else return 30;
//  return (- sin(pi*x[0]));
}


double InitalValueB ( const std::vector < double >& x ) {
  return 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
}


bool SetBoundaryCondition ( const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time ) {
  bool dirichlet = false;
  if ( !strcmp ( SolName, "HT" ) ) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
  if ( !strcmp ( SolName, "T" ) ) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
  return dirichlet;
}


void ETDRosenbrock ( MultiLevelProblem& ml_prob );

void RK4 ( MultiLevelProblem& ml_prob );


int main ( int argc, char** args ) {

  SlepcInitialize ( &argc, &args, PETSC_NULL, PETSC_NULL );

  // init Petsc-MPI communicator
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = static_cast<unsigned> ( floor ( pow ( 2.,/*11*/4 ) + 0.5 ) ); //Grid cell size = 3.90625 m
  nx += 3;

  double length = 10.; //2 * 1465700.;

  //mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh" );
  mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, 0, length, 0., 0., 0., 0., EDGE3, "seventh" );

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol ( &mlMsh );

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    mlSol.AddSolution ( name, DISCONTINOUS_POLYNOMIAL, ZERO, 2 );
    sprintf ( name, "v%d", i );
    mlSol.AddSolution ( name, LAGRANGE, FIRST, 2 );
    sprintf ( name, "T%d", i );
    mlSol.AddSolution ( name, DISCONTINOUS_POLYNOMIAL, ZERO, 2 );
    sprintf ( name, "HT%d", i );
    mlSol.AddSolution ( name, DISCONTINOUS_POLYNOMIAL, ZERO, 2 );
  }

  mlSol.AddSolution ( "b", DISCONTINOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.AddSolution ( "eta", DISCONTINOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.Initialize ( "All" );

//   mlSol.Initialize("v0",InitalValueV0);
//   mlSol.Initialize("v1",InitalValueV1);
//   mlSol.Initialize("v2",InitalValueV2);
//   mlSol.Initialize("v3",InitalValueV3);
//   mlSol.Initialize("v4",InitalValueV4);
//   mlSol.Initialize("v5",InitalValueV5);
//   mlSol.Initialize("v6",InitalValueV6);
//   mlSol.Initialize("v7",InitalValueV7);
//   mlSol.Initialize("v8",InitalValueV8);
//   mlSol.Initialize("v9",InitalValueV9);
//   mlSol.Initialize("v10",InitalValueV10);
//   mlSol.Initialize("v11",InitalValueV11);
//   mlSol.Initialize("v12",InitalValueV12);
//   mlSol.Initialize("v13",InitalValueV13);
//   mlSol.Initialize("v14",InitalValueV14);
//   mlSol.Initialize("v15",InitalValueV15);
//   mlSol.Initialize("v16",InitalValueV16);
//   mlSol.Initialize("v17",InitalValueV17);
//   mlSol.Initialize("v18",InitalValueV18);
//   mlSol.Initialize("v19",InitalValueV19);

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    mlSol.Initialize ( name, InitalValueH );
  }

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "T%d", i );
    mlSol.Initialize ( name, InitalValueT );
  }

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "v%d", i );
    mlSol.Initialize ( name, InitalValueV );
  }

  mlSol.Initialize ( "b", InitalValueB );

  mlSol.AttachSetBoundaryConditionFunction ( SetBoundaryCondition );
  mlSol.GenerateBdc ( "All" );

  MultiLevelProblem ml_prob ( &mlSol );

  // ******* Add FEM system to the MultiLevel problem *******
  //TransientLinearImplicitSystem& system = ml_prob.add_system < TransientLinearImplicitSystem > ( "SWhv" );
  TransientLinearImplicitSystem& system2 = ml_prob.add_system < TransientLinearImplicitSystem > ( "SWt" );
  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
//     sprintf ( name, "h%d", i );
//     system.AddSolutionToSystemPDE ( name );
//     sprintf ( name, "v%d", i );
//     system.AddSolutionToSystemPDE ( name );
    sprintf ( name, "HT%d", i );
    system2.AddSolutionToSystemPDE ( name );
  }
  //system.init();
  system2.init();

  mlSol.SetWriter ( VTK );
  std::vector<std::string> print_vars;
  print_vars.push_back ( "All" );
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, 0 );

  unsigned numberOfTimeSteps = 1000; //17h=1020 with dt=60, 17h=10200 with dt=6
  for ( unsigned i = 0; i < numberOfTimeSteps; i++ ) {
      
    assembly = true;  
    system2.CopySolutionToOldSolution();
//     dt = 60.;
//     ETD ( ml_prob );
//     dt = 60.;
    ETDRosenbrock ( ml_prob );
    //RK4 ( ml_prob );
    mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, ( i + 1 ) / 1 );
  }
  return 0;
}


void ETDRosenbrock ( MultiLevelProblem& ml_prob ) {

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ( "SWt" ); // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)
  
  NumericVector* EPS = pdeSys->_EPS; // pointer to the global residual vector object in pdeSys (level)

  NumericVector* RES2; 
  RES2 = NumericVector::build().release();
  RES2->init(*RES);
  
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh ( NLayers );
  //std::vector < unsigned > solPdeIndexh ( NLayers );

  std::vector < unsigned > solIndexv ( NLayers );
  //std::vector < unsigned > solPdeIndexv ( NLayers );

  std::vector < unsigned > solIndexHT ( NLayers );
  std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT ( NLayers );

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for ( unsigned i = 0; i < NLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    solIndexh[i] = mlSol->GetIndex ( name ); // get the position of "hi" in the sol object
    //solPdeIndexh[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "hi" in the pdeSys object

    sprintf ( name, "v%d", i );
    solIndexv[i] = mlSol->GetIndex ( name ); // get the position of "vi" in the sol object
    //solPdeIndexv[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "vi" in the pdeSys object

    sprintf ( name, "HT%d", i );
    solIndexHT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf ( name, "T%d", i );
    solIndexT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType ( solIndexh[0] ); // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType ( solIndexv[0] ); // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType ( solIndexHT[0] ); // get the finite element type for "Ti"

  if(assembly) KK->zero();
  RES->zero();

  MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
    }
    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW ( NLayers, -1.e6 );
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
  for ( unsigned i =  start; i <  end; i++ ) {

    vector < double > solhm ( NLayers );
    vector < double > solh ( NLayers ); // local coordinates
    vector < double > solhp ( NLayers );
    vector < double > solvm ( NLayers ); // local coordinates
    vector < double > solvp ( NLayers ); // local coordinates
    vector < adept::adouble > solHTm ( NLayers ); // local coordinates
    vector < adept::adouble > solHT ( NLayers ); // local coordinates
    vector < adept::adouble > solHTp ( NLayers ); // local coordinates

    vector < adept::adouble > solHTmm ( NLayers ); // local coordinates
    vector < adept::adouble > solHTpp ( NLayers ); // local coordinates

    //vector< adept::adouble > aResh ( NLayers );
    //vector< adept::adouble > aResv ( NLayers );
    vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = ( i == start ) ? 0 : 1;
    unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

    unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
    unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize ( NLayers );
    l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );
    //l2GMapColumn.resize ( ( 1 + bc1 + bc2 + bc3 + bc4) * NLayers );

    //std::fill ( aResh.begin(), aResh.end(), 0 ); //set aRes to zero
    std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for ( unsigned j = 0; j < NLayers; j++ ) {

      solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
      solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
      //l2GMapRow[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      //l2GMapColumn[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
      solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

      //l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 0, i );
      //l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 1, i );

      if ( i > start ) {
        solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
        solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

        //l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i - 1 );
        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

      }

      if ( i < end - 1 ) {
        solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
        solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

        //l2GMapColumn[ ( 4 + 2 * bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i + 1 );
        l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

    if(assembly) s.new_recording();

    vector < double > x ( 2 ); // local coordinates
    for ( unsigned j = 0; j < 2; j++ ) {
      unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
      x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
    }
    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    double hTot = 0.;
    for ( unsigned k = 0; k < NLayers; k++ ) {
      hTot += solh[k]/*.value()*/;
    }

    std::vector < double > hALE ( NLayers, 0. );

    hALE[0] = hRest[0] + ( hTot - b );
    for ( unsigned k = 1; k < NLayers; k++ ) {
      hALE[k] = hRest[k];
    }

    std::vector < double > w ( NLayers + 1, 0. );
    w[0] = 1.;

//     NEW w
//     for ( unsigned k = NLayers; k > 1; k-- ) {
//       w[k - 1] = w[k] - ( hALE[k - 1] - solh[k - 1]) / dt;
//       if(bc2){
//         w[k - 1] -=   0.5 * ( solh[k - 1] + solhp[k - 1] ) * solvp[k - 1] /dx;
//       }
//       else{
//         w[k - 1] -=   solh[k - 1] * 1 /dx;
//       }
//       if(bc1){
//         w[k - 1] +=   0.5 * ( solh[k - 1] + solhm[k - 1] ) * solvm[k - 1] /dx;
//       }
//       else{
//         w[k - 1] +=   solh[k - 1] * 1 /dx;
//       }
//       //std::cout<< w[k-1] << " ";
//     }
//     //std::cout<<std::endl;

//     OLD w
//     for(unsigned k = 1; k < NLayers; k++){
//       w[k] = w[k-1] + (  0.5 * ( solh[k-1].value() + solhp[k-1].value() ) * solvp[k-1].value()
// 		       - 0.5 * ( solh[k-1].value() + solhm[k-1].value() ) * solvm[k-1].value() )/dx
// 		    + ( hALE[k-1] - solh[k-1].value()) / dt;//TODO
// 		    //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<w[k-1]<<std::endl;
//     }

    std::vector < double > zMid ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      zMid[k] = -b + solh[k] / 2.;
      for ( unsigned i = k + 1; i < NLayers; i++ ) {
        zMid[k] += solh[i];
      }
    }

    std::vector < double > psi2 ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
      //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
    }

    for ( unsigned k = NLayers; k > 1; k-- ) {
      //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
      w[k - 1] = 1.;
      if ( maxW[k - 1] < w[k - 1] ) {
        maxW[k - 1] = w[k - 1];
      }
    }


    for ( unsigned k = 0; k < NLayers; k++ ) {

      //BEGIN FIRST ORDER
      if ( i > start ) {
        //aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; //second order
        if ( solvm[k] > 0 ) {
          aResHT[k] += solHTm[k].value() * solvm[k] / dx;
        }
        else {
          aResHT[k] += solHT[k].value() * solvm[k] / dx;
        }
      }
      if ( i < end - 1 ) {
        //aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; //second order
        if ( solvp[k] > 0 ) {
          aResHT[k] -= solHT[k].value() * solvp[k] / dx; //first order upwind
        }
        else {
          aResHT[k] -= solHTp[k].value() * solvp[k] / dx; //first order upwind
        }
      }
//       else{
//         aResHT[k] -= solHT[k] /*.value()*/  / dx; //first order upwind
//       }
      //END

      //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
      //END

      if ( k < NLayers - 1 ) {
        //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
        aResHT[k] += w[k + 1] * ( solHT[k + 1] / solh[k + 1] );
      }
      if ( k >= 0 ) {
        //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
        aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
      }

      adept::adouble deltaZt = 0.;
      adept::adouble deltaZb = 0.;
      adept::adouble ht = 0.;
      adept::adouble hb = 0.;
      if ( k > 0 ) {
        ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
        deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
        //aResv[k] -= 0.5 * w[k] * deltaZt;
      }
      else {
        ht = 0.5 * ( solhm[k] + solhp[k] );
        deltaZt = 0.* ( 0. - solHT[k] ) / ht;
      }
      if ( k < NLayers - 1 ) {
        hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
        deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
        //aResv[k] -= 0.5 * w[k+1] * deltaZb;
      }
      else {
        hb = 0.5 * ( solhm[k] + solhp[k] );
        deltaZb = 0.* ( solHT[k] - 0. ) / hb;
      }

      //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

//       aResHT[k] += solhm[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2. ); // vertical diffusion
//
//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

    }

    vector< double > Res ( NLayers ); // local redidual vector
    vector< double > solht ( NLayers ); // local redidual vector
    for ( unsigned k = 0; k < NLayers; k++ ) {
      //Res[k] =  aResh[k].value();
      Res[k] =  aResHT[k].value();
      solht[k] = solHT[k].value();
      //std::cout<< "Res["<<k<<"] = " << Res[k] <<std::endl;
      //std::cout<< "Res["<<NLayers+k<<"] = " << Res[NLayers+k] <<std::endl;
    }

    RES->add_vector_blocked ( Res, l2GMapRow );
    
    if(assembly){
      //s.dependent ( &aResh[0], NLayers );
      s.dependent ( &aResHT[0], NLayers );

      // define the independent variables
      //s.independent ( &solh[0], NLayers );
      s.independent ( &solHT[0], NLayers );
      //s.independent ( &solvm[0], NLayers );
      //s.independent ( &solvp[0], NLayers );
      if ( i > start ) {
        //s.independent ( &solhm[0], NLayers );
        s.independent ( &solHTm[0], NLayers );
      }
      if ( i < end - 1 ) {
        //s.independent ( &solhp[0], NLayers );
        s.independent ( &solHTp[0], NLayers );
      }
      /*    if ( i > start + 1) {
        s.independent ( &solHTmm[0], NLayers );
        }
        if ( i < end - 2 ) {
          s.independent ( &solHTpp[0], NLayers );
        } */

      // get the jacobian matrix (ordered by row major )
      vector < double > Jac ( NLayers * NLayers * ( 1 + bc1 + bc2 ) );
      //vector < double > Jac ( NLayers * NLayers * ( 1 + bc1 + bc2 + bc3 +bc4 ) );
      s.jacobian ( &Jac[0], true );

      //store K in the global matrix KK
      KK->add_matrix_blocked ( Jac, l2GMapRow, l2GMapColumn );

      s.clear_independents();
      s.clear_dependents();
    }
  }

  RES->close();
  if(assembly) KK->close();
 

  for ( unsigned k = 0; k < NLayers; k++ ) {
    std::cout << "layer " << k << " " << maxW[k] << std::endl;
  }

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;
// //

//  abort();
  MFN mfn;
  Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
  FN f, f1, f2, f3 , f4;

  //std::cout << "dt = " << dt << " dx = "<< dx << " maxWaveSpeed = "<<maxWaveSpeed << std::endl;
  std::cout << "dt = " << dt << std::endl;

  //dt = 100.;

  Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
  Vec y = ( static_cast< PetscVector* > ( EPS ) )->vec();

  MFNCreate ( PETSC_COMM_WORLD, &mfn );

  MFNSetOperator ( mfn, A );
  MFNGetFN ( mfn, &f );

  FNPhiSetIndex ( f, 1 );
  FNSetType ( f, FNPHI );
// FNView(f,PETSC_VIEWER_STDOUT_WORLD);

  FNSetScale ( f, dt, dt );
  MFNSetFromOptions ( mfn );

  MFNSolve ( mfn, v, y );
  MFNDestroy ( &mfn );

  sol->UpdateSol ( mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset );
    
  if ( twostage == true ) {

    RES2->zero();
   
    for ( unsigned i =  start; i <  end; i++ ) {

      vector < double > solhm ( NLayers );
      vector < double > solh ( NLayers ); // local coordinates
      vector < double > solhp ( NLayers );
      vector < double > solvm ( NLayers ); // local coordinates
      vector < double > solvp ( NLayers ); // local coordinates
      vector < adept::adouble > solHTm ( NLayers ); // local coordinates
      vector < adept::adouble > solHT ( NLayers ); // local coordinates
      vector < adept::adouble > solHTp ( NLayers ); // local coordinates

      vector < adept::adouble > solHTmm ( NLayers ); // local coordinates
      vector < adept::adouble > solHTpp ( NLayers ); // local coordinates

      //vector< adept::adouble > aResh ( NLayers );
      //vector< adept::adouble > aResv ( NLayers );
      vector< adept::adouble > aResHT ( NLayers );

      unsigned bc1 = ( i == start ) ? 0 : 1;
      unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

      unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
      unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

      l2GMapRow.resize ( NLayers );
      l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );
      //l2GMapColumn.resize ( ( 1 + bc1 + bc2 + bc3 + bc4) * NLayers );

      //std::fill ( aResh.begin(), aResh.end(), 0 ); //set aRes to zero
      std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

      for ( unsigned j = 0; j < NLayers; j++ ) {

        solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
        solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
        //l2GMapRow[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
        l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

        //l2GMapColumn[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
        l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

        solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
        solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

        //l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 0, i );
        //l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 1, i );

        if ( i > start ) {
          solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
          solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

          //l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i - 1 );
          l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

        }

        if ( i < end - 1 ) {
          solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
          solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

          //l2GMapColumn[ ( 4 + 2 * bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i + 1 );
          l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
        }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

      }

     // s.new_recording();

      vector < double > x ( 2 ); // local coordinates
      for ( unsigned j = 0; j < 2; j++ ) {
        unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
        x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
      }
      double dx = x[1] - x[0];

      double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

      double hTot = 0.;
      for ( unsigned k = 0; k < NLayers; k++ ) {
        hTot += solh[k]/*.value()*/;
      }

      std::vector < double > hALE ( NLayers, 0. );

      hALE[0] = hRest[0] + ( hTot - b );
      for ( unsigned k = 1; k < NLayers; k++ ) {
        hALE[k] = hRest[k];
      }

      std::vector < double > w ( NLayers + 1, 0. );
      w[0] = 1.;

//     NEW w
//     for ( unsigned k = NLayers; k > 1; k-- ) {
//       w[k - 1] = w[k] - ( hALE[k - 1] - solh[k - 1]) / dt;
//       if(bc2){
//         w[k - 1] -=   0.5 * ( solh[k - 1] + solhp[k - 1] ) * solvp[k - 1] /dx;
//       }
//       else{
//         w[k - 1] -=   solh[k - 1] * 1 /dx;
//       }
//       if(bc1){
//         w[k - 1] +=   0.5 * ( solh[k - 1] + solhm[k - 1] ) * solvm[k - 1] /dx;
//       }
//       else{
//         w[k - 1] +=   solh[k - 1] * 1 /dx;
//       }
//       //std::cout<< w[k-1] << " ";
//     }
//     //std::cout<<std::endl;

//     OLD w
//     for(unsigned k = 1; k < NLayers; k++){
//       w[k] = w[k-1] + (  0.5 * ( solh[k-1].value() + solhp[k-1].value() ) * solvp[k-1].value()
// 		       - 0.5 * ( solh[k-1].value() + solhm[k-1].value() ) * solvm[k-1].value() )/dx
// 		    + ( hALE[k-1] - solh[k-1].value()) / dt;//TODO
// 		    //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<w[k-1]<<std::endl;
//     }

      std::vector < double > zMid ( NLayers );
      for ( unsigned k = 0; k < NLayers; k++ ) {
        zMid[k] = -b + solh[k] / 2.;
        for ( unsigned i = k + 1; i < NLayers; i++ ) {
          zMid[k] += solh[i];
        }
      }

      std::vector < double > psi2 ( NLayers );
      for ( unsigned k = 0; k < NLayers; k++ ) {
        psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
        //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
      }

      for ( unsigned k = NLayers; k > 1; k-- ) {
        //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
        w[k - 1] = 1.;
        if ( maxW[k - 1] < w[k - 1] ) {
          maxW[k - 1] = w[k - 1];
        }
      }


      for ( unsigned k = 0; k < NLayers; k++ ) {

        //BEGIN FIRST ORDER
        if ( i > start ) {
          //aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; //second order
          if ( solvm[k] > 0 ) {
            aResHT[k] += solHTm[k].value() * solvm[k] / dx;
          }
          else {
            aResHT[k] += solHT[k].value() * solvm[k] / dx;
          }
        }
        if ( i < end - 1 ) {
          //aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; //second order
          if ( solvp[k] > 0 ) {
            aResHT[k] -= solHT[k].value() * solvp[k] / dx; //first order upwind
          }
          else {
            aResHT[k] -= solHTp[k].value() * solvp[k] / dx; //first order upwind
          }
        }
//       else{
//         aResHT[k] -= solHT[k] /*.value()*/  / dx; //first order upwind
//       }
        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if ( k < NLayers - 1 ) {
          //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
          aResHT[k] += w[k + 1] * ( solHT[k + 1] / solh[k + 1] );
        }
        if ( k >= 0 ) {
          //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
          aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
        }

        adept::adouble deltaZt = 0.;
        adept::adouble deltaZb = 0.;
        adept::adouble ht = 0.;
        adept::adouble hb = 0.;
        if ( k > 0 ) {
          ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
          //aResv[k] -= 0.5 * w[k] * deltaZt;
        }
        else {
          ht = 0.5 * ( solhm[k] + solhp[k] );
          deltaZt = 0.* ( 0. - solHT[k] ) / ht;
        }
        if ( k < NLayers - 1 ) {
          hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
          //aResv[k] -= 0.5 * w[k+1] * deltaZb;
        }
        else {
          hb = 0.5 * ( solhm[k] + solhp[k] );
          deltaZb = 0.* ( solHT[k] - 0. ) / hb;
        }

        //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

//       aResHT[k] += solhm[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2. ); // vertical diffusion
//
//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

      }

      vector< double > Res ( NLayers ); // local redidual vector
      for ( unsigned k = 0; k < NLayers; k++ ) {
        //Res[k] =  aResh[k].value();
        Res[k] =  aResHT[k].value();
        //std::cout<< "Res["<<k<<"] = " << Res[k] <<std::endl;
        //std::cout<< "Res["<<NLayers+k<<"] = " << Res[NLayers+k] <<std::endl;
      }

      RES2->add_vector_blocked ( Res, l2GMapRow );
      
    }
    
    //come Konstantin fa su Matlab: R2 = F(U2) - fn - An*(U2 - u);
    //noi dobbiamo fare: R2 = RES2 - RES - KK * EPS;

    RES2->scale(-1.);
    RES2->add(*RES);
    RES2->add_vector(*EPS,*KK);
    RES2->scale(-1.);
    
  }

  //PARAVIEW
  unsigned solIndexeta = mlSol->GetIndex ( "eta" );
  unsigned solIndexb = mlSol->GetIndex ( "b" );
  sol->_Sol[solIndexeta]->zero();
  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    sol->_Sol[solIndexeta]->add ( *sol->_Sol[solIndexh[k]] );
  }
  sol->_Sol[solIndexeta]->add ( -1, *sol->_Sol[solIndexb] );

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueT = valueHT / valueH;
      //if (i == 0) valueT = 0.;
      //if (i == msh->_dofOffset[solTypeHT][iproc + 1] - 1 ) valueT = 0.;


      sol->_Sol[solIndexT[k]]->set ( i, valueT );
    }

    sol->_Sol[solIndexT[k]]->close();

  }

  
  
   delete RES2;
}


void RK4 ( MultiLevelProblem& ml_prob ) {

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ( "SWt" ); // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  //SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  //NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)
  //NumericVector* EPS = pdeSys->_EPS; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh ( NLayers );
  //std::vector < unsigned > solPdeIndexh ( NLayers );

  std::vector < unsigned > solIndexv ( NLayers );
  //std::vector < unsigned > solPdeIndexv ( NLayers );

  std::vector < unsigned > solIndexHT ( NLayers );
  std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT ( NLayers );

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for ( unsigned i = 0; i < NLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    solIndexh[i] = mlSol->GetIndex ( name ); // get the position of "hi" in the sol object
    //solPdeIndexh[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "hi" in the pdeSys object

    sprintf ( name, "v%d", i );
    solIndexv[i] = mlSol->GetIndex ( name ); // get the position of "vi" in the sol object
    //solPdeIndexv[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "vi" in the pdeSys object

    sprintf ( name, "HT%d", i );
    solIndexHT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf ( name, "T%d", i );
    solIndexT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType ( solIndexh[0] ); // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType ( solIndexv[0] ); // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType ( solIndexHT[0] ); // get the finite element type for "Ti"

  //KK->zero();
  //RES->zero();

  //MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
    }
    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW ( NLayers, -1.e6 );
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
  for ( unsigned i =  start; i <  end; i++ ) {

    vector < double > solhm ( NLayers );
    vector < double > solh ( NLayers ); // local coordinates
    vector < double > solhp ( NLayers );
    vector < double > solvm ( NLayers ); // local coordinates
    vector < double > solvp ( NLayers ); // local coordinates
    vector < double > solHTm ( NLayers ); // local coordinates
    vector < double > solHT ( NLayers ); // local coordinates
    vector < double > solHTp ( NLayers ); // local coordinates

    vector < adept::adouble > solHTmm ( NLayers ); // local coordinates
    vector < adept::adouble > solHTpp ( NLayers ); // local coordinates

    //vector< adept::adouble > aResh ( NLayers );
    //vector< adept::adouble > aResv ( NLayers );
    vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = ( i == start ) ? 0 : 1;
    unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

    unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
    unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize ( NLayers );
    l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );
    //l2GMapColumn.resize ( ( 1 + bc1 + bc2 + bc3 + bc4) * NLayers );

    //std::fill ( aResh.begin(), aResh.end(), 0 ); //set aRes to zero
    std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for ( unsigned j = 0; j < NLayers; j++ ) {

      solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
      solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
      //l2GMapRow[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      //l2GMapColumn[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
      solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

      //l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 0, i );
      //l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 1, i );

      if ( i > start ) {
        solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
        solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

        //l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i - 1 );
        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

      }

      if ( i < end - 1 ) {
        solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
        solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

        //l2GMapColumn[ ( 4 + 2 * bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i + 1 );
        l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

 //   s.new_recording();

    vector < double > x ( 2 ); // local coordinates
    for ( unsigned j = 0; j < 2; j++ ) {
      unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
      x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
    }
    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    double hTot = 0.;
    for ( unsigned k = 0; k < NLayers; k++ ) {
      hTot += solh[k]/*.value()*/;
    }

    std::vector < double > hALE ( NLayers, 0. );

    hALE[0] = hRest[0] + ( hTot - b );
    for ( unsigned k = 1; k < NLayers; k++ ) {
      hALE[k] = hRest[k];
    }

    std::vector < double > w ( NLayers + 1, 0. );
    w[0] = 1.;

//     NEW w
//     for ( unsigned k = NLayers; k > 1; k-- ) {
//       w[k - 1] = w[k] - ( hALE[k - 1] - solh[k - 1]) / dt;
//       if(bc2){
//         w[k - 1] -=   0.5 * ( solh[k - 1] + solhp[k - 1] ) * solvp[k - 1] /dx;
//       }
//       else{
//         w[k - 1] -=   solh[k - 1] * 1 /dx;
//       }
//       if(bc1){
//         w[k - 1] +=   0.5 * ( solh[k - 1] + solhm[k - 1] ) * solvm[k - 1] /dx;
//       }
//       else{
//         w[k - 1] +=   solh[k - 1] * 1 /dx;
//       }
//       //std::cout<< w[k-1] << " ";
//     }
//     //std::cout<<std::endl;

//     OLD w
//     for(unsigned k = 1; k < NLayers; k++){
//       w[k] = w[k-1] + (  0.5 * ( solh[k-1].value() + solhp[k-1].value() ) * solvp[k-1].value()
// 		       - 0.5 * ( solh[k-1].value() + solhm[k-1].value() ) * solvm[k-1].value() )/dx
// 		    + ( hALE[k-1] - solh[k-1].value()) / dt;//TODO
// 		    //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<w[k-1]<<std::endl;
//     }

    std::vector < double > zMid ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      zMid[k] = -b + solh[k] / 2.;
      for ( unsigned i = k + 1; i < NLayers; i++ ) {
        zMid[k] += solh[i];
      }
    }

    std::vector < double > psi2 ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
      //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
    }

    for ( unsigned k = NLayers; k > 1; k-- ) {
      //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
      w[k - 1] = 1.;
      if ( maxW[k - 1] < w[k - 1] ) {
        maxW[k - 1] = w[k - 1];
      }
    }


    std::vector < double > k1_RK ( NLayers, 0. );
    std::vector < double > k2_RK ( NLayers, 0. );
    std::vector < double > k3_RK ( NLayers, 0. );
    std::vector < double > k4_RK ( NLayers, 0. );

    for ( unsigned RK_step = 0; RK_step < 4; RK_step++ ) {
      for ( unsigned k = 0; k < NLayers; k++ ) {
        double LHS = 0.;
        double addition = 0.;
        if ( RK_step == 1 ) {
          addition = k1_RK[k] * 0.5;
        }
        else if ( RK_step == 2 ) {
          addition = k2_RK[k] * 0.5;
        }

        else if ( RK_step == 3 ) {
          addition = k3_RK[k];
        }

        //BEGIN FIRST ORDER
        if ( i > start ) {
          //aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; //second order
          if ( solvm[k] > 0 ) {
            LHS += ( solHTm[k] + addition ) * solvm[k] / dx;
          }
          else {
            LHS += ( solHT[k] + addition ) * solvm[k] / dx;
          }
        }
        if ( i < end - 1 ) {
          //aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; //second order
          if ( solvp[k] > 0 ) {
            LHS -= ( solHT[k] + addition ) * solvp[k] / dx; //first order upwind
          }
          else {
            LHS -= ( solHTp[k] + addition ) * solvp[k] / dx; //first order upwind
          }
        }
//       else{
//         aResHT[k] -= solHT[k] /*.value()*/  / dx; //first order upwind
//       }
        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if ( k < NLayers - 1 ) {
          //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
          LHS += w[k + 1] * ( ( solHT[k + 1] + addition ) / solh[k + 1] );
        }
        if ( k >= 0 ) {
          //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
          LHS -= w[k] * ( ( solHT[k] + addition ) / solh[k] );
        }

        if ( RK_step == 0 ) {
          k1_RK[k] = LHS * dt;
        }
        else if ( RK_step == 1 ) {
          k2_RK[k] = LHS * dt;
        }
        else if ( RK_step == 2 ) {
          k3_RK[k] = LHS * dt;
        }
        else {
          k4_RK[k] = LHS * dt;
        }

//         adept::adouble deltaZt = 0.;
//         adept::adouble deltaZb = 0.;
//         adept::adouble ht = 0.;
//         adept::adouble hb = 0.;
//         if ( k > 0 ) {
//           ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
//           deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
//           //aResv[k] -= 0.5 * w[k] * deltaZt;
//         }
//         else {
//           ht = 0.5 * ( solhm[k] + solhp[k] );
//           deltaZt = 0.* ( 0. - solHT[k] ) / ht;
//         }
//         if ( k < NLayers - 1 ) {
//           hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
//           deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
//           //aResv[k] -= 0.5 * w[k+1] * deltaZb;
//         }
//         else {
//           hb = 0.5 * ( solhm[k] + solhp[k] );
//           deltaZb = 0.* ( solHT[k] - 0. ) / hb;
//         }

        //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

//       aResHT[k] += solhm[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2. ); // vertical diffusion
//
//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

      }

    }

    for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
      double valueHT = solHT[k] + 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] );
      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
      sol->_Sol[solIndexHT[k]]->close();
    }

  }

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueT = valueHT / valueH;

      sol->_Sol[solIndexT[k]]->set ( i, valueT );
    }

    sol->_Sol[solIndexT[k]]->close();

  }

  for ( unsigned k = 0; k < NLayers; k++ ) {
    std::cout << "layer " << k << " " << maxW[k] << std::endl;
  }

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;
// //

//  abort();

//   unsigned solIndexeta = mlSol->GetIndex ( "eta" );
//   unsigned solIndexb = mlSol->GetIndex ( "b" );
//   sol->_Sol[solIndexeta]->zero();
//   for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//     sol->_Sol[solIndexeta]->add ( *sol->_Sol[solIndexh[k]] );
//   }
//   sol->_Sol[solIndexeta]->add ( -1, *sol->_Sol[solIndexb] );


}


