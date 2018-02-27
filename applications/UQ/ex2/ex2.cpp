
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/sfem_assembly.hpp"

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}

void GetEigenPair(MultiLevelProblem& ml_prob, Mat &CCSLEPc, Mat &MMSLEPc);

unsigned numberOfUniformLevels = 2;
double sigma = 0.4;
double L = 0.4;

int main() {

  //BEGIN eigenvalue problem instances
  
  int argc = 0;
  char** argv = NULL;
  char** args = NULL; //TODO  not sure about this
  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  EPS eps;
  EPSType eigSolverType = EPSKRYLOVSCHUR;
  EPSProblemType problemType = EPS_HEP;
  EPSWhich eigSorting = EPS_LARGEST_MAGNITUDE;
  PetscInt numberOfEigPairs = 1;
  PetscInt dimWorkingSpace = 20 * numberOfEigPairs;
  PetscInt convergedSolns, numberOfIterations;
  PetscInt  maxIterations = 20;
  PetscReal  error = 1.e-10;
  PetscReal  tol = 1.e-10;
  PetscScalar  kr, ki;
  Mat CCSLEPc;
  Mat MMSLEPc;
  Vec  xr, xi;
  
  //END


  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfSelectiveLevels = 0;

  mlMsh.ReadCoarseMesh("../input/square.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND, 2);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("All");

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("UQ");
  system.AddSolutionToSystemPDE("u");


  // ******* System FEM Assembly *******
  system.SetAssembleFunction(AssembleUQSys);
  system.SetMaxNumberOfLinearIterations(5);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance(1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetMgSmoother(GMRES_SMOOTHER);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 4);


  //BEGIN solve the eigenvalue problem

  GetEigenPair(ml_prob, CCSLEPc, MMSLEPc);

  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
  CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, CCSLEPc, MMSLEPc);
  CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, problemType);
  CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, numberOfEigPairs, dimWorkingSpace, 600);
  CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, eigSorting);
  CHKERRQ(ierr);
  ierr = EPSSetType(eps, eigSolverType);
  CHKERRQ(ierr);
  ierr = EPSGetConverged(eps, &convergedSolns);
  CHKERRQ(ierr);


  ierr = EPSSolve(eps);
  CHKERRQ(ierr);
// 
//   ierr = EPSGetIterationNumber(eps, &numberOfIterations);
//   CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of iterations of the method: %D\n", numberOfIterations);
//   CHKERRQ(ierr);
//   ierr = EPSGetType(eps, &eigSolverType);
//   CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", eigSolverType);
//   CHKERRQ(ierr);
//   ierr = EPSGetDimensions(eps, &numberOfEigPairs, NULL, NULL);
//   CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", numberOfEigPairs);
//   CHKERRQ(ierr);
//   ierr = EPSGetTolerances(eps, &tol, &maxIterations);
//   CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD, " Stopping condition: tol=%.4g, maxit=%D\n", (double)tol, maxIterations);
//   CHKERRQ(ierr);
// 
// // Display solution and clean up
// 
// 
// // Get number of converged approximate eigenpairs
// 
//   ierr = EPSGetConverged(eps, &convergedSolns);
//   CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", convergedSolns);
//   CHKERRQ(ierr);
//   if(convergedSolns > 0) {
// 
// // Display eigenvalues and relative errors
// 
//     ierr = PetscPrintf(PETSC_COMM_WORLD, " k ||Ax-kx||/||kx||\n" " ----------------- ------------------\n");
//     CHKERRQ(ierr);
// 
//     for(unsigned i = 0; i < convergedSolns; i++) {
// 
// // Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)
// 
//       ierr = EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
//       CHKERRQ(ierr);
// 
// // Compute the relative error associated to each eigenpair
// 
//       ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
//       CHKERRQ(ierr);
// 
//     }
//   }
// 
//   ierr = EPSDestroy(&eps);
//   CHKERRQ(ierr);
//   ierr = SlepcFinalize();
//   CHKERRQ(ierr);

  //END


  //system.MGsolve();


  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  return 0;

} //end main

void GetEigenPair(MultiLevelProblem & ml_prob, Mat &CCSLEPc, Mat &MMSLEPc) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("UQ");   // pointer to the linear implicit system named "Poisson"

  unsigned level = numberOfUniformLevels - 1;

  double sigma2 = sigma * sigma;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             MM = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  vector <double> phi1;  // local test function
  vector <double> phi2;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight1; // gauss point weight
  double weight2; // gauss point weight

  phi1.reserve(maxSize);
  phi2.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve(maxSize);
  l2GMap2.reserve(maxSize);

  vector < double > MMlocal;
  MMlocal.reserve(maxSize * maxSize);

  vector < double > CClocal;
  CClocal.reserve(maxSize * maxSize);

  MM->zero(); // Set to zero all the entries of the Global Matrix

  int MM_size = msh->_dofOffset[solType][nprocs];
  int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];

  SparseMatrix *CC;
  CC = SparseMatrix::build().release();
  CC->init(MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size);
  CC->zero();

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom1 = msh->GetElementType(iel);
    unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap1.resize(nDof1);

    for(int k = 0; k < dim; k++) {
      x1[k].resize(nDofx1);
    }

    MMlocal.resize(nDof1 * nDof1);    //resize
    std::fill(MMlocal.begin(), MMlocal.end(), 0);    //set Jac to zero

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDof1; i++) {
      l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofx1; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x1[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    vector < double > *nullDoublePointer = NULL;

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x, *nullDoublePointer);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      vector < double > xg1(dim, 0.);

      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg1[k] += x1[k][i] * phi1[i];
        }
      }

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned i1 = 0; i1 < nDof1; i1++) {
          MMlocal[ i * nDof1 + i1 ] += phi1[i] * phi1[i1] * weight1;
        }

        for(int jel = msh->_elementOffset[iproc]; jel < msh->_elementOffset[iproc + 1]; jel++) {

          short unsigned ielGeom2 = msh->GetElementType(jel);
          unsigned nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
          unsigned nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs

          // resize local arrays
          l2GMap2.resize(nDof2);

          for(int k = 0; k < dim; k++) {
            x2[k].resize(nDofx2);
          }

          CClocal.resize(nDof1 * nDof2);    //resize
          std::fill(CClocal.begin(), CClocal.end(), 0);    //set Jac to zero

          // local storage of global mapping and solution
          for(unsigned j = 0; j < nDof2; j++) {
            l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);    // global to global mapping between solution node and pdeSys dof
          }

          // local storage of coordinates
          for(unsigned j = 0; j < nDofx2; j++) {
            unsigned xDof  = msh->GetSolutionDof(j, jel, xType);    // global to global mapping between coordinates node and coordinate dof
            for(unsigned k = 0; k < dim; k++) {
              x2[k][j] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
            }
          }

          for(unsigned jg = 0; jg < msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber(); jg++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2, phi2, phi_x, *nullDoublePointer);

            // evaluate the solution, the solution derivatives and the coordinates in the gauss point
            vector < double > xg2(dim, 0.);

            for(unsigned j = 0; j < nDof2; j++) {
              for(unsigned k = 0; k < dim; k++) {
                xg2[k] += x2[k][j] * phi2[j];
              }
            }
            double dist = 0.;
            for(unsigned k = 0; k < dim; k++) {
              dist += fabs(xg1[k] - xg2[k]);
            }
            double C = sigma2 * exp(-dist / L);
            for(unsigned j = 0; j < nDof2; j++) {
              CClocal[i * nDof2 + j] += weight1 * phi1[i] * C * phi2[j] * weight2;
            }
          }
          CC->add_matrix_blocked(CClocal, l2GMap1, l2GMap2);
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    MM->add_matrix_blocked(MMlocal, l2GMap1, l2GMap1);
  } //end element loop for each process

  MM->close();
  CC->close();


  CCSLEPc = (static_cast<PetscMatrix*>(CC))->mat();
  MMSLEPc = (static_cast<PetscMatrix*>(MM))->mat();
  
//   delete CC; //TODO with this delete it cannot read CCSLEPc anymore in the main

  // ***************** END ASSEMBLY *******************
}



