#ifndef __femus_utils_FE_convergence_hpp__
#define __femus_utils_FE_convergence_hpp__


#include "MultiLevelSolution.hpp"
#include "Math.hpp"
#include "Assemble_unknown.hpp"




namespace femus {

class MultiLevelProblem;
class MultiLevelSolution;
class MultiLevelMesh;



class Solution_generation_single_level {
    
public: 

// What this function does is to yield a MultiLevelSolution object at each level, to construct a level hierarchy to be used for error analysis 
virtual const MultiLevelSolution  run_on_single_level(MultiLevelProblem & ml_prob,
                                                      MultiLevelMesh & ml_mesh,
                                                      const unsigned lev,
                                                      const std::vector< Unknown > & unknowns,
                                                      const std::vector< Math::Function< double > * > &  exact_sol,
                                                      const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                      const MultiLevelSolution::BoundaryFuncMLProb  SetBoundaryCondition,
                                                  const bool equation_solve
                                                     ) const = 0;
                                                     
};


/// The idea of this class is to generate a certain amount of MultilevelSolutions and to compute the convergence rate.
/// It may be used either to study Approximation Theory or to study Convergence of FE Solution of PDEs
template < class type = double >
class FE_convergence {
 
    
    
public: 
 

  void  convergence_study(MultiLevelProblem & ml_prob,
                          MultiLevelMesh & ml_mesh,
                          MultiLevelMesh & ml_mesh_all_levels,   //auxiliary
                          const unsigned max_number_of_meshes,   //auxiliary
                          const unsigned norm_flag,
                          const unsigned conv_rate_computation_method,
                          const unsigned volume_or_boundary,
                          const bool equation_solve,                                          // true only if I have a System inside
                          const Solution_generation_single_level & main_in,
                          const std::vector< Unknown > & unknowns,                            //vector of Solutions
                          const std::vector< Math::Function< double > * >  & exact_sol,       // analytical function to Approximate
                          const MultiLevelSolution::InitFuncMLProb SetInitialCondition,       // routine that does the Approximation - needed in all cases
                          const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition  // routine for the BC - needed only if I have a System inside
                         );    //vector of Exact Solutions, if available

    
private: 
    
static   std::vector < std::vector < std::vector < std::vector < type > > > >  initialize_vector_of_norms(const unsigned volume_or_boundary,
                                                                                                          const unsigned unknowns_size, 
                                                                                                          const unsigned max_number_of_meshes,                                                                                              const unsigned norm_flag);

    
   
static   const MultiLevelSolution  initialize_convergence_study(MultiLevelProblem & ml_prob,
                                                                const std::vector< Unknown > &  unknowns,  
                                                                const std::vector< Math::Function< double > * >  & exact_sol,
                                                                MultiLevelMesh & ml_mesh_all_levels, 
                                                                const unsigned max_number_of_meshes, 
                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition,
                                                                const bool equation_solve
                                                               );  


//   print the error and the order of convergence between different levels
static  void output_convergence_order(const std::vector < std::vector < std::vector < type > > > &  norm,
                                    const unsigned int u,
                                    const unsigned int i,
                                    const unsigned int n);



static  void output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < std::vector < std::vector < type > > > > &  norms, 
                                        const unsigned norm_flag, 
                                        const unsigned volume_or_boundary);


static  void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes , int loop_i); 


static  std::vector< type > compute_error_norms_volume_or_boundary(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                const std::vector<Gauss> & quad_rules,
                                                const MultiLevelSolution* ml_sol, 
                                                const MultiLevelSolution* ml_sol_all_levels,
                                                const std::string & unknown,
                                                const unsigned current_level,
                                                const unsigned norm_flag,
                                                const unsigned conv_rate_computation_method,
                                                const unsigned volume_or_boundary,
                                                const Math::Function< type > * ex_sol_in
                                             );


static  void compute_error_norms_per_unknown_per_level(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                       const std::vector<Gauss> & quad_rules,
                                                       const MultiLevelSolution* ml_sol_single_level,
                                                       MultiLevelSolution* ml_sol_all_levels,
                                                       const std::vector< Unknown > &  unknowns,
                                                       const unsigned i,
                                                       const unsigned norm_flag,
                                                       std::vector < std::vector < std::vector < std::vector < type > > > > &  norms,
                                                       const unsigned conv_rate_computation_method,
                                                       const unsigned volume_or_boundary,
                                                       const  std::vector< Math::Function< type > * > & ex_sol_in);


 
    
};






    

template < class type>
  void  FE_convergence< type >::convergence_study(MultiLevelProblem & ml_prob,
                                                  MultiLevelMesh & ml_mesh,
                                                  MultiLevelMesh & ml_mesh_all_levels,
                                                  const unsigned max_number_of_meshes,
                                                  const unsigned norm_flag,
                                                  const unsigned conv_rate_computation_method,
                                                  const unsigned volume_or_boundary,
                                                  const bool equation_solve,
                                                  const Solution_generation_single_level & main_in,
                                                  const std::vector< Unknown > & unknowns,
                                                  const std::vector< Math::Function< double > * >  & exact_sol,
                                                  const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                  const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition
                                                 ) {
   
      // Controls - BEGIN 
      if ( equation_solve == true &&  SetBoundaryCondition == NULL)  { std::cout << "Must provide a BC function" << std::endl; abort(); }
      
       if (conv_rate_computation_method == 1 && exact_sol.size() == 0)  { std::cout << "Must provide exact sol" << std::endl; abort(); }
       // Controls - END


  
    
    std::vector < std::vector < std::vector < std::vector < double > > > > norms = FE_convergence::initialize_vector_of_norms (  volume_or_boundary,
                                                                                                                                 unknowns.size(), 
                                                                                                                                 max_number_of_meshes,                                                                                                  norm_flag);
    
     MultiLevelSolution         ml_sol_all_levels = FE_convergence::initialize_convergence_study(ml_prob,
                                                                                                 unknowns,
                                                                                                 exact_sol,
                                                                                                 ml_mesh_all_levels, 
                                                                                                 max_number_of_meshes, 
                                                                                                 SetInitialCondition, 
                                                                                                 SetBoundaryCondition,
                                                                                                 equation_solve);
    
  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
            
       for (int lev = 0; lev < max_number_of_meshes; lev++) {
                  
            const MultiLevelSolution ml_sol_single_level = main_in.run_on_single_level(ml_prob,
                                                                                       ml_mesh,
                                                                                       lev,
                                                                                       unknowns,
                                                                                       exact_sol,
                                                                                       SetInitialCondition, 
                                                                                       SetBoundaryCondition,
                                                                                       equation_solve
                                                                                      );
            

            FE_convergence::compute_error_norms_per_unknown_per_level ( elem_all,
                                                                        ml_prob.GetQuadratureRuleAllGeomElems(),
                                                                        & ml_sol_single_level,
                                                                        & ml_sol_all_levels,
                                                                        unknowns,
                                                                        lev,
                                                                        norm_flag,
                                                                        norms,
                                                                        conv_rate_computation_method,
                                                                        volume_or_boundary,
                                                                        exact_sol);
        
      }
   
       FE_convergence::output_convergence_order_all(unknowns, norms, norm_flag, volume_or_boundary);
   
}

    
    
template < class type>
/*static*/   std::vector < std::vector < std::vector < std::vector < type > > > >  FE_convergence< type >::initialize_vector_of_norms(const unsigned volume_or_boundary,
                                                                                                                                      const unsigned unknowns_size, 
                                                                                                                                      const unsigned max_number_of_meshes,                                                                                              const unsigned norm_flag) {
   
       // VB, how many Unknowns, how many mesh levels, how many norms
       
   std::vector < std::vector < std::vector < std::vector < type > > > > norms( volume_or_boundary + 1 );
  
    for (unsigned int vb = 0; vb < norms.size(); vb++) { //0: volume, 1: boundary
          norms[vb].resize( unknowns_size );     
      for (unsigned int u = 0; u < norms[vb].size(); u++) {
              norms[vb][u].resize( max_number_of_meshes );
       for (int i = 0; i < norms[vb][u].size(); i++) {   // loop on the mesh level
               norms[vb][u][i].resize(norm_flag + 1);
               std::fill(norms[vb][u][i].begin(), norms[vb][u][i].end(), 0.);
           }   
       }
    }
    
    return norms; 
       
}

    
   
template < class type>
/*static*/   const MultiLevelSolution  FE_convergence< type >::initialize_convergence_study(MultiLevelProblem & ml_prob, 
                                                                                            const std::vector< Unknown > &  unknowns,
                                                                                            const std::vector< Math::Function< double > * > &  exact_sol,
                                                                                            MultiLevelMesh & ml_mesh_all_levels,
                                                                                            const unsigned max_number_of_meshes,
                                                                                            const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                                            const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                                                                            const bool equation_solve
                                                                                           )  {

 //Mesh: construct all levels - BEGIN  ==================
        unsigned numberOfUniformLevels_finest = max_number_of_meshes;
        ml_mesh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
 //Mesh: construct all levels - END ==================

 
 //Solution - BEGIN ==================
//         std::vector < MultiLevelSolution * >   ml_sol_all_levels(unknowns.size());
//                ml_sol_all_levels[u] = new MultiLevelSolution (& ml_mesh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes
               MultiLevelSolution ml_sol_all_levels(& ml_mesh_all_levels);
 //Solution - END ==================

 //Problem - BEGIN ==================
  MultiLevelProblem  ml_prob_aux(ml_prob);
  ml_prob_aux.SetMultiLevelMeshAndSolution(& ml_sol_all_levels);
 //Problem - END ==================
               
               
               for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);  //We have to do so to avoid buildup of AddSolution with different FE families
               ml_sol_all_levels.set_analytical_function(unknowns[u]._name.c_str(), exact_sol[u]);   
               ml_sol_all_levels.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, & ml_prob_aux);
               }
               
               
 // Only for an EQUATION - BEGIN  ==================
         if (equation_solve) {
             
         ml_sol_all_levels.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
               for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels.GenerateBdc(unknowns[u]._name.c_str(), "Steady", & ml_prob_aux);
               }
               
         }
 // Only for an EQUATION - END  ==================

 return ml_sol_all_levels;
} 
    


//   print the error and the order of convergence between different levels
template < class type>
/*static*/  void FE_convergence< type >::output_convergence_order(const std::vector < std::vector < std::vector < type > > > &  norm,
                                    const unsigned int u,
                                    const unsigned int i,
                                    const unsigned int n) {

   if(i < norm[u].size() - 2)  {
//   std::cout << norm_name << " ERROR and ORDER OF CONVERGENCE: " << fam << " " << ord << "\n\n";

    std::cout << i + 1 << "\t\t";
    std::cout.precision(14);

    std::cout << norm[u][i][n] << "\t";

    std::cout << std::endl;
  
      std::cout.precision(3);
      std::cout << "\t\t";

        std::cout << log( norm[u][i][n] / norm[u][i + 1][n] ) / log(2.) << "  \t\t\t\t";

      std::cout << std::endl;
    }
    
  
  
}


template < class type>
/*static */ void FE_convergence< type >::output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < std::vector < std::vector < type > > > > &  norms, 
                                        const unsigned norm_flag,
                                        const unsigned volume_or_boundary) {
    
    
     std::cout << std::endl;
    
    const std::vector< std::string > norm_names = {"L2-NORM", "H1-SEMINORM"};
  
    for (unsigned int vb = 0; vb < volume_or_boundary + 1; vb++) {
        
     std::cout << "==== Volume = 0, boundary = 1: here we have " << vb << std::endl;
     
    assert( unknowns.size() == norms[vb].size() );
     for (unsigned int u = 0; u < unknowns.size(); u++) {
       for (int n = norm_flag; n >= 0; n--) {
            std::cout << unknowns[u]._name << " : " << norm_names[n] << " ERROR and ORDER OF CONVERGENCE"  << std::endl;
         for (int i = 0; i < norms[vb][u].size(); i++) {
                output_convergence_order(norms[vb], u, i, n);
            }
            std::cout << std::endl;
         }
       }
      }
      
 }
     


template < class type>
/*static */ void FE_convergence< type >::output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes , int loop_i) {

    std::cout << loop_i + 1 << "\t\t" <<  std::setw(11) << std::setprecision(10) << norm_i << "\t\t\t\t" ;
  
    if (loop_i < maxNumberOfMeshes/*norm.size()*/ - 2) {
      std::cout << std::setprecision(3) << log( norm_i/ norm_ip1 ) / log(2.) << std::endl;
    }
  
} 

 
 

template < class type>
/*static*/  std::vector< type > FE_convergence< type >::compute_error_norms_volume_or_boundary(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                                            const std::vector<Gauss> & quad_rules,
                                                                            const MultiLevelSolution* ml_sol,
                                                                            const MultiLevelSolution* ml_sol_all_levels,
                                                                            const std::string & unknown,
                                                                            const unsigned current_level,
                                                                            const unsigned norm_flag,
                                                                            const unsigned conv_rate_computation_method,
                                                                            const unsigned volume_or_boundary,
                                                                            const Math::Function< type > * ex_sol_in
                                             ) {
     
  // (//0 = only L2: //1 = L2 + H1)
  
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

   if (conv_rate_computation_method == 1 && ex_sol_in == NULL) { std::cout << "Please provide analytical solution" << std::endl; abort(); }
   
    
  const unsigned num_norms = norm_flag + 1;
  //norms that we are computing here //first L2, then H1 ============
  std::vector< type > norms_exact_function(num_norms);                  std::fill(norms_exact_function.begin(), norms_exact_function.end(), 0.);   
  std::vector< type > norms_exact_dofs(num_norms);       std::fill(norms_exact_dofs.begin(), norms_exact_dofs.end(), 0.);
  std::vector< type > norms_inexact_dofs(num_norms);     std::fill(norms_inexact_dofs.begin(), norms_inexact_dofs.end(), 0.);
  //norms that we are computing here //first L2, then H1 ============
  
  
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);
  const Solution* sol = ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 unsigned iproc = msh->processor_id();

  
 //***************************************************  
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = 3;

  // ======================================
  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
   
  //solution variable
  unsigned sol_uIndex = ml_sol->GetIndex(unknown.c_str()); // ml_sol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned sol_uType = ml_sol->GetSolutionType(sol_uIndex);    // get the finite element type for "u"

  

  

  CurrentElem < double > geom_element_iel(dim, msh);

  
  //-----------------  V initialization - BEGIN
 //***************************************************  
   std::vector < std::vector < double > >  JacI_qp(space_dim);
   std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    type detJac_qp = 0.;
  type weight = 0.; // gauss point weight
 //***************************************************  

  vector < type > phi_coords;
  vector < type > phi_coords_x;

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * dim_offset_grad);

  vector < type > phi;
  vector < type > phi_x;
  
  phi.reserve(max_size);
  phi_x.reserve(max_size * dim_offset_grad);
  


  vector < type >  sol_u;                               sol_u.reserve(max_size);
  vector < type >  sol_u_exact_at_dofs;   sol_u_exact_at_dofs.reserve(max_size);
  vector < type >  sol_u_inexact_prolongated_from_coarser_level;     sol_u_inexact_prolongated_from_coarser_level.reserve(max_size);

  unsigned solType_coords = BIQUADR_FE;
  vector < vector < type > > x(dim_offset_grad);    // local coordinates
  for (unsigned i = 0; i < x.size(); i++)   x[i].reserve(max_size);

//-----------------  V initialization - END

  
  
  
//-----------------  B initialization - BEGIN
    
 //***************************************************  
     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    type detJac_iqp_bdry = 0.;
  type weight_iqp_bdry = 0.;
 //***************************************************  
    
  std::vector <double> phi_u_bdry;  
  std::vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);

//-----------------  B initialization - END

  
  
  
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

//----------------
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
//----------------
    
    
    unsigned nDofu  = msh->GetElementDofNumber(iel, sol_uType);
    unsigned nDofx = msh->GetElementDofNumber(iel, solType_coords);

    // resize local arrays
    sol_u.resize(nDofu);
    sol_u_exact_at_dofs.resize(nDofu);
    sol_u_inexact_prolongated_from_coarser_level.resize(nDofu);

    for (int i = 0; i < x.size(); i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, solType_coords);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim_offset_grad/*dim*/; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
    

    
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        
        std::vector< type > x_at_dof(dim_offset_grad, 0.);
        for (unsigned jdim = 0; jdim < x_at_dof.size(); jdim++) x_at_dof[jdim] = x[jdim][i];
        
      unsigned solDof = msh->GetSolutionDof(i, iel, sol_uType);
                   sol_u[i]  =                                        (*sol->_Sol[sol_uIndex])(solDof);
      sol_u_inexact_prolongated_from_coarser_level[i]  = (*ml_sol_all_levels->GetSolutionLevel(current_level)->_Sol[sol_uIndex])(solDof);
      if (ex_sol_in != NULL) sol_u_exact_at_dofs[i] = ex_sol_in->value(x_at_dof);
    }


    // *** BOUNDARY - BEGIN ***
if (volume_or_boundary == 1 )	{  
	       
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
// ----------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
// ----------

		            const int bdry_index_j = msh->el->GetFaceElementIndex(iel, iface);

	    if(  bdry_index_j < 0 ) {

	
		
		for(unsigned ig_bdry = 0; ig_bdry < quad_rules[ielGeom_bdry].GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_iqp_bdry, space_dim);
    weight_iqp_bdry = detJac_iqp_bdry * quad_rules[ielGeom_bdry].GetGaussWeightsPointer()[ig_bdry];
    elem_all[ielGeom_bdry][sol_uType] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_u_bdry, phi_u_x_bdry, boost::none, space_dim);

    elem_all[ielGeom_bdry][solType_coords]->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_coords, phi_coords_x,  boost::none, space_dim);

		  
		 //========== QP eval - BEGIN ===============================================
      type		  sol_u_bdry_gss = 0.;
      type exactSol_from_dofs_gss = 0.;
      type sol_u_inexact_prolongated_from_coarser_level_gss = 0.;
          
      std::vector< type > sol_u_x_bdry_gss(space_dim, 0.);
      std::vector < type > gradSolu_exact_at_dofs_bdry_gss(dim_offset_grad, 0.);
      std::vector < type > gradSolu_inexact_prolongated_from_coarser_level_bdry_gss(dim_offset_grad, 0.);

                  
		      for (int i_bdry = 0; i_bdry < phi_u_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_u_bdry_gss                                    +=  sol_u[i_vol]              * phi_u_bdry[i_bdry];
            exactSol_from_dofs_gss                            += sol_u_exact_at_dofs[i_vol] * phi_u_bdry[i_bdry];
            
            sol_u_inexact_prolongated_from_coarser_level_gss  += sol_u_inexact_prolongated_from_coarser_level[i_vol]              * phi_u_bdry[i_bdry];
        
            for (int d = 0; d < space_dim; d++) {
			      sol_u_x_bdry_gss[d] += sol_u[i_vol] * phi_u_x_bdry[i_bdry * space_dim + d];
			      gradSolu_inexact_prolongated_from_coarser_level_bdry_gss[d] += sol_u_inexact_prolongated_from_coarser_level[i_vol] * phi_u_x_bdry[i_bdry * space_dim + d];
			    }
		    }
		      
		      double laplace_u_bdry = 0.;  for (int d = 0; d < space_dim; d++) { laplace_u_bdry += sol_u_x_bdry_gss[d] * sol_u_x_bdry_gss[d]; }


       std::vector < type > x_gss_bdry(dim_offset_grad, 0.);
       
    for (unsigned i = 0; i < phi_coords.size(); i++) {
          for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
           x_gss_bdry[jdim] += geom_element_iel.get_coords_at_dofs_bdry_3d()[jdim][i] * phi_coords[i];
          }
       }
		 //========== QP eval - END ===============================================

// H^0 ==============      
      type exactSol_bdry = 0.; if (ex_sol_in != NULL) exactSol_bdry = ex_sol_in->value(x_gss_bdry);
                 norms_exact_function[0] += (sol_u_bdry_gss - exactSol_bdry) * (sol_u_bdry_gss - exactSol_bdry) * weight_iqp_bdry; 
//                  norms_exact_dofs[0]     += (sol_u_gss - exactSol_from_dofs_gss)  * (sol_u_gss - exactSol_from_dofs_gss) * weight;

      norms_inexact_dofs[0]   += (sol_u_bdry_gss - sol_u_inexact_prolongated_from_coarser_level_gss)   * (sol_u_bdry_gss - sol_u_inexact_prolongated_from_coarser_level_gss)  * weight_iqp_bdry;
                 norms_inexact_dofs[0]   += (sol_u_bdry_gss - exactSol_bdry) * (sol_u_bdry_gss - exactSol_bdry) * weight_iqp_bdry; 
                  
                  
                  
// H^1 ==============      
    /*else*/ if (norm_flag == 1) {
      std::vector < type > exactGradSol_bdry(dim_offset_grad, 0.);    if (ex_sol_in != NULL) exactGradSol_bdry = ex_sol_in->gradient(x_gss_bdry);  
      ///@todo THIS IS WRONG! It is NOT the SURFACE GRADIENTTTT, so we have to be careful when we have more than one boundary face!!!

      for (unsigned d = 0; d < dim_offset_grad ; d++) {
        norms_exact_function[1] += ( (sol_u_x_bdry_gss[d] - exactGradSol_bdry[d])               * (sol_u_x_bdry_gss[d]  - exactGradSol_bdry[d]) ) * weight_iqp_bdry;
        
//         norms_exact_dofs[1]     += ((gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d]) * (gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d])) * weight_iqp_bdry;
        norms_inexact_dofs[1]   += ((sol_u_x_bdry_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_bdry_gss[d])  * (sol_u_x_bdry_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_bdry_gss[d])) * weight_iqp_bdry;
      }
   }
   
                  
        }
            
        
	      } //end face
	      
	  }  // loop over element faces   
	  

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
     
    
}   
    // *** BOUNDARY - END ***
  
    // *** VOLUME - BEGIN ***
  else if (volume_or_boundary == 0) {  
      
      
    const short unsigned ielGeom = geom_element_iel.geom_type();
    
   
    for (unsigned ig = 0; ig < quad_rules[ielGeom].GetGaussPointsNumber(); ig++) {

      // *** get gauss point weight, test function and test function partial derivatives ***
	elem_all[ielGeom][solType_coords]->JacJacInv(x, ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    elem_all[ielGeom][sol_uType]->shape_funcs_current_elem(ig, JacI_qp, phi, phi_x, boost::none, space_dim);
    
    elem_all[ielGeom][solType_coords]->shape_funcs_current_elem(ig, JacI_qp, phi_coords, phi_coords_x,  boost::none, space_dim);
    
    weight = detJac_qp * quad_rules[ielGeom].GetGaussWeightsPointer()[ig];

    
		 //========== QP eval - BEGIN ===============================================
      type sol_u_gss = 0.;
      type exactSol_from_dofs_gss = 0.;
      type sol_u_inexact_prolongated_from_coarser_level_gss = 0.;
      
      std::vector < type > gradSolu_gss(dim_offset_grad, 0.);
      std::vector < type > gradSolu_exact_at_dofs_gss(dim_offset_grad, 0.);
      std::vector < type > gradSolu_inexact_prolongated_from_coarser_level_gss(dim_offset_grad, 0.);


      for (unsigned i = 0; i < phi.size(); i++) {
        sol_u_gss                += phi[i] * sol_u[i];
        exactSol_from_dofs_gss  += phi[i] * sol_u_exact_at_dofs[i];
        sol_u_inexact_prolongated_from_coarser_level_gss   += phi[i] * sol_u_inexact_prolongated_from_coarser_level[i];

        for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
          gradSolu_gss[jdim]                += phi_x[i * dim_offset_grad + jdim] * sol_u[i];
          gradSolu_exact_at_dofs_gss[jdim]  += phi_x[i * dim_offset_grad + jdim] * sol_u_exact_at_dofs[i];
          gradSolu_inexact_prolongated_from_coarser_level_gss[jdim]   += phi_x[i * dim_offset_grad + jdim] * sol_u_inexact_prolongated_from_coarser_level[i];
        }
        
 
      }
      
      
      std::vector < type > x_gss(dim_offset_grad, 0.);
      
       for (unsigned i = 0; i < phi_coords.size(); i++) {
          for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
           x_gss[jdim] += x[jdim][i] * phi_coords[i];
          }
       }
		 //========== QP eval - END ===============================================
      
      
      
      

// H^0 ==============      
//     if (norm_flag == 0) {
      type exactSol = 0.; if (ex_sol_in != NULL) exactSol = ex_sol_in->value(x_gss);

      norms_exact_function[0] += (sol_u_gss - exactSol)                * (sol_u_gss - exactSol)       * weight;
      norms_exact_dofs[0]     += (sol_u_gss - exactSol_from_dofs_gss)  * (sol_u_gss - exactSol_from_dofs_gss) * weight;
      norms_inexact_dofs[0]   += (sol_u_gss - sol_u_inexact_prolongated_from_coarser_level_gss)   * (sol_u_gss - sol_u_inexact_prolongated_from_coarser_level_gss)  * weight;
//     }
    
// H^1 ==============      
    /*else*/ if (norm_flag == 1) {
      std::vector < type > exactGradSol(dim_offset_grad,0.);    if (ex_sol_in != NULL) exactGradSol = ex_sol_in->gradient(x_gss);

      for (unsigned d = 0; d < dim_offset_grad ; d++) {
        norms_exact_function[1] += ((gradSolu_gss[d] - exactGradSol[d])               * (gradSolu_gss[d]  - exactGradSol[d])) * weight;
        norms_exact_dofs[1]     += ((gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d]) * (gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d])) * weight;
        norms_inexact_dofs[1]   += ((gradSolu_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_gss[d])  * (gradSolu_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_gss[d]))  * weight;
      }
   }
      
   } // end gauss point loop
  }
    // *** VOLUME - END ***
  
  
  
  } //end element loop for each process

  
  // add the norms of all processes
  NumericVector* norm_vec_exact_using_qp;
                 norm_vec_exact_using_qp = NumericVector::build().release();
                 norm_vec_exact_using_qp->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_exact_using_qp->set(iproc, norms_exact_function[0]);  norm_vec_exact_using_qp->close();  norms_exact_function[0] = norm_vec_exact_using_qp->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_exact_using_qp->set(iproc, norms_exact_function[1]);  norm_vec_exact_using_qp->close();  norms_exact_function[1] = norm_vec_exact_using_qp->l1_norm(); }

          delete norm_vec_exact_using_qp;

   // add the norms of all processes
  NumericVector* norm_vec_exact_dofs;
                 norm_vec_exact_dofs = NumericVector::build().release();
                 norm_vec_exact_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_exact_dofs->set(iproc, norms_exact_dofs[0]);  norm_vec_exact_dofs->close();  norms_exact_dofs[0] = norm_vec_exact_dofs->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_exact_dofs->set(iproc, norms_exact_dofs[1]);  norm_vec_exact_dofs->close();  norms_exact_dofs[1] = norm_vec_exact_dofs->l1_norm(); }

          delete norm_vec_exact_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
                 norm_vec_inexact = NumericVector::build().release();
                 norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_inexact->set(iproc, norms_inexact_dofs[0]);  norm_vec_inexact->close();  norms_inexact_dofs[0] = norm_vec_inexact->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_inexact->set(iproc, norms_inexact_dofs[1]);  norm_vec_inexact->close();  norms_inexact_dofs[1] = norm_vec_inexact->l1_norm(); }

          delete norm_vec_inexact;

          
    for (int n = 0; n < norms_exact_function.size(); n++) { 
   norms_exact_function[n] = sqrt(norms_exact_function[n]);                                            
       norms_exact_dofs[n] = sqrt(norms_exact_dofs[n]);                                            
     norms_inexact_dofs[n] = sqrt(norms_inexact_dofs[n]);
    }
    
    
if (conv_rate_computation_method == 0)  return norms_inexact_dofs;
if (conv_rate_computation_method == 1)  return norms_exact_function;
 //   return norms_exact_dofs;

 
} 
 
 
 

     
template < class type>
/*static*/  void FE_convergence< type >::compute_error_norms_per_unknown_per_level(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                                                   const std::vector<Gauss> & quad_rules,
                                                                                   const MultiLevelSolution* ml_sol_single_level,
                                                                                   MultiLevelSolution* ml_sol_all_levels,
                                                                                   const std::vector< Unknown > &  unknowns,
                                                                                   const unsigned i,
                                                                                   const unsigned norm_flag,
                                                                                   std::vector < std::vector < std::vector < std::vector < type > > > > &  norms,
                                                                                   const unsigned conv_rate_computation_method,
                                                                                   const unsigned volume_or_boundary,
                                                                                   const  std::vector< Math::Function< type > * > &  ex_sol_in
                                         ) {
     
    
              if (volume_or_boundary == 0) {
            }
            else if (volume_or_boundary == 1) {
            }
            else { std::cout << "Boundary of boundary not implemented in function " << __func__; abort(); }
   
        if ( i > 0 ) {

            // ======= prolongate to the current level (i) from the coarser level (i-1) (so that you can compare the two) ========================
            ml_sol_all_levels->RefineSolution(i);
            
            // =======  compute the error norm at the current level (i) ========================
    for (unsigned int vb = 0; vb < volume_or_boundary + 1; vb++) { //0: volume, 1: boundary
            for (unsigned int u = 0; u < unknowns.size(); u++) {  //this loop could be inside the below function
                
            std::vector< type > norm_out;
            
            norm_out = FE_convergence::compute_error_norms_volume_or_boundary (elem_all, quad_rules, ml_sol_single_level, ml_sol_all_levels, unknowns[u]._name, i, norm_flag, conv_rate_computation_method, vb,  ex_sol_in[u]);
            
            
              for (int n = 0; n < norms[vb][u][i-1].size(); n++)      norms[vb][u][i-1][n] = norm_out[n];
                                       
                   }
                   
          }
        }         
                 
              // ======= store the last computed solution to prepare the next iteration (the current level i is now overwritten) ========================
            ml_sol_all_levels->fill_at_level_from_level(i, ml_sol_single_level->_mlMesh->GetNumberOfLevels() - 1, *ml_sol_single_level);
        
                 
                 
                 
 } 
 







}  //end namespace









#endif
