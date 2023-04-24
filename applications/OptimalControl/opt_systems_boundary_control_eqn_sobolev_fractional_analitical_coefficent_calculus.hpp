#ifndef OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_
#define OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_



namespace femus {

namespace ctrl {


       //calculation c parameter - BEGIN
       double calculation_c_parameter(double first_component_in_normal_direction_to_straigh_line){
           if(first_component_in_normal_direction_to_straigh_line > 0) { return - 0.75;}
           else { return  - 0.25;}
       }
       //calculation c parameter - END



      // give a vector return tur if the vector components are the same , false if not - BEGIN
       template < class component_type >
       bool vector_components_are_equal (std::vector<component_type> my_vector){
                    bool are_equalal = true;
                    for(int d = 0; d < my_vector.size() - 1; d++){
                       if( my_vector[d] != my_vector[d + 1] ) { are_equalal = false; break; }
                    }
                    return are_equalal;
       }
      // give a vector return tur if the vector components are the same , false if not - END



      //calculus of parameter for analitical solution of mixed integral - BEGIN
      void calculation_of_coefficent_of_analitical_solution(std::vector  < std::vector  <  double> > radius_centered_at_x_qp_of_iface_bdry_bdry,
                                                            std::vector< unsigned > global_dirs_for_atan,
                                                            double& a,
                                                            double& b,
                                                            double& c){
            if( vector_components_are_equal<double>( radius_centered_at_x_qp_of_iface_bdry_bdry[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ] ) ) {
                a = 1.;
                b = 0.;
                c = calculation_c_parameter( radius_centered_at_x_qp_of_iface_bdry_bdry[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][ 0 ] );
            }
            else if( vector_components_are_equal<double>( radius_centered_at_x_qp_of_iface_bdry_bdry[ global_dirs_for_atan[ 1 /*global_dir_second */] ] ) ) {
                a = 0.;
                b = 1.;
                c = calculation_c_parameter( radius_centered_at_x_qp_of_iface_bdry_bdry[ global_dirs_for_atan[ 1 /*global_dir_second*/ ] ][ 0 ] );
            }
      }

      //calculus of parameter for analitical solution of mixed integral - END




 } //end namespace ctrl
} //end namespace femus


#endif
