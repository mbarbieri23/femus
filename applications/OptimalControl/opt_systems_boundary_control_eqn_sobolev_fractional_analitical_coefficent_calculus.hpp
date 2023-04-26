#ifndef OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_
#define OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_

    // 1) togliere estremi 0.25 e 0.75


namespace femus {

namespace ctrl {


//************************* c parameter calculation - BEGIN *************************
       double calculation_c_parameter(double first_nodes_of_bdryt_bdry_component_in_normal_direction_to_straigh_line){    //the normal components (respect the integration line) of the points along the line are equal !!!!!!
           if(first_nodes_of_bdryt_bdry_component_in_normal_direction_to_straigh_line > 0.25) { return - 0.75;}
           else { return  - 0.25;}
       }
//************************* c parameter calculation - END *************************


//************************* give a vector return tur if the vector components are the same , false if not - BEGIN *************************
       template < class component_type >
       bool vector_components_are_equal (std::vector<component_type> my_vector){
                    bool are_equalal = true;
                    for(int d = 0; d < my_vector.size() - 1; d++){
                       if( my_vector[d] != my_vector[d + 1] ) { are_equalal = false; break; }
                    }
                    return are_equalal;
       }
//************************* give a vector return tur if the vector components are the same , false if not - END *************************


//************************* give a vector return if the vector has the component decreasing along specific direction - BEGIN *************************
bool is_vector_oriented_with_cords_decreasing_along_specific_direction(std::vector< std::vector< double > > vector_of_cords, unsigned int direction){

    if( vector_of_cords[direction][1] < vector_of_cords[direction][0]) {return true;}
    else {return false;}

}
//************************* give a vector return if the vector has the component decreasing along specific direction - END *************************


//************************* integration_extreme_cords - BEGIN *************************
     std::vector<std::vector<double>> integration_extreme_cords(const unsigned dim,
                                    const unsigned boundary_dim,
                                    //-----vector & matrix -------
                                    std::vector< std::vector< double > >face_index_bndry_region_vertex_cords,
                                    std::vector< std::vector< double > > nodes_on_line_of_bndry_bndry,
                                    //------- direction of tangential vector -------
                                    unsigned int direction_of_bndry_bndry_integration_line,
                                    unsigned int normal_tangential_direction_to_bndry_bndry_integration_line,
                                    //----------- orientation ------------
                                    bool are_nodes_coordinates_decreasing_along_integration_directin,
                                    //--------- matrix of extreeme of control region --------------
                                    std::vector<std::vector< double > > & cords_of_analitical_integer_extreme){


    constexpr unsigned first_node_along_line_of_bndry_bndry = 0;

    unsigned int integration_extreme_point_index = 0;


    //----------------- build matrix - BEGIN -----------------
    for(unsigned pt = 0; pt <  face_index_bndry_region_vertex_cords[normal_tangential_direction_to_bndry_bndry_integration_line].size(); pt++) {

        if( face_index_bndry_region_vertex_cords[normal_tangential_direction_to_bndry_bndry_integration_line][pt] ==
            nodes_on_line_of_bndry_bndry[normal_tangential_direction_to_bndry_bndry_integration_line][first_node_along_line_of_bndry_bndry]) {

            for(unsigned dir = 0; dir <  cords_of_analitical_integer_extreme.size(); dir++) {

                cords_of_analitical_integer_extreme[dir][integration_extreme_point_index] = face_index_bndry_region_vertex_cords[dir][pt];

            }
         integration_extreme_point_index++; //add +1 only if you add something in the matrix cords_of_analitical_integer_extreme

        }

    }
    //----------------- build matrix - END -----------------

    //----------------- reverse matrix - BEGIN -----------------
    bool is_integration_extreme_matrix_decreasing = is_vector_oriented_with_cords_decreasing_along_specific_direction(cords_of_analitical_integer_extreme, direction_of_bndry_bndry_integration_line);
    bool is_reserve_vector_needed = ( are_nodes_coordinates_decreasing_along_integration_directin ^ is_integration_extreme_matrix_decreasing );

    if( is_reserve_vector_needed ){
            for(unsigned dir = 0; dir <  cords_of_analitical_integer_extreme.size(); dir++) {
             std::reverse( cords_of_analitical_integer_extreme[dir].begin(), cords_of_analitical_integer_extreme[dir].end() );
            }
    }
    //----------------- reverse matrix - END -----------------

    return cords_of_analitical_integer_extreme;

}
//************************* integration_extreme_cords - END *************************


//************************* calculus of parameter for analitical solution of mixed integral - BEGIN *************************
      void calculation_of_coefficent_of_analitical_solution(const unsigned int dim,
                                                            const unsigned int dim_bdry,
                                                            //-----vector & matrix -------
                                                            std::vector< std::vector< double > >face_index_bndry_region_vertex_cords,
                                                            std::vector< std::vector< double > > nodes_on_line_of_bndry_bndry,
                                                            //------- tangent vector -------
                                                            std::vector< unsigned int> global_dirs_for_atan,
                                                            //------- output -------
                                                            std::vector< std::vector< double > >& analitical_integral_extreme_cords,
                                                            double& a, double& b, double& c){

                constexpr unsigned global_dir_first = 0;
                constexpr unsigned global_dir_second = 1;

                constexpr unsigned first_node_along_line_of_bndry_bndry = 0;

                bool are_nodes_coordinates_decreasing_along_integration_directin = false;


            if( vector_components_are_equal<double>( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ global_dir_first ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_first ];
                unsigned int direction_of_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_second ];

                a = 1.;
                b = 0.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_specific_direction(nodes_on_line_of_bndry_bndry, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(dim,
                                                                             dim_bdry,
                                                                             //-----vector & matrix -------
                                                                             face_index_bndry_region_vertex_cords,
                                                                             nodes_on_line_of_bndry_bndry,
                                                                             //------- direction of tangential vector -------
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line,
                                                                             //----------- orientation ------------
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             //--------- matrix of extreeme of control region --------------
                                                                             analitical_integral_extreme_cords);

                c = calculation_c_parameter( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[  global_dir_first ] ][ first_node_along_line_of_bndry_bndry ] );

            }
            else if( vector_components_are_equal<double>( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[  global_dir_second ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_second ];
                unsigned int direction_of_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_first ];

                a = 0.;
                b = 1.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_specific_direction(nodes_on_line_of_bndry_bndry, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(dim,
                                                                             dim_bdry,
                                                                             //-----vector & matrix -------
                                                                             face_index_bndry_region_vertex_cords,
                                                                             nodes_on_line_of_bndry_bndry,
                                                                             //------- direction of tangential vector -------
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line,
                                                                             //----------- orientation ------------
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             //--------- matrix of extreeme of control region --------------
                                                                             analitical_integral_extreme_cords);

                c = calculation_c_parameter( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ global_dir_second ] ][ first_node_along_line_of_bndry_bndry ] );
            }
      }
//************************* calculus of parameter for analitical solution of mixed integral - END *************************

// geom_element_iel.get_elem_center_bdry_3d()[direction]
//************************* Extreme of analitical integration - BEGIN *************************
//      std::vector<std::vector<double>> extreme_of_analitical_integration(//------- direction of tangential vector -------
//                                                                         const unsigned first_tangential_direction,
//                                                                         const unsigned second_tangential_direction,
//                                                                         //------- vector of nodes along element face
//                                                                         vector< vector< double > > nodes_along_line,
//                                                                         //------- extreme fo the integration
//                                                                         std::vector<std::vector<double>> cords_of_analitical_integer_extreme){
//
//
//
//     //----------------- build matrix - BEGIN -----------------
//     for(unsigned pt = 0; pt <  face_index_bndry_region_vertex_cords[normal_tangential_direction_to_bndry_bndry_integration_line].size(); pt++) {
//
//         if( face_index_bndry_region_vertex_cords[normal_tangential_direction_to_bndry_bndry_integration_line][pt] ==
//             nodes_on_line_of_bndry_bndry[normal_tangential_direction_to_bndry_bndry_integration_line][first_node_along_line_of_bndry_bndry]) {
//
//             for(unsigned dir = 0; dir <  cords_of_analitical_integer_extreme.size(); dir++) {
//
//                 cords_of_analitical_integer_extreme[dir][integration_extreme_point_index] = face_index_bndry_region_vertex_cords[dir][pt];
//
//             }
//          integration_extreme_point_index++; //add +1 only if you add something in the matrix cords_of_analitical_integer_extreme
//
//         }
//
//     }
//     //----------------- build matrix - END -----------------
//
//     //----------------- reverse matrix - BEGIN -----------------
//     bool is_integration_extreme_matrix_decreasing = is_vector_oriented_with_cords_decreasing_along_specific_direction(cords_of_analitical_integer_extreme, direction_of_bndry_bndry_integration_line);
//     bool is_reserve_vector_needed = ( are_nodes_coordinates_decreasing_along_integration_directin ^ is_integration_extreme_matrix_decreasing );
//
//     if( is_reserve_vector_needed ){
//             for(unsigned dir = 0; dir <  cords_of_analitical_integer_extreme.size(); dir++) {
//              std::reverse( cords_of_analitical_integer_extreme[dir].begin(), cords_of_analitical_integer_extreme[dir].end() );
//             }
//     }
//     //----------------- reverse matrix - END -----------------
//
//     return cords_of_analitical_integer_extreme;
//
// }
//************************* Extreme of analitical integration - END *************************

 } //end namespace ctrl

} //end namespace femus


#endif
