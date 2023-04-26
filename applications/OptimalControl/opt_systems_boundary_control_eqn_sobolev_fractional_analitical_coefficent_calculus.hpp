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

//************************* give a vector return if the vector has the component increasing along the integration direction - BEGIN *************************
bool is_vector_oriented_with_cords_decreasing_along_integration_line_direction(std::vector< std::vector< double > > vector_of_cords, unsigned int direction_of_bndry_bndry_integration_line){

    if( vector_of_cords[direction_of_bndry_bndry_integration_line][1] < vector_of_cords[direction_of_bndry_bndry_integration_line][0]) {return true;}
    else {return false;}

}
//************************* give a vector return if the vector has the component increasing along the integration direction - END *************************





//************************* integration_extreme_cords - BEGIN *************************
      std::vector<std::vector<double>> integration_extreme_cords(const unsigned dim,
                                                                 const unsigned boundary_dim,
                                                                 std::vector< std::vector< double > >face_index_bndry_region_vertex_cords,
                                                                 std::vector< std::vector< double > > nodes_on_line_of_bndry_bndry,
//                                                                  std::vector< unsigned > global_dirs_for_atan,
                                                                 bool are_nodes_coordinates_decreasing_along_integration_directin,
                                                                 unsigned int direction_of_bndry_bndry_integration_line,
                                                                 unsigned int normal_tangential_direction_to_bndry_bndry_integration_line
//                                                                  ,const unsigned int face_index,
                                                                 ){


    constexpr unsigned first_node_along_line_of_bndry_bndry = 0;

    //----------------- matrix preparation - BEGIN -----------------
    std::vector<std::vector<double>> cords_of_analitical_integer_extreme(dim);
    for(unsigned f = 0; f <  cords_of_analitical_integer_extreme.size(); f++) {
             cords_of_analitical_integer_extreme[f].resize(boundary_dim);
         }
    //----------------- matrix preparation - END -----------------

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
    bool is_integration_extreme_matrix_decreasing = is_vector_oriented_with_cords_decreasing_along_integration_line_direction(cords_of_analitical_integer_extreme, direction_of_bndry_bndry_integration_line);
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
                                                            std::vector< std::vector< double > >face_index_bndry_region_vertex_cords,
                                                            std::vector< std::vector< double > > nodes_on_line_of_bndry_bndry,
                                                            std::vector< unsigned int> global_dirs_for_atan,
//                                                             std::vector< std::vector< double > > extreeme_bndry_bndry_points_coords,
                                                            double& a,
                                                            double& b,
                                                            double& c){

//             double first_direction_point_one;
//             double second_direction_point_one;
//             double first_direction_point_two;
//             double second_direction_point_two;

                constexpr unsigned global_dir_first = 0;
                constexpr unsigned global_dir_second = 1;

                constexpr unsigned first_node_along_line_of_bndry_bndry = 0;

                bool are_nodes_coordinates_decreasing_along_integration_directin = false;

            std::vector<std::vector<double>> analitical_integral_extreme_cords;

            if( vector_components_are_equal<double>( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ global_dir_first ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_first ];
                unsigned int direction_of_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_second ];

                a = 1.;
                b = 0.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_integration_line_direction(nodes_on_line_of_bndry_bndry, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(dim,
                                                                             dim_bdry,
                                                                             face_index_bndry_region_vertex_cords,
                                                                             nodes_on_line_of_bndry_bndry,
//                                                                              global_dirs_for_atan,
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line
//                                                                              ,face_index,
                                                                             );

                c = calculation_c_parameter( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[  global_dir_first ] ][ first_node_along_line_of_bndry_bndry ] );

//                 extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][0] = nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][ 0 ];
//                 extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][1] = extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][0];
//
//                 extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 1 /*global_dir_second*/ ] ][0] = nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][ 0 ];
//                 extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 1 /*global_dir_second*/ ] ][1] = extreeme_bndry_bndry_points_coords[ global_dirs_for_atan[ 0 /*global_dir_first*/ ] ][0];

            }
            else if( vector_components_are_equal<double>( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[  global_dir_second ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_second ];
                unsigned int direction_of_bndry_bndry_integration_line = global_dirs_for_atan[  global_dir_first ];

                a = 0.;
                b = 1.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_integration_line_direction(nodes_on_line_of_bndry_bndry, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(dim,
                                                                             dim_bdry,
                                                                             face_index_bndry_region_vertex_cords,
                                                                             nodes_on_line_of_bndry_bndry,
//                                                                              global_dirs_for_atan,
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line
//                                                                              ,face_index,
                                                                             );

                c = calculation_c_parameter( nodes_on_line_of_bndry_bndry[ global_dirs_for_atan[ global_dir_second ] ][ first_node_along_line_of_bndry_bndry ] );
            }
      }
//************************* calculus of parameter for analitical solution of mixed integral - END *************************




 } //end namespace ctrl
} //end namespace femus


#endif
