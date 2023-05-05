#ifndef OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_
#define OPT_SYSTEMS_BOUNDARY_CONTROL_EQN_SOBOLEV_FRACTIONAL_ANALITICAL_COEFFICENT_CALCULUS_HPP_

    // 1) togliere estremi 0.25 e 0.75


namespace femus {

namespace ctrl {

////////LIST_OF_CTRL_FACES stuff - BEGIN
//************************* get_position_on_list_faces_clas - BEGIN *************************
           template < class LIST_OF_CTRL_FACES >
           unsigned int get_position_on_list_faces_clas( unsigned int face_index){
               unsigned int number_of_face = LIST_OF_CTRL_FACES:: _face_with_extremes_index_size;
               for(int t = 0; t < number_of_face; t++){
                   if(face_index ==  LIST_OF_CTRL_FACES::_face_with_extremes_index[t]){return t;}
               }
           }
//************************* get_position_on_list_faces_clas - END *************************

//************************* get_extreme_for_each_direction - BEGIN *************************
           template < class LIST_OF_CTRL_FACES >
           std::vector<std::vector< double > > get_extreme_for_each_direction( unsigned int face_index,
                                                                               const std::vector< unsigned int> tangential_vector,
                                                                               unsigned int normal_dir){

               unsigned int position = get_position_on_list_faces_clas< LIST_OF_CTRL_FACES >(face_index);
               unsigned int dim = LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face + 1;
               unsigned int number_of_extreme_for_tng_direction = LIST_OF_CTRL_FACES :: _num_of_control_extremes_per_tang_comp_per_face;

//                std::vector< unsigned int > tangential_vector = LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control( face_index, dim-1);
               bool is_face_on_maximal_coordinate = LIST_OF_CTRL_FACES ::is_face_on_maximal_coordinate_by_tangential_component( tangential_vector );

               //************** declaration - BEGIN **************
               std::vector< std::vector< double > > face_index_extreme_of_contrl_region( dim, std::vector< double >(number_of_extreme_for_tng_direction) );
//                for(int dir = 0; dir < dim; dir++){
//                    face_index_extreme_of_contrl_region[dir].resize(number_of_extreme_for_tng_direction,0);
//                }
               //************** declaration - END **************

               //************** built - BEGIN **************

               //----------------- filling the matrix of bndry bndry cords along normal direction - BEGIN -----------------
//                unsigned int normal_dir =  LIST_OF_CTRL_FACES :: normal_direction_to_Gamma_control( face_index );
               for (int p = 0; p < face_index_extreme_of_contrl_region[ normal_dir ].size(); p++){
                   face_index_extreme_of_contrl_region[ normal_dir ][ p ] =LIST_OF_CTRL_FACES :: normal_coordinate(is_face_on_maximal_coordinate);
               }
              //----------------- filling the matrix of bndry bndry cords along normal direction - END -----------------

              //----------------- filling the matrix of bndry bndry cords along tamgemtial direction - BEGIN -----------------
              for(int t = 0; t < tangential_vector.size(); t++){
                  for (int p = 0; p < face_index_extreme_of_contrl_region[ normal_dir ].size(); p++){
                      face_index_extreme_of_contrl_region[ tangential_vector[t] ][p] = LIST_OF_CTRL_FACES::_face_with_extremes_extremes_on_tang_surface[position][t][p];
                }
            }
               //----------------- filling the matrix of bndry bndry cords along tamgemtial direction - END -----------------

               //************** built - END **************
               return face_index_extreme_of_contrl_region;
           }
//************************* get_extreme_for_each_direction - END *************************
////////LIST_OF_CTRL_FACES stuff - END


//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//********************************************* PREVIOUS FUNCTION - BEGIN *****************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************

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


//************************* integration_extreme_cords 3 x 3 - BEGIN *************************
     std::vector<std::vector<double>> integration_extreme_cords(unsigned int face_index,
                                                                //const unsigned dim,
                                                                //const unsigned boundary_dim,
                                                                //-----element face center
                                                                std::vector<double> element_face_center_3d,
                                                                //-----vector & matrix -------
                                                                std::vector< std::vector< double > >face_index_extreme_3d_along_direction,
                                                                std::vector< std::vector< double > > nodes_on_element_boundary_face_line,
                                                                //------- direction of tangential vector -------
                                                                unsigned int direction_of_bndry_bndry_integration_line,
                                                                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line,
                                                                //------- normal dir -------
                                                                unsigned int normal_dir,
                                                                //----------- orientation ------------
                                                                bool are_nodes_coordinates_decreasing_along_integration_directin,
                                                                //--------- matrix of extreeme of control region --------------
                                                                std::vector<std::vector< double > > & cords_of_analitical_integer_extreme){


    constexpr unsigned first_node_along_element_face_line = 0;

    //----------------- filling extreme along normal tangential direction - BEGIN -----------------
    for (int p = 0; p < cords_of_analitical_integer_extreme[normal_dir].size(); p++){

        if(nodes_on_element_boundary_face_line[normal_tangential_direction_to_bndry_bndry_integration_line][first_node_along_element_face_line] > element_face_center_3d[normal_tangential_direction_to_bndry_bndry_integration_line] ){

            std::vector<double>::iterator result = std::max_element( face_index_extreme_3d_along_direction[normal_tangential_direction_to_bndry_bndry_integration_line].begin(),
                                                                     face_index_extreme_3d_along_direction[normal_tangential_direction_to_bndry_bndry_integration_line].end()   );

            cords_of_analitical_integer_extreme[normal_tangential_direction_to_bndry_bndry_integration_line][p] = *result;
        }
        else{

            std::vector<double>::iterator result = std::min_element( face_index_extreme_3d_along_direction[normal_tangential_direction_to_bndry_bndry_integration_line].begin(),
                                                                     face_index_extreme_3d_along_direction[normal_tangential_direction_to_bndry_bndry_integration_line].end()    );

            cords_of_analitical_integer_extreme[normal_tangential_direction_to_bndry_bndry_integration_line][p] = *result;
        }

    }
    //----------------- filling extreme along normal tangential direction - END -------------------

    //----------------- filling normal direction -----------------
    cords_of_analitical_integer_extreme[normal_dir] = face_index_extreme_3d_along_direction[normal_dir];
    //------------------------------------------------------------

    //----------------- filling integration line direction -----------------
    cords_of_analitical_integer_extreme[direction_of_bndry_bndry_integration_line] = face_index_extreme_3d_along_direction[direction_of_bndry_bndry_integration_line];
    //----------------------------------------------------------------------

    //----------------- reverse matrix - BEGIN -----------------
    bool is_integration_extreme_matrix_decreasing = is_vector_oriented_with_cords_decreasing_along_specific_direction(cords_of_analitical_integer_extreme, direction_of_bndry_bndry_integration_line);
    bool is_reserve_vector_needed = ( are_nodes_coordinates_decreasing_along_integration_directin ^ is_integration_extreme_matrix_decreasing );

    if( is_reserve_vector_needed ){
            for(unsigned dir = 0; dir <  cords_of_analitical_integer_extreme.size(); dir++) {
             std::reverse( cords_of_analitical_integer_extreme[dir].begin(), cords_of_analitical_integer_extreme[dir].end() );
            }
    }
    //----------------- reverse matrix - END -------------------

    return cords_of_analitical_integer_extreme;

}
//************************* integration_extreme_cords 3 x 3 - END *************************


//************************* calculus of parameter for analitical solution of mixed integral - BEGIN *************************
           template < class LIST_OF_CTRL_FACES >
      void calculation_of_coefficent_of_analitical_solution(unsigned int face_index,
                                                            //const unsigned int dim,
                                                            //const unsigned int dim_bdry,
                                                            //-----element face center
                                                            std::vector<double> element_face_center_3d,
                                                            //-----vector & matrix -------
                                                            std::vector< std::vector< double > > nodes_on_element_boundary_face_line,
                                                            //------- tangent vector -------
                                                            const std::vector< unsigned int> tangential_face_index_vector,
                                                            //------- normal dir -------
                                                            unsigned int normal_dir,
                                                            //------- output -------
                                                            std::vector< std::vector< double > >& analitical_integral_extreme_cords,
                                                            double& a, double& b, double& c){

                constexpr unsigned first_tang_dir = 0;
                constexpr unsigned second_tang_dir = 1;

                constexpr unsigned first_node_along_element_face_line = 0;

                bool are_nodes_coordinates_decreasing_along_integration_directin = false;

                std::vector<std::vector< double > > face_index_direction_extreme_bdry_region = get_extreme_for_each_direction< LIST_OF_CTRL_FACES >( face_index,
                                                                                                                                                     tangential_face_index_vector,
                                                                                                                                                     normal_dir);


            if( vector_components_are_equal<double>( nodes_on_element_boundary_face_line[ tangential_face_index_vector[ first_tang_dir ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = tangential_face_index_vector[  first_tang_dir ];// normal is the components (respect the integration line) where points  are equal !!!!!!
                unsigned int direction_of_bndry_bndry_integration_line = tangential_face_index_vector[  second_tang_dir ];

                a = 1.;
                b = 0.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_specific_direction(nodes_on_element_boundary_face_line, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(face_index,
                                                                             //dim,
                                                                             //dim_bdry,
                                                                             //-----element face center
                                                                             element_face_center_3d,
                                                                             //-----vector & matrix -------
                                                                             face_index_direction_extreme_bdry_region,
                                                                             nodes_on_element_boundary_face_line,
                                                                             //------- direction of tangential vector -------
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line,
                                                                             //------- normal direction-------
                                                                             normal_dir,
                                                                             //----------- orientation ------------
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             //--------- matrix of extreeme of control region --------------
                                                                             analitical_integral_extreme_cords);

                c = - analitical_integral_extreme_cords[ normal_tangential_direction_to_bndry_bndry_integration_line ][ first_node_along_element_face_line ] ;

            }
            else if( vector_components_are_equal<double>( nodes_on_element_boundary_face_line[ tangential_face_index_vector[  second_tang_dir ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = tangential_face_index_vector[  second_tang_dir ];
                unsigned int direction_of_bndry_bndry_integration_line = tangential_face_index_vector[  first_tang_dir ];

                a = 0.;
                b = 1.;

                are_nodes_coordinates_decreasing_along_integration_directin = is_vector_oriented_with_cords_decreasing_along_specific_direction(nodes_on_element_boundary_face_line, direction_of_bndry_bndry_integration_line);

                analitical_integral_extreme_cords =integration_extreme_cords(face_index,
                                                                             //dim,
                                                                             //dim_bdry,
                                                                             //-----element face center
                                                                             element_face_center_3d,
                                                                             //-----vector & matrix -------
                                                                             face_index_direction_extreme_bdry_region,
                                                                             nodes_on_element_boundary_face_line,
                                                                             //------- direction of tangential vector -------
                                                                             direction_of_bndry_bndry_integration_line,
                                                                             normal_tangential_direction_to_bndry_bndry_integration_line,
                                                                             //------- normal direction-------
                                                                             normal_dir,
                                                                             //----------- orientation ------------
                                                                             are_nodes_coordinates_decreasing_along_integration_directin,
                                                                             //--------- matrix of extreeme of control region --------------
                                                                             analitical_integral_extreme_cords);

                c = - analitical_integral_extreme_cords[ normal_tangential_direction_to_bndry_bndry_integration_line ][ first_node_along_element_face_line ] ;
            }
      }
//************************* calculus of parameter for analitical solution of mixed integral - END *************************

//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//********************************************* PREVIOUS FUNCTION - END *****************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************


//*************************  parameter for analitical solution of mixed integral - BEGIN *************************
      void coefficent_of_analitical_solution(const std::vector< std::vector< double > > nodes_on_element_boundary_face_line,
                                             //------- i_face qd_point ----------
                                             const std::vector< double> center_of_polar_coords,
                                             //------- tangent vector -------
                                             const std::vector< unsigned int> tangential_face_index_vector,
                                             //------- output -------
                                             double& a, double& b, double& c){

                constexpr unsigned first_tang_dir = 0;
                constexpr unsigned second_tang_dir = 1;

                constexpr unsigned first_node_along_element_face_line = 0;


            if( vector_components_are_equal<double>( nodes_on_element_boundary_face_line[ tangential_face_index_vector[ first_tang_dir ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = tangential_face_index_vector[  first_tang_dir ];// normal is the components (respect the integration line) where points  are equal !!!!!!
                unsigned int direction_of_bndry_bndry_integration_line = tangential_face_index_vector[  second_tang_dir ];

                a = 1.;
                b = 0.;
                c = - ( nodes_on_element_boundary_face_line[ normal_tangential_direction_to_bndry_bndry_integration_line ][ first_node_along_element_face_line ]
                        -center_of_polar_coords[ normal_tangential_direction_to_bndry_bndry_integration_line ]
                       );

            }
            else if( vector_components_are_equal<double>( nodes_on_element_boundary_face_line[ tangential_face_index_vector[  second_tang_dir ] ] ) ) {

                unsigned int normal_tangential_direction_to_bndry_bndry_integration_line = tangential_face_index_vector[  second_tang_dir ];
                unsigned int direction_of_bndry_bndry_integration_line = tangential_face_index_vector[  first_tang_dir ];

                a = 0.;
                b = 1.;
                c = - ( nodes_on_element_boundary_face_line[ normal_tangential_direction_to_bndry_bndry_integration_line ][ first_node_along_element_face_line ]
                         -center_of_polar_coords[ normal_tangential_direction_to_bndry_bndry_integration_line ]
                       );
            }
      }
//*************************  parameter for analitical solution of mixed integral - END *************************






 } //end namespace ctrl

} //end namespace femus


#endif
