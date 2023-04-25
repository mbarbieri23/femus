#ifndef CONTROL_FACES_SQUARE_OR_CUBE_HPP
#define CONTROL_FACES_SQUARE_OR_CUBE_HPP


namespace femus {

namespace ctrl {


   class Domain_of_PDE {

  public:

static constexpr double domain_length = 1.;
 };
 
 

 namespace square_or_cube {
     

   class List_of_faces : public Domain_of_PDE {

  public:

       static constexpr unsigned _num_of_control_extremes_per_tang_comp_per_face = 2;
       

   static  const  int sign_function_for_delimiting_region(const unsigned int face_index) {

   int  target_line_sign;

        if (face_index == 1 || face_index == 3 || face_index == 5) { target_line_sign = 1;  }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { target_line_sign = -1; }

   return target_line_sign;

   }


 //direction of the line that contains \Gamma_c
      static const unsigned int normal_direction_to_Gamma_control(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;

}



//-------------- my change - BEGIN
     static std::vector<  unsigned int > tangential_direction_to_Gamma_control(const unsigned int face_index, const unsigned boundary_dim) {

         std::vector<  unsigned int > axis_dir( boundary_dim, 0);

         if( boundary_dim == 1){

             if (face_index == 1 || face_index == 2) { axis_dir[ 0 ] = 1; }
        else if (face_index == 3 || face_index == 4) { axis_dir[ 0 ] = 0; }

        }
        else{

              switch(face_index) {

             case(1): { axis_dir[ 0 ] = 2; axis_dir[ 1 ]  = 1; return axis_dir; }
             case(2): { axis_dir[ 0 ] = 1; axis_dir[ 1 ]  = 2; return axis_dir; }

             case(3): { axis_dir[ 0 ] = 0; axis_dir[ 1 ]  = 2; return axis_dir; }
             case(4): { axis_dir[ 0 ] = 2; axis_dir[ 1 ]  = 0; return axis_dir; }

             case(5): { axis_dir[ 0 ] = 1; axis_dir[ 1 ]  = 0; return axis_dir; }
             case(6): { axis_dir[ 0 ] = 0; axis_dir[ 1 ]  = 1; return axis_dir; }

             }

//             if (face_index == 1 || face_index == 2) { axis_dir[ 0 ] = 1; axis_dir[ 1 ]  = 2; }
//        else if (face_index == 3 || face_index == 4) { axis_dir[ 0 ] = 2; axis_dir[ 1 ]  = 0; }
//        else if (face_index == 5 || face_index == 6) { axis_dir[ 0 ] = 0; axis_dir[ 1 ]  = 1; }  ///@todo this depends on the mesh file
//
        }

    return axis_dir;

}
//-------------- my change - END

static bool is_face_on_maximal_coordinate_by_tangential_component(std::vector<  unsigned int > axis_dir){


        if ( ( axis_dir[ 0 ] == 1 && axis_dir[ 1 ]  == 2 ) ||
             ( axis_dir[ 0 ] == 2 && axis_dir[ 1 ]  == 0 ) ||
             ( axis_dir[ 0 ] == 0 && axis_dir[ 1 ]  == 1 )    ) {return true;}

      else { return false; }


}




      static const unsigned int tangential_direction_to_Gamma_control(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 1; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 0; }
   else if (face_index == 5 || face_index == 6) { /*abort();*/ axis_dir = 1; }  ///@todo this depends on the mesh file

    return axis_dir;

}


      static const unsigned int opposite_face(const unsigned int face_index) {

   unsigned int face_opposite = 0;

        if (face_index == 1 || face_index == 3 || face_index == 5) { face_opposite = face_index + 1; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { face_opposite = face_index - 1; }

    return face_opposite;

}


      static const double normal_coordinate(const bool is_face_on_maximal_coordinate) {

   const double square_or_cube__control_face_normal_coord_min = 0.;
   const double square_or_cube__control_face_normal_coord_max = domain_length;
   
       return (is_face_on_maximal_coordinate)? square_or_cube__control_face_normal_coord_max : square_or_cube__control_face_normal_coord_min;

      }
      

// static std::vector< double>  exteeme_value_along_some_direction (){}



//-------------- create_matrix_of_cords_of_bndry_bndry_points_from_face_index- BEGIN --------------

     static std::vector< std::vector< double > > create_matrix_of_cords_of_bndry_bndry_points_from_face_index(const unsigned int face_index,
                                                                                                              const unsigned boundary_dim,
                                                                                                              const unsigned numbers_of_nodes,
                                                                                                              const double first_direction_value_one,
                                                                                                              const double first_direction_value_two,
                                                                                                              const double second_direction_value_one,
                                                                                                              const double second_direction_value_two){



         // vectors of extreeme along tangential direction - BEGIN
         std::vector< double > exteeme_value_along_first_tangential_direction( boundary_dim);
         exteeme_value_along_first_tangential_direction[0]  = first_direction_value_one;
         exteeme_value_along_first_tangential_direction[1]  = first_direction_value_two;

         std::vector< double > exteeme_value_along_second_tangential_direction( boundary_dim);
         exteeme_value_along_second_tangential_direction[0] = second_direction_value_one;
         exteeme_value_along_second_tangential_direction[1] = second_direction_value_two;
         // vectors of extreeme along tangential direction - END

         std::vector< unsigned int > tangetntial_vector = tangential_direction_to_Gamma_control( face_index, boundary_dim );
         bool is_face_on_maximal_coordinate = is_face_on_maximal_coordinate_by_tangential_component( tangetntial_vector );
         unsigned first_tangential_direction = tangetntial_vector[0];
         unsigned second_tangential_direction = tangetntial_vector[1];

         //********************* Matrix of nodes of the bndry bndry - BEGIN *********************

         //creation of the matrix - BEGIN
         std::vector< std::vector< double > > bndry_bdnry_matrix_of_cords(boundary_dim + 1);
         for(unsigned f = 0; f <  bndry_bdnry_matrix_of_cords.size(); f++) {
             bndry_bdnry_matrix_of_cords[f].resize(numbers_of_nodes);
         }
         //creation of the matrix - END

         // filling the matrix of bndry bndry cords along normal direction - BEGIN
         unsigned int normal_dir = normal_direction_to_Gamma_control( face_index );
         for (int p = 0; p < bndry_bdnry_matrix_of_cords[ normal_dir ].size(); p++){
             bndry_bdnry_matrix_of_cords[ normal_dir ][ p ] = normal_coordinate(is_face_on_maximal_coordinate);
         }
         // filling the matrix of bndry bndry cords along normal direction - END

         // filling the matrix of bndry bndry cords along tangential direction - BEGIN
         for (int t = 0; t < bndry_bdnry_matrix_of_cords[ normal_dir ].size(); t++){
             bndry_bdnry_matrix_of_cords[ first_tangential_direction ][ t ]  = exteeme_value_along_first_tangential_direction[ ( ( ( 3 - t )*( t % 3 ) ) / 2 )];
             bndry_bdnry_matrix_of_cords[ second_tangential_direction ][ t ] = exteeme_value_along_second_tangential_direction[ ( t / 2 )];
         }
         // filling the matrix of bndry bndry cords along tangential direction - END

         //********************* Matrix of nodes of the bndry bndry - END *********************

//          bndry_bdnry_matrix_of_cords[first_tangential_direction][0]  = first_direction_value_one;
//          bndry_bdnry_matrix_of_cords[second_tangential_direction][0] = second_direction_value_one;
//
//          bndry_bdnry_matrix_of_cords[first_tangential_direction][1]  = first_direction_value_two;
//          bndry_bdnry_matrix_of_cords[second_tangential_direction][1] = second_direction_value_one;
//
//          bndry_bdnry_matrix_of_cords[first_tangential_direction][2]  = first_direction_value_two;
//          bndry_bdnry_matrix_of_cords[second_tangential_direction][2] = second_direction_value_two;
//
//          bndry_bdnry_matrix_of_cords[first_tangential_direction][3]  = first_direction_value_one;
//          bndry_bdnry_matrix_of_cords[second_tangential_direction][3] = second_direction_value_two;

         return bndry_bdnry_matrix_of_cords;

    }

//-------------- create_matrix_of_cords_of_bndry_bndry_points_from_face_index- END --------------

 };
 
// Here, the assumption is that each face of the Domain contains at most 1 interval of Control (hence, two extremes)
template < unsigned N_TANG_COMPONENTS >
      class List_of_Gamma_control_faces : public square_or_cube::List_of_faces {

  public:

       static constexpr unsigned  _num_of_tang_components_per_face  = N_TANG_COMPONENTS;
       
       static  std::map<unsigned int, unsigned int >  from_ctrl_faces_to_their_boundaries() {
          
         // ctrl face to node-node association - BEGIN
          std::map<unsigned int, unsigned int >  ctrl_faces_VS_their_nodes;
          
          ctrl_faces_VS_their_nodes[1] =  7  ;
          ctrl_faces_VS_their_nodes[2] =  8  ;
          ctrl_faces_VS_their_nodes[3] =  9  ;
          ctrl_faces_VS_their_nodes[4] =  10  ;
          ctrl_faces_VS_their_nodes[5] =  11  ;
          ctrl_faces_VS_their_nodes[6] =  12  ;
         // ctrl face to node-node association - END
 
 return ctrl_faces_VS_their_nodes;

      }
       

   };

 
 }
 
 
namespace boundary_control_between_extreme {

    
  namespace square {



//*********************** Single - BEGIN *****************************************


  class List_of_Gamma_control_faces_Single : public square_or_cube::List_of_Gamma_control_faces<1> {

  public:

     static constexpr unsigned _face_with_extremes_index_size = 1 ;

     static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
         true
    };

     static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


  };

      const double   List_of_Gamma_control_faces_Single::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }
    };



  class List_of_Gamma_control_faces_One :  public  List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

     const unsigned List_of_Gamma_control_faces_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    1    };



  class List_of_Gamma_control_faces_Two : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    2    };


  class List_of_Gamma_control_faces_Three : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    3    };



     class List_of_Gamma_control_faces_Four : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    4    };

//*********************** Single - END *****************************************


//*********************** Double - BEGIN *****************************************


  class List_of_Gamma_control_faces_Double : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 2 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Double::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = { 
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  } },
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  } }
      };

      class List_of_Gamma_control_faces_One_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  2   };

    class List_of_Gamma_control_faces_One_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  3   };




     class List_of_Gamma_control_faces_One_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  4   };


    class List_of_Gamma_control_faces_Two_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  1   };


    class List_of_Gamma_control_faces_Two_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  3   };



    class List_of_Gamma_control_faces_Two_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  4   };




    class List_of_Gamma_control_faces_Three_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  1   };



    class List_of_Gamma_control_faces_Three_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  2   };





    class List_of_Gamma_control_faces_Three_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  4   };




    class List_of_Gamma_control_faces_Four_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  1   };



    class List_of_Gamma_control_faces_Four_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  2   };





    class List_of_Gamma_control_faces_Four_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  3   };




//*********************** Double - END *****************************************


//*********************** Triple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Triple : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 3 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Triple::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } },
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } },
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }

      };


      class List_of_Gamma_control_faces_One_Three_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    2  };



      class List_of_Gamma_control_faces_One_Four_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    2  };



      class List_of_Gamma_control_faces_One_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    4  };



      class List_of_Gamma_control_faces_One_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    3  };





      class List_of_Gamma_control_faces_Two_Three_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    1  };



      class List_of_Gamma_control_faces_Two_Four_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    1  };



      class List_of_Gamma_control_faces_Two_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    4  };


      class List_of_Gamma_control_faces_Two_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    3  };



      class List_of_Gamma_control_faces_Three_One_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    4  };



      class List_of_Gamma_control_faces_Three_Two_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    4  };


      class List_of_Gamma_control_faces_Three_One_Two:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    2  };


      class List_of_Gamma_control_faces_Three_Two_One:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    1  };


      class List_of_Gamma_control_faces_Four_One_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    3  };



      class List_of_Gamma_control_faces_Four_Two_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  2  ,    3  };



      class List_of_Gamma_control_faces_Four_One_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    2  };

//*********************** Triple - END *****************************************


//*********************** Quadruple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Quadruple : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 4 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Quadruple::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
        { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }
       ,{ { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }
       ,{ { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }
       ,{ { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX } }
      };




      class List_of_Gamma_control_faces_One_Three_Two_Four :  public  List_of_Gamma_control_faces_Quadruple {

        public:

           static const unsigned _face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two_Four::_face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ] = {    1  ,  3  , 4  ,   2  };



//*********************** Quadruple - END *****************************************


} // end namespace square



  namespace cube {
      


//*********************** Single - BEGIN *****************************************


  class List_of_Gamma_control_faces_Single : public square_or_cube::List_of_Gamma_control_faces<2> {

  public:

     static constexpr unsigned _face_with_extremes_index_size = 1 ;

     static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
         true
    };

     static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


  };


  const double   List_of_Gamma_control_faces_Single::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX },
            { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX /*0., 1.*/ } }
    };




    class List_of_Gamma_control_faces_One :  public  List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

     const unsigned List_of_Gamma_control_faces_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    1    };



  class List_of_Gamma_control_faces_Two : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    2    };


      class List_of_Gamma_control_faces_Three : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    3    };



      class List_of_Gamma_control_faces_Four : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    4    };
      

      
     class List_of_Gamma_control_faces_Five : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Five::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    5    };
      


     class List_of_Gamma_control_faces_Six : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Six::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    6    };
      

//*********************** Single - END *****************************************


//*********************** Double - BEGIN *****************************************


  class List_of_Gamma_control_faces_Double : public square_or_cube::List_of_Gamma_control_faces<2> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 2 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][_num_of_tang_components_per_face][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Double::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][_num_of_tang_components_per_face][_num_of_control_extremes_per_tang_comp_per_face] = {
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  },
            { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  } },
          { { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  },
            { SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN, SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX  } }
      };

      class List_of_Gamma_control_faces_One_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  2   };

    class List_of_Gamma_control_faces_One_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  3   };




     class List_of_Gamma_control_faces_One_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  4   };


    class List_of_Gamma_control_faces_Two_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  1   };


    class List_of_Gamma_control_faces_Two_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  3   };



    class List_of_Gamma_control_faces_Two_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  4   };




    class List_of_Gamma_control_faces_Three_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  1   };



    class List_of_Gamma_control_faces_Three_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  2   };





    class List_of_Gamma_control_faces_Three_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  4   };




    class List_of_Gamma_control_faces_Four_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  1   };



    class List_of_Gamma_control_faces_Four_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  2   };





    class List_of_Gamma_control_faces_Four_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  3   };




//*********************** Double - END *****************************************



    } // end namespace cube

      
} //end namespace boundary_control_between_extreme

 
 namespace  boundary_control_full_face {
     
  namespace square {



//*********************** Single - BEGIN *****************************************


  class List_of_Gamma_control_faces_Single : public square_or_cube::List_of_Gamma_control_faces<1> {

  public:

     static constexpr unsigned _face_with_extremes_index_size = 1 ;

     static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
         true
    };

     static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


  };

      const double   List_of_Gamma_control_faces_Single::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
       { { 0, domain_length } }
    };



  class List_of_Gamma_control_faces_One :  public  List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

     const unsigned List_of_Gamma_control_faces_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    1    };



  class List_of_Gamma_control_faces_Two : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    2    };


  class List_of_Gamma_control_faces_Three : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    3    };



     class List_of_Gamma_control_faces_Four : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    4    };

//*********************** Single - END *****************************************


//*********************** Double - BEGIN *****************************************


  class List_of_Gamma_control_faces_Double : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 2 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Double::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
         { { 0 , domain_length  } }
       , { { 0 , domain_length  } }

      };

      class List_of_Gamma_control_faces_One_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  2   };

    class List_of_Gamma_control_faces_One_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  3   };




     class List_of_Gamma_control_faces_One_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  4   };


    class List_of_Gamma_control_faces_Two_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  1   };


    class List_of_Gamma_control_faces_Two_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  3   };



    class List_of_Gamma_control_faces_Two_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  4   };




    class List_of_Gamma_control_faces_Three_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  1   };



    class List_of_Gamma_control_faces_Three_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  2   };





    class List_of_Gamma_control_faces_Three_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  4   };




    class List_of_Gamma_control_faces_Four_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  1   };



    class List_of_Gamma_control_faces_Four_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  2   };





    class List_of_Gamma_control_faces_Four_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  3   };




//*********************** Double - END *****************************************


//*********************** Triple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Triple : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 3 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Triple::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
         { { 0, domain_length } } 
       , { { 0, domain_length } }
       , { { 0, domain_length } }

      };


      class List_of_Gamma_control_faces_One_Three_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    2  };



      class List_of_Gamma_control_faces_One_Four_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    2  };



      class List_of_Gamma_control_faces_One_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    4  };



      class List_of_Gamma_control_faces_One_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    3  };





      class List_of_Gamma_control_faces_Two_Three_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    1  };



      class List_of_Gamma_control_faces_Two_Four_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    1  };



      class List_of_Gamma_control_faces_Two_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    4  };


      class List_of_Gamma_control_faces_Two_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    3  };



      class List_of_Gamma_control_faces_Three_One_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    4  };



      class List_of_Gamma_control_faces_Three_Two_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    4  };


      class List_of_Gamma_control_faces_Three_One_Two:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    2  };


      class List_of_Gamma_control_faces_Three_Two_One:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    1  };


      class List_of_Gamma_control_faces_Four_One_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    3  };



      class List_of_Gamma_control_faces_Four_Two_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  2  ,    3  };



      class List_of_Gamma_control_faces_Four_One_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    2  };

//*********************** Triple - END *****************************************


//*********************** Quadruple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Quadruple : public square_or_cube::List_of_Gamma_control_faces<1> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 4 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Quadruple::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face] = {
         { { 0, domain_length } }
       , { { 0, domain_length } }
       , { { 0, domain_length } }
       , { { 0, domain_length } }

      };




      class List_of_Gamma_control_faces_One_Three_Two_Four :  public  List_of_Gamma_control_faces_Quadruple {

        public:

           static const unsigned _face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two_Four::_face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ] = {    1  ,  3  , 4  ,   2  };



//*********************** Quadruple - END *****************************************



} // end namespace square



 

  namespace cube {

//*********************** Single - BEGIN *****************************************

  class List_of_Gamma_control_faces_Single : public square_or_cube::List_of_Gamma_control_faces<2> {

  public:

     static constexpr unsigned _face_with_extremes_index_size = 1 ;

     static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
         true
    };

     static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ][_num_of_control_extremes_per_tang_comp_per_face];


  };

      const double   List_of_Gamma_control_faces_Single::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][ _num_of_tang_components_per_face ]  [_num_of_control_extremes_per_tang_comp_per_face] = {
       { { 0, domain_length },
         { 0, domain_length } }
    };



  class List_of_Gamma_control_faces_One :  public  List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

     const unsigned List_of_Gamma_control_faces_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    1    };



  class List_of_Gamma_control_faces_Two : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    2    };


  class List_of_Gamma_control_faces_Three : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    3    };



     class List_of_Gamma_control_faces_Four : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    4    };


      
     class List_of_Gamma_control_faces_Five : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Five::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    5    };
      


     class List_of_Gamma_control_faces_Six : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Six::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    6    };
      

//*********************** Single - END *****************************************


//*********************** Double - BEGIN *****************************************


  class List_of_Gamma_control_faces_Double : public square_or_cube::List_of_Gamma_control_faces<2> {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 2 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
      };

       static const double   _face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][_num_of_tang_components_per_face][_num_of_control_extremes_per_tang_comp_per_face];


    };


      const double   List_of_Gamma_control_faces_Double::_face_with_extremes_extremes_on_tang_surface[ _face_with_extremes_index_size ][_num_of_tang_components_per_face][_num_of_control_extremes_per_tang_comp_per_face] = {
         { { 0 , domain_length  }, { 0, domain_length } }
       , { { 0 , domain_length  }, { 0, domain_length } }

      };

      class List_of_Gamma_control_faces_One_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  2   };

    class List_of_Gamma_control_faces_One_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  3   };




     class List_of_Gamma_control_faces_One_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  4   };


    class List_of_Gamma_control_faces_Two_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  1   };


    class List_of_Gamma_control_faces_Two_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  3   };



    class List_of_Gamma_control_faces_Two_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  4   };




    class List_of_Gamma_control_faces_Three_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  1   };



    class List_of_Gamma_control_faces_Three_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  2   };





    class List_of_Gamma_control_faces_Three_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  4   };




    class List_of_Gamma_control_faces_Four_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  1   };



    class List_of_Gamma_control_faces_Four_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  2   };





    class List_of_Gamma_control_faces_Four_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  3   };




//*********************** Double - END *****************************************





     } // end namespace cube

  } // end  namespace  boundary_control_full_face
  
  
  
 
  }// end namespace ctrl

} // end namespace femus



#endif

