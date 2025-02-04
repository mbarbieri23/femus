CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(femus)

#############################################################################################
### Source files
#############################################################################################

SET(femus_src 
00_utils/FemusInit.cpp
00_utils/Files.cpp
00_utils/input_parser/InputParser.cpp
00_utils/input_parser/JsonInputParser.cpp
00_utils/parallel/MyMatrix.cpp
00_utils/parallel/MyVector.cpp
01_problem/MultiLevelProblem.cpp
02_mesh/gencase/Box.cpp
02_mesh/gencase/Domain.cpp
02_mesh/gencase/ElemSto.cpp
02_mesh/gencase/MultiLevelMeshTwo.cpp
02_mesh/gencase/GenCase.cpp
02_mesh/00_geom_elements/GeomElemBase.cpp
02_mesh/00_geom_elements/GeomElemHex8.cpp
02_mesh/00_geom_elements/GeomElemHex27.cpp
02_mesh/00_geom_elements/GeomElemQuad4.cpp
02_mesh/00_geom_elements/GeomElemQuad9.cpp
02_mesh/00_geom_elements/GeomElemTet4.cpp
02_mesh/00_geom_elements/GeomElemTet10.cpp
02_mesh/00_geom_elements/GeomElemTri3.cpp
02_mesh/00_geom_elements/GeomElemTri6.cpp
02_mesh/00_geom_elements/GeomElemTri7.cpp
02_mesh/00_geom_elements/GeomElemEdge2.cpp
02_mesh/00_geom_elements/GeomElemEdge3.cpp
02_mesh/01_input/MeshGeneration.cpp
02_mesh/01_input/GambitIO.cpp
02_mesh/01_input/MED_IO.cpp
02_mesh/02_partitioning/MeshPartitioning.cpp
02_mesh/02_partitioning/MeshMetisPartitioning.cpp
02_mesh/02_partitioning/MeshASMPartitioning.cpp
02_mesh/03_refinement/MeshRefinement.cpp
02_mesh/Elem.cpp
02_mesh/Mesh.cpp
02_mesh/MultiLevelMesh.cpp
03_solution/MultiLevelSolution.cpp
03_solution/Solution.cpp
03_solution/01_output/Writer.cpp
03_solution/01_output/VTKWriter.cpp
03_solution/01_output/GMVWriter.cpp
03_solution/01_output/XDMFWriter.cpp
04_equations/assemble/Quantity.cpp
04_equations/assemble/DofMap.cpp
04_equations/assemble/Assemble_jacobian.cpp
04_equations/assemble/Assemble_unknown.cpp
04_equations/assemble/BoundaryConditions.cpp
04_equations/assemble/CurrentElem.cpp
04_equations/assemble/CurrentQuantity.cpp
04_equations/assemble/SystemTwo.cpp
04_equations/assemble/TimeLoop.cpp
04_equations/time_dependent/TransientSystem.cpp
04_equations/ExplicitSystem.cpp
04_equations/ImplicitSystem.cpp
04_equations/LinearImplicitSystem.cpp
04_equations/MonolithicFSINonLinearImplicitSystem.cpp
04_equations/NonLinearImplicitSystem.cpp
04_equations/System.cpp
algebra/DenseMatrixBase.cpp
algebra/DenseMatrix.cpp
algebra/DenseSubmatrix.cpp
algebra/DenseVectorBase.cpp
algebra/FieldSplitTree.cpp
algebra/Graph.cpp
algebra/LinearEquation.cpp
algebra/LinearEquationSolver.cpp
algebra/NumericVector.cpp
algebra/LinearEquationSolverPetsc.cpp
algebra/LinearEquationSolverPetscAsm.cpp
algebra/LinearEquationSolverPetscFieldSplit.cpp
algebra/PetscMatrix.cpp
algebra/PetscPreconditioner.cpp
algebra/PetscVector.cpp
algebra/Preconditioner.cpp
algebra/SparseMatrix.cpp
algebra/FunctionBase.cpp
algebra/ParsedFunction.cpp
algebra/SlepcSVD.cpp
fe/ElemType.cpp
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
fe_quadrature/ElemType_template.cpp
fe_quadrature/ElemType_template_base.cpp
geom_quadrature/GaussPoints.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
optimal_control/00_cost_functional.cpp
optimal_control/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Solid.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

# IF(HAVE_HDF5)
#   SET(femus_src ${femus_src} 
#   02_mesh/input/MED_IO.cpp
#   03_solution/01_output/XDMFWriter.cpp
#   02_mesh/gencase/GenCase.cpp
#   04_equations/TimeLoop.cpp
#   )
# ENDIF(HAVE_HDF5) 

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femus_src})

IF(APPLE)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libhdf5.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmpi.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.dylib)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libjsoncpp.dylib)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libadept.a)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libb64.a)

#target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${B64_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${JSONCPP_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${ADEPT_LIBRARIES})

ENDIF(APPLE)
