
function match_problem(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number

    ! class
    if (type_string == "PROBLEM_ELASTICITY_CLASS") then
      type_number = CMFE_PROBLEM_ELASTICITY_CLASS
    elseif (type_string == "PROBLEM_FLUID_MECHANICS_CLASS") then
      type_number = CMFE_PROBLEM_FLUID_MECHANICS_CLASS
    elseif (type_string == "PROBLEM_ELECTROMAGNETICS_CLASS") then
      type_number = CMFE_PROBLEM_ELECTROMAGNETICS_CLASS
    elseif (type_string == "CLASSICAL_FIELD_CLASS") then
      type_number = CMFE_PROBLEM_CLASSICAL_FIELD_CLASS
    elseif (type_string == "PROBLEM_BIOELECTRICS_CLASS") then
      type_number = CMFE_PROBLEM_BIOELECTRICS_CLASS
    elseif (type_string == "PROBLEM_MODAL_CLASS") then
      type_number = CMFE_PROBLEM_MODAL_CLASS
    elseif (type_string == "PROBLEM_FITTING_CLASS") then
      type_number = CMFE_PROBLEM_FITTING_CLASS
    elseif (type_string == "PROBLEM_OPTIMISATION_CLASS") then
      type_number = CMFE_PROBLEM_OPTIMISATION_CLASS
    elseif (type_string == "PROBLEM_MULTI_PHYSICS_CLASS") then
      type_number = CMFE_PROBLEM_MULTI_PHYSICS_CLASS
    elseif (type_string == "PROBLEM_NO_CLASS") then
      type_number = CMFE_PROBLEM_NO_CLASS
    elseif (type_string == "PROBLEM_FINITE_ELASTICITY_TYPE") then
      type_number = CMFE_PROBLEM_FINITE_ELASTICITY_TYPE
    ! TODO add more types.. TODO
    ! subtype
    elseif (type_string == "U_VARIABLE") then

      type_number = CMFE_FIELD_U_VARIABLE_TYPE
    elseif (type_string == "SEPARATED") then
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    elseif (type_string == "PROBLEM_NO_SUBTYPE") then
      type_number = CMFE_PROBLEM_NO_SUBTYPE

    ! TODO add more subtypes.. TODO
    else
      call handle_error("Invalid string "// type_string)
    endif

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.DEFINES MESH TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function match_generated_mesh(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number

    ! number of dimensions
    if (type_string == "GENERATED_MESH_1D") then
      type_number = 1
    elseif (type_string == "GENERATED_MESH_2D") then
      type_number = 2
    elseif (type_string == "GENERATED_MESH_3D") then
      type_number = 3
    ! generated mesh type
    elseif (type_string == "REGULAR_MESH_TYPE") then
      type_number = CMFE_GENERATED_MESH_REGULAR_MESH_TYPE
     
    elseif (type_string == "POLAR_MESH_TYPE") then
      type_number = CMFE_GENERATED_MESH_POLAR_MESH_TYPE
    elseif (type_string == "TREE_MESH_TYPE") then
      type_number = CMFE_GENERATED_MESH_FRACTAL_TREE_MESH_TYPE
    elseif (type_string == "CYLINDER_MESH_TYPE") then
      type_number = CMFE_GENERATED_MESH_CYLINDER_MESH_TYPE
    elseif (type_string == "ELLIPSOID_MESH_TYPE") then
      type_number = CMFE_GENERATED_MESH_ELLIPSOID_MESH_TYPE
    else
      print *, 2222
    endif

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.REMOVES INTERPOLATION TYPE FOR STATE VARIABLE!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function match_basis(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number
    SELECT CASE (type_string)
    case("LINEAR_LAGRANGE_INTERPOLATION")
    type_number=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    case("QUADRATIC_LAGRANGE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    case("CUBIC_LAGRANGE_INTERPOLATION ")
    type_number=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    case("QUADRATIC1_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC1_HERMITE_INTERPOLATION
    case("CUBIC_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION
    case("QUADRATIC2_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC2_HERMITE_INTERPOLATION
    case("LINEAR_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION
    case("QUADRATIC_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION
    case("CUBIC_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION
    END SELECT
  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.DEFINES COORDINATE SYSTEM !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function match_coordinate_system(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number

    if (type_string == "COORDINATE_RECTANGULAR_CARTESIAN_TYPE") then
      type_number = CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    elseif (type_string == "COORDINATE_CYLINDRICAL_POLAR_TYPE") then
      type_number = CMFE_COORDINATE_CYLINDRICAL_POLAR_TYPE
    elseif (type_string == "COORDINATE_SPHERICAL_POLAR_TYPE") then
      type_number = CMFE_COORDINATE_SPHERICAL_POLAR_TYPE
    elseif (type_string == "COORDINATE_PROLATE_SPHEROIDAL_TYPE") then
      type_number = CMFE_COORDINATE_PROLATE_SPHEROIDAL_TYPE
    elseif (type_string == "COORDINATE_OBLATE_SPHEROIDAL_TYPE") then
      type_number = CMFE_COORDINATE_OBLATE_SPHEROIDAL_TYPE
    elseif (type_string == "COORDINATE_SYSTEM_1D") then
      type_number = 1
    elseif (type_string == "COORDINATE_SYSTEM_2D") then
      type_number = 2
    elseif (type_string == "COORDINATE_SYSTEM_3D") then
      type_number = 3
    else
print *, "i am here"
    print *, "COORDINATE_RECTANGULAR_CARTESIAN_TYPE"
    print *,  type_string
print *, "i am here"

   !    call handle_error("Invalid string "//type_string)
    endif

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.REMOVES INTERPOLATION TYPE FOR STATE VARIABLE!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function match_equations_set(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                   :: type_number
    ! sparsity type
    if (type_string == "SPARSE_MATRIX") then
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    elseif (type_string == "EQUATIONS_FULL_MATRICES") then
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    elseif (type_string == "SEPARATED") then
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ! output type
    elseif (type_string == "NO_OUTPUT") then
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    elseif (type_string == "EQUATIONS_TIMING_OUTPUT") then
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    elseif (type_string == "EQUATIONS_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    elseif (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    elseif (type_string == "U_VARIABLE") then
      type_number = CMFE_FIELD_U_VARIABLE_TYPE
    ! class
    elseif (type_string == "EQUATIONS_SET_ELASTICITY_CLASS") then
      type_number = CMFE_EQUATIONS_SET_ELASTICITY_CLASS
    ! TODO insert more classes.. TODO
    ! type
    elseif (type_string == "EQUATIONS_SET_FINITE_ELASTICITY_TYPE") then
      type_number = CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE
    ! TODO insert more types.. TODO
    ! subtype
    elseif (type_string == "EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE") then
      type_number = CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE
    elseif (type_string == "EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE") then
      type_number = CMFE_EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE
    elseif (type_string == "EQUATIONS_SET_LAPLACE_EQUATION_TYPE") then
      type_number = CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE
    elseif (type_string == "EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE") then
      type_number = CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE
    elseif (type_string == "GUESS_TYPE") then
      type_number = CMFE_FIELD_VALUES_SET_TYPE
    elseif (type_string == "BOUNDARY_CONDITIONS_SET") then
     type_number =  CMFE_FIELD_BOUNDARY_CONDITIONS_SET_TYPE
    elseif (type_string == "UX_VALUE") then
     type_number =  1.0
    elseif (type_string == "UY_VALUE") then
     type_number =  2.0
    elseif (type_string == "UZ_VALUE") then
     type_number =  3.0
    ! TODO insert more subtypes.. TODO
    else
      call handle_error("Invalid string "//type_string)
    endif

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE the OUTPUT!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!OF THE ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   function output_type(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number
    ! sparsity type
    if (type_string == "EQUATIONS_SPARSE_MATRICES") then
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    elseif (type_string == "EQUATIONS_FULL_MATRICES") then
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    elseif (type_string == "SEPARATED") then
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ! output type
    elseif (type_string == "EQUATIONS_NO_OUTPUT") then
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    elseif (type_string == "EQUATIONS_TIMING_OUTPUT") then
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    elseif (type_string == "EQUATIONS_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    elseif (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    else
      call handle_error("Invalid string "//type_string)
    endif
   end function 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE the dependent field!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function match_dependent_field(type_string) result (type_number)
    character(*), intent(in) :: type_string
    integer                  :: type_number
    if (type_string == "FIELD_U_VARIABLE_TYPE") then
    type_number = CMFE_FIELD_U_VARIABLE_TYPE
    end if      
end function match_dependent_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE BCs!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function bc_def(type_string) result (type_number)
    character(*), intent(in) :: type_string
    integer                :: type_number
    if (type_string == "SURFACE") then
    
    elseif (type_string == "RIGHT") then 
    type_number = CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE
    elseif (type_string == "LEFT") then 
    type_number = CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE
    elseif (type_string == "FRONT") then 
    type_number = CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE
    elseif (type_string == "BOTTOM") then 
    type_number = CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE
    end if      
end function 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE SOLVER!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function solver_def(type_string) result (type_number)
    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number
    if (type_string == "SPARSE_MATRIX") then
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    elseif (type_string == "FULL_MATRICES") then
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    elseif (type_string == "SEPARATED") then
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ! output type
    elseif (type_string == "NO_OUTPUT") then
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    elseif (type_string == "EQUATIONS_TIMING_OUTPUT") then
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    elseif (type_string == "EQUATIONS_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    elseif (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") then
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    elseif (type_string == "SOLVER_LINEAR_DIRECT_SOLVE_TYPE") then
      type_number = CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE
    elseif (type_string == "ITERATIVE") then
      type_number = CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
    elseif (type_string == "INCREMENTAL") then
      type_number = CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE
    else
      call handle_error("Invalid string "//type_string)
    endif

  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE SOLVER!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function control_loop_def(type_string) result (type_number)
    implicit none
    character(*), intent(in) :: type_string
    integer                  :: type_number
    if (type_string == "CONTROL_LOAD_INCREMENT_LOOP_TYPE") then
      type_number = CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE
 !!! add more options later
    else
      call handle_error("Invalid string "//type_string)
    endif

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE MATERIALS!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function material_parameters(type_string) result (type_number)

    implicit none
    character(*), intent(in) :: type_string
    integer                 :: type_number
    if (type_string == "EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE") then
      type_number = 5
    else if (type_string == "EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE") then
       type_number =2
    elseif (type_string == "EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE")  then
      type_number = 2
    else
      call handle_error("Invalid string "//type_string)
    endif

  end function material_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING function  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!..used for error handling !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE HANDLE_ERROR(ERROR_STRING)

         CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING
         WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
         STOP

 END SUBROUTINE HANDLE_ERROR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

