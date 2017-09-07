!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Monodomain equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example Bioelectrics/Monodomain/src/MonodomainExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS calls. Example-0401
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM MONODOMAINEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP) :: PhysicalStimulationLength =0.03125_CMISSRP !0.0016_CMISSRP   ! X-direction   ### PAPERBRANCH SETTING: a value small enough, such that ONLY ONE CELL is stimulated. !NMJ area: 200 (um)² -> NMJ diameter: 16 um = 0.0016cm. Based on Tse et al., 2014, The Neuromuscular Junction: Measuring Synapse Size, Fragmentation and Changes in Synaptic Protein Density Using Confocal Fluorescence Microscopy

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=25500) :: COMMAND_ARGUMENT,CellmlFile,Filename,SplittingOrder
  LOGICAL :: fileExist

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS
  INTEGER(CMISSIntg) :: SOLVER_TYPE, INTERPOLATION_TYPE, NUMBER_OF_GAUSS_XI
  LOGICAL :: SLOW_TWITCH

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: CellMLModelIndex

  INTEGER(CMISSIntg) :: gNacomponent,StimComponent,node_idx

  REAL(CMISSRP) :: X,Y,DISTANCE,gNa_VALUE
  
  INTEGER(CMISSIntg) :: OUTPUT_FREQUENCY = 1
  REAL(CMISSRP) :: STIM_VALUE
  REAL(CMISSRP), PARAMETER :: STIM_STOP = 0.10_CMISSRP
  REAL(CMISSRP) :: TIME_STOP = 10.0_CMISSRP
  REAL(CMISSRP) :: ODE_TIME_STEP
  REAL(CMISSRP) :: PDE_TIME_STEP
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY = 3.828_CMISSRP

  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: JunctionNodeIdx,NumberStimulatedNodesPerFibre,StimulatedNodeBegin,StimulatedNodeEnd,StimulationNodeIdx
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
  INTEGER(CMISSIntg) :: Err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


  ! process command line arguments before getting started.
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 3) THEN
    ! first argument: X elements
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    WRITE(*, '("Number X elements: ", I0)') NUMBER_GLOBAL_X_ELEMENTS

    ! 2rd argument: INTERPOLATION_TYPE
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid interpolation specification.")
    
    ! 3th argument: solver type
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SOLVER_TYPE
    IF((SOLVER_TYPE<0).OR.(SOLVER_TYPE>1)) CALL HANDLE_ERROR("Invalid solver type specification.")

    ! 4th argument pde time step
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) PDE_TIME_STEP
    WRITE(*, '("PDE Step Size: ", E14.7)') PDE_TIME_STEP

    ! 5th argument: time stop
    CALL GET_COMMAND_ARGUMENT(5,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) TIME_STOP
    WRITE(*, '("Stop Time: ", E14.7)') TIME_STOP

    ! 6th argument: output frequency
    CALL GET_COMMAND_ARGUMENT(6,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) OUTPUT_FREQUENCY
    WRITE(*, '("Output Frequency: ", I10)') OUTPUT_FREQUENCY

    ! 7th argument cellml file
    CALL GET_COMMAND_ARGUMENT(7,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    CellmlFile = adjustl(COMMAND_ARGUMENT)
    WRITE(*, '("CellML File: ", A)') TRIM(CellmlFile)
    
    ! 8th argument fast or slow twitch parameters of cellml model
    CALL GET_COMMAND_ARGUMENT(8,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SLOW_TWITCH
    WRITE(*, '("Slow Twitch: ", L)') SLOW_TWITCH
    
    ! 9th argument ode time step
    CALL GET_COMMAND_ARGUMENT(9,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 9.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) ODE_TIME_STEP
    WRITE(*, '("ODE Step Size: ", E14.7)') ODE_TIME_STEP
    
    ! 10th argument ode time step
    CALL GET_COMMAND_ARGUMENT(10,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    SplittingOrder=adjustl(COMMAND_ARGUMENT)
    WRITE(*, '("Splitting order: ", A)') TRIM(SplittingOrder)

    inquire(file=CellmlFile, exist=fileExist)
    if (.not. fileExist) then
      write(*, '(">>ERROR: File does not exist")')
      stop
    endif
  ELSE
    !If there are not enough arguments default the problem specification 
    WRITE(*,'(">>USAGE: ",A)') "MonodomainExample <number elements X> <interpolation type> <solver type>" // &
      & "<PDE step size> <stop time> <output frequency> <CellML Model URL>"
    
    NUMBER_GLOBAL_X_ELEMENTS=16
    SOLVER_TYPE=0
    !INTERPOLATION_TYPE=1
    
    INTERPOLATION_TYPE=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    !INTERPOLATION_TYPE=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    !INTERPOLATION_TYPE=CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION 

    ! BASIS_LINEAR_LAGRANGE_INTERPOLATION=1 !<Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_QUADRATIC_LAGRANGE_INTERPOLATION=2 !<Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_CUBIC_LAGRANGE_INTERPOLATION=3 !<Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_CUBIC_HERMITE_INTERPOLATION=4 !<Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_QUADRATIC1_HERMITE_INTERPOLATION=5 !<Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_QUADRATIC2_HERMITE_INTERPOLATION=6 !<Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_LINEAR_SIMPLEX_INTERPOLATION=7 !<Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_QUADRATIC_SIMPLEX_INTERPOLATION=8 !<Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    ! BASIS_CUBIC_SIMPLEX_INTERPOLATION=9 !<Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
    
    PDE_TIME_STEP = 0.0005_CMISSRP
    ODE_TIME_STEP = 0.0001_CMISSRP
    
    TIME_STOP=10.00
    OUTPUT_FREQUENCY=1
    CellmlFile="slow_TK_2014_12_08.xml"
    !CellmlFile="hodgkin_huxley_1952.cellml"   
    SLOW_TWITCH=.TRUE.
  ENDIF

  ! determine file name for output files
  WRITE(filename, "(A,I0,A,I1,A,I0,A)") &
    & "results/current_run/l1x1_n", &
    & NUMBER_GLOBAL_X_ELEMENTS, &
    & "_i",INTERPOLATION_TYPE,"_s",SOLVER_TYPE,"/Example"
  filename=trim(filename)
  PRINT *, "Filename: """ // TRIM(Filename) // """"
  
  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
!  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)

  WRITE(Filename,'(A,"_",I0,"x",I0)') "Monodomain",NUMBER_GLOBAL_X_ELEMENTS
  
  CALL cmfe_OutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
  
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Set the region label
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is linear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  
  !Set the basis to be a linear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(Basis,1,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE],Err)
  IF(NUMBER_OF_GAUSS_XI>0) THEN
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI],Err)
  ENDIF
  
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
        
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  !Set the equations set to be a Monodomain equations set
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,&
    & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS,CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],&
    & EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 500.0_CMISSRP,Err)
  
  IF(CellmlFile .EQ. "hodgkin_huxley_1952.cellml") THEN
    IF(SLOW_TWITCH) THEN
      !Set Cm, slow-twitch
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & 0.58_CMISSRP,Err)

      IF(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE > 8) THEN
          STIM_VALUE=(225.0_CMISSRP/8.0_CMISSRP)*(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE)
      ELSE
          STIM_VALUE=225.0_CMISSRP
      ENDIF
    ELSE  
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & 1.0_CMISSRP,Err)
      IF(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE > 16) THEN
        STIM_VALUE=(375.0_CMISSRP/16.0_CMISSRP)*(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE)
      ELSE
        STIM_VALUE=375.0_CMISSRP
      ENDIF
    ENDIF    
  ELSEIF(CellmlFile .EQ. "slow_TK_2014_12_08.xml") THEN
    IF(SLOW_TWITCH) THEN
      !Set Cm, slow-twitch
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & 0.58_CMISSRP,Err)
      !!!IF(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE > 10) THEN
      !!!    STIM_VALUE=(375.0_CMISSRP/10.0_CMISSRP)*(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE)
      !!!ELSE
      !!!    STIM_VALUE=375.0_CMISSRP
      !!!ENDIF
      STIM_VALUE=1200.0_CMISSRP
    ELSE  
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & 1.0_CMISSRP,Err)
      IF(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE > 12) THEN
        STIM_VALUE=(375.0_CMISSRP/12.0_CMISSRP)*(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE)
      ELSE
        STIM_VALUE=375.0_CMISSRP
      ENDIF
    ENDIF
  ENDIF
  
  !Set conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & CONDUCTIVITY,Err)

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL cmfe_CellML_ModelImport(CellML,CellmlFile,CellMLModelIndex,Err)
   
  !---------------------------------------------------------------------------------------------------------------------------------
  ! To Include Shorten
  IF(CellmlFile .EQ. "hodgkin_huxley_1952.cellml") THEN
    CALL cmfe_CellML_VariableSetAsKnown(CellML,CellMLModelIndex,"membrane/i_Stim ",Err)
    CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_Na",Err)
  ELSEIF(CellmlFile .EQ. "slow_TK_2014_12_08.xml") THEN
    CALL cmfe_CellML_VariableSetAsKnown(CellML,CellMLModelIndex,"wal_environment/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"wal_environment/I_Na",Err)
  ENDIF
  !---------------------------------------------------------------------------------------------------------------------------------
  
  !CALL cmfe_CellML_VariableSetAsKnown(CellML,CellMLModelIndex,"fast_sodium_current/g_Na ",Err)
  !CALL cmfe_CellML_VariableSetAsKnown(CellML,CellMLModelIndex,"membrane/IStim",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_K1",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_to",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_K",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_K_ATP",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_Ca_L_K",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_b_K",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_NaK",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_Na",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_b_Na",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_Ca_L_Na",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"membrane/i_NaCa",Err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  IF(CellmlFile .EQ. "hodgkin_huxley_1952.cellml") THEN
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
      & CellMLModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,CellMLModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE, &
      & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  ELSEIF(CellmlFile .EQ. "slow_TK_2014_12_08.xml") THEN
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
      & CellMLModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,CellMLModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
      & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)  
  ENDIF
  
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !Set initial value for Vm 
  IF(CellmlFile .EQ. "hodgkin_huxley_1952.cellml") THEN
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & -75.0_CMISSRP,Err)
  ELSEIF(CellmlFile .EQ. "slow_TK_2014_12_08.xml") THEN
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & -79.974_CMISSRP,Err)
  ENDIF
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the CellML models field
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Find the domains of the first and last nodes
  FirstNodeNumber=1
  LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE+1)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  
  IF(CellmlFile .EQ. "hodgkin_huxley_1952.cellml") THEN
    CALL cmfe_CellML_FieldComponentGet(CellML,CellMLModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"membrane/i_Stim",StimComponent,Err)
  ELSEIF(CellmlFile .EQ. "slow_TK_2014_12_08.xml") THEN
    CALL cmfe_CellML_FieldComponentGet(CellML,CellMLModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"wal_environment/I_HH", &
      & StimComponent,Err)
  ENDIF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! compute number of bioelectric nodes that will be stimulated
  !---------------------------------------------------------------------------------------------------------------------------------
  NumberStimulatedNodesPerFibre = MAX(1, NINT(DBLE(PhysicalStimulationLength) * &
                                                 &((NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE+1) / WIDTH)))
  !PRINT *, "Number of Nodes which are stimulated per fibre: ", NumberStimulatedNodesPerFibre
  
  ! get middle point of fibre
  JunctionNodeIdx = INT(CEILING(DBLE(NUMBER_GLOBAL_X_ELEMENTS*INTERPOLATION_TYPE+1)/2))
  !PRINT *, "JunctionNodeIdx: ", JunctionNodeIdx
  
  ! compute first node for stimulation
  StimulatedNodeBegin = JunctionNodeIdx - NumberStimulatedNodesPerFibre/2
  !PRINT *, "StimulatedNodeBegin: ", StimulatedNodeBegin
  !compute first node for stimulation
  StimulatedNodeEnd = StimulatedNodeBegin + NumberStimulatedNodesPerFibre-1
  !PRINT *, "StimulatedNodeEnd: ", StimulatedNodeEnd

  ! loop over nodes of the fibre to be stimulated
  DO StimulationNodeIdx = StimulatedNodeBegin, StimulatedNodeEnd 
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,StimulationNodeIdx,1,NodeDomain,Err)
    !turn stimulus on at the stimulation point
    IF(NodeDomain==ComputationalNodeNumber) THEN
       CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
            & StimulationNodeIdx,StimComponent,STIM_VALUE,Err)
    ENDIF
  ENDDO
  !---------------------------------------------------------------------------------------------------------------------------------
  !Set up the g_Na gradient
  !CALL cmfe_CellML_FieldComponentGet(CellML,CellMLModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
  !  & gNacomponent,Err)
  !Loop over the nodes
  !DO node_idx=1,LastNodeNumber
  !  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
  !  IF(NodeDomain==ComputationalNodeNumber) THEN
  !    CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
  !      & X,Err)
  !    CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
  !      & Y,Err)
  !    DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSRP)
  !    gNa_VALUE=2.0_CMISSRP*(DISTANCE+0.5_CMISSRP)*385.5e-3_CMISSRP
  !    CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
  !      & node_idx,gNacomponent,gNa_VALUE,Err)
  !  ENDIF
  !ENDDO
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a problem.
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Problem_Initialise(Problem,Err)
  
  
  IF(SplittingOrder .EQ. "O1") THEN
    CALL cmfe_Problem_CreateStart(ProblemUserNumber,&
    & [CMFE_PROBLEM_BIOELECTRICS_CLASS,CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE,CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE],&
    & Problem,Err)
    PRINT *,"Splitting: ", CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE
  ELSEIF(SplittingOrder.EQ. "O2") THEN
    CALL cmfe_Problem_CreateStart(ProblemUserNumber,&
    & [CMFE_PROBLEM_BIOELECTRICS_CLASS,CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE,&
    & CMFE_PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE],Problem,Err)
    PRINT *, "Splitting: ",CMFE_PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE
  ENDIF
    
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

   
  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,STIM_STOP,PDE_TIME_STEP,Err)
  !Set the output
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
  !Set the output frequency (0 for no output, n for output every n time steps)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solvers
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Get the first (DAE) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  !Set the DAE time step
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err) 
  !changing to the Heun method for the Strang-Splitting
  IF(SplittingOrder .EQ. "O2") THEN
    CALL cmfe_Solver_DAEEulerSolverTypeSet(Solver,CMFE_SOLVER_DAE_EULER_IMPROVED,Err)
  ENDIF
  
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  !Get the second (Parabolic) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  
  IF(SplittingOrder .EQ. "O2") THEN
    !!! changing to CRANK_NICOLSON
    CALL cmfe_Solver_DynamicSchemeSet(Solver,CMFE_SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,Err)
  ENDIF
  
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  !!! The third solver is created here(!)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

! TODO: this was copied from laplace example
!  IF(SOLVER_TYPE==0) THEN
!    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!    CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
!  ELSE
!    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
!    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!  ENDIF
  
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solver CellML equations
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  IF(SplittingOrder .EQ. "O2") THEN
    !Get the first solver  
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
    !Get the CellML equations
    CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
    CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
    !Add in the CellML environement
    CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)  
  ENDIF
  !---------------------------------------------------------------------------------------------------------------------------------
  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solver parabolic equations
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the first STIM_STOP
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_Problem_Solve(Problem,Err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  ! StimulatedNodeBegin, StimulatedNodeEnd are set previously
  DO StimulationNodeIdx = StimulatedNodeBegin, StimulatedNodeEnd
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,StimulationNodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & StimulationNodeIdx,StimComponent,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !---------------------------------------------------------------------------------------------------------------------------------
  !Set the time loop from STIM_STOP to TIME_STOP
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)  
  !Solve the problem for the rest time
  CALL cmfe_Problem_Solve(Problem,Err)
  
  !Export results
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM MONODOMAINEXAMPLE
