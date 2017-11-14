!> \file
!> \author Adam Reeve, Benjamin Maier, Nehzat Emamy, Aaron Krämer, Thomas Klotz
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Contributor(s): Adam Reeve, Thomas Heidlauf
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

!> \example FiniteElasticity/LargeUniAxialExtension/src/LargeUniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE

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


  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
  LOGICAL :: DebuggingOnlyRunShortPartOfSimulation = .FALSE.    ! only run one timestep of MAIN_LOOP with stimulus
  LOGICAL, PARAMETER :: DEBUGGING_PARALLEL_BARRIER = .FALSE.   ! execute a barrier where all processes wait when running in parallel, this is helpful to only debug the execution of single processes in a parallel scenario with gdb
  LOGICAL, PARAMETER :: DEBUGGING_PROBLEM_OUTPUT = .FALSE.     ! output the 'solver' object after it is created
  LOGICAL :: DebuggingOutput = .FALSE.    ! enable information from solvers
  INTEGER(CMISSIntg) :: ModelType = 0 ! ### PAPERBRANCH SETTING     ! type of the model (was OldTomoMechanics): 0 = "3a","MultiPhysStrain", old version of tomo that works in parallel, 1 = "3","MultiPhysStrain", new version of tomo that is more stable in numerical sense, 2 = "4","Titin"
  INTEGER(CMISSIntg) :: SplittingType = 0   ! 0 = godunov splitting, 1 = strang splitting
  LOGICAL :: EnableExportEMG = .FALSE.
  LOGICAL :: ElasticityDisabled = .FALSE. ! do create the elasticity control loop
  
  ! physical dimensions in [cm]
  REAL(CMISSRP) :: PhysicalLength=1.0_CMISSRP ! ### PAPERBRANCH SETTING !    X-direction
  REAL(CMISSRP) :: PhysicalWidth =1.0_CMISSRP ! ### PAPERBRANCH SETTING !    Y-direction
  REAL(CMISSRP) :: PhysicalHeight=1.0_CMISSRP ! ### PAPERBRANCH SETTING !    Z-direction
  REAL(CMISSRP) :: PhysicalStimulationLength = 0.03125_CMISSRP  ! X-direction   ### PAPERBRANCH SETTING: value 0.03125 is unphysical (roughly 20 times too large). !NMJ area: 200 (um)² -> NMJ diameter: 16 um = 0.0016cm. Based on Tse et al., 2014, The Neuromuscular Junction: Measuring Synapse Size, Fragmentation and Changes in Synaptic Protein Density Using Confocal Fluorescence Microscopy
  
  !all times in [ms]
  REAL(CMISSRP) :: time
  REAL(CMISSRP) :: StimPeriod=10.0_CMISSRP ! ### PAPERBRANCH SETTING
  REAL(CMISSRP) :: TimeStop=10.0_CMISSRP ! ### PAPERBRANCH SETTING

  ! time discretizatzion settings
  REAL(CMISSRP) :: ODETimeStep = 0.0001_CMISSRP ! ### PAPERBRANCH SETTING: 0.0001
  REAL(CMISSRP) :: PDETimeStep = 0.0005_CMISSRP ! ### PAPERBRANCH SETTING: 0.0005
  REAL(CMISSRP) :: ElasticityTimeStep = 0.1_CMISSRP ! ### PAPERBRANCH SETTING
  INTEGER(CMISSIntg) :: OdeNSteps = -1 ! can be used to set ODETimeStep implicitly.
  INTEGER(CMISSIntg) :: PdeNSteps = -1  ! overrides PDETimeStep

  REAL(CMISSRP) :: StimDuration=0.1_CMISSRP ! ### PAPERBRANCH SETTING ! should be the same as ElasticityTimeStep

  INTEGER(CMISSIntg)  :: OutputTimestepStride=10  ! (10)

  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------

  !stimulation current in [uA/cm^2]
  REAL(CMISSRP) :: StimValue = 1200.0_CMISSRP ! ### PAPERBRANCH SETTING ! will be applied to all nodes in PhysicalStimulationLength.

  REAL(CMISSRP) :: PMax=7.3_CMISSRP ! ### PAPERBRANCH SETTING ! N/cm^2

  !conductivity in [mS/cm] - This is sigma
  REAL(CMISSRP) :: Conductivity=3.828_CMISSRP ! ### PAPERBRANCH SETTING

  !surface area to volume ratio in [cm^-1]
  REAL(CMISSRP) :: Am=500.0_CMISSRP ! ### PAPERBRANCH SETTING

  !membrane capacitance in [uF/cm^2]
  REAL(CMISSRP) :: CmFast=0.58_CMISSRP ! ### PAPERBRANCH SETTING: Which one is needed? fast/slow? I guess slow.
  REAL(CMISSRP) :: CmSlow=0.58_CMISSRP ! ### PAPERBRANCH SETTING: Which one is needed? fast/slow? I guess slow.

  !Material-Parameters C=[mu_1, mu_2, mu_3, alpha_1, alpha_2, alpha_3, mu_0, XB_stiffness (= e or \eta in different papers)]
  REAL(CMISSRP), PARAMETER, DIMENSION(8) :: C = &
    & [0.0085_CMISSRP*5.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP, & ! [N/cm^2 = 10^4 J/m^3]            ! ### NOT RELEVANT FOR PAPER (other scenario)
    &  11.0_CMISSRP,1.0_CMISSRP,6.0_CMISSRP, &                                                  ! ### NOT RELEVANT FOR PAPER
    &  1.0_CMISSRP,2.2e-9_CMISSRP]                                                              ! ### NOT RELEVANT FOR PAPER

  !maximum contraction velocity in [cm/ms] not relevant for paper
  REAL(CMISSRP) :: VMax=-0.02_CMISSRP ! =0.2 m/s, rat GM                                        ! ### NOT RELEVANT FOR PAPER (isometric conditions)

  !CAUTION - what are the units???   ![N/cm^2]?
!  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE= &
!    &[0.0000000000635201_CMISSRP,0.3626712895523322_CMISSRP,0.0000027562837093_CMISSRP,43.372873938671383_CMISSRP] ! ### SKIPPED due to citation differences

  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE = [0.00000000006352_CMISSRP,0.3627_CMISSRP,0.000002756_CMISSRP,43.373_CMISSRP] ! ### PAPERBRANCH SETTING

  REAL(CMISSRP) :: TkLinParam=1.0_CMISSRP ! 1: With Actin-Tintin Interaction 0: No Actin-Titin Interactions    ! ### NOT RELEVANT FOR PAPER

  !Inital Conditions
  REAL(CMISSRP) :: InitialStretch=1.0_CMISSRP   ! previous value in new mechanical description: 1.2_CMISSRP
  
  INTEGER(CMISSIntg) :: ElasticityLoopMaximumNumberOfIterations = 5     ! only relevant if InitialStretch /= 1.0
  INTEGER(CMISSIntg) :: NewtonMaximumNumberOfIterations = 500           ! 
  REAL(CMISSRP) :: NewtonTolerance = 1.E-8_CMISSRP                      ! 
  REAL(CMISSRP) :: DAERelativeTolerance = 1.E-7_CMISSRP, DAEAbsoluteTolerance = 1.E-7_CMISSRP   ! ### Something we dont have to tell about.

  !--------------------------------------------------------------------------------------------------------------------------------

  ! ### PAPERBRANCH SETTING DESCRIPTION:
  ! we have 36 fibers spread over 2X2X2 FEM elements. Each fiber has 32 elements, 33 nodes (cells). The 16nth cell is stimulated (counting from 1 to 32):
  INTEGER(CMISSIntg) :: NumberGlobalXElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberGlobalYElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberGlobalZElements      ! ### PAPERBRANCH SETTING, set later!
  INTEGER(CMISSIntg) :: NumberOfInSeriesFibres = 1 ! ### has no good effect. (set /=1 to break code) TO BE DISCUSSED? Only Benni will need this. Default 1 makes sense.
  INTEGER(CMISSIntg) :: NumberOfNodesInXi1 = 16    ! ### PAPERBRANCH SETTING, number of 1D elements per 3D element, i.e. number of nodes is +1 (but nodes on edges of 3d elements are shared)
  INTEGER(CMISSIntg) :: NumberOfNodesInXi2 = 3     ! ### PAPERBRANCH SETTING
  INTEGER(CMISSIntg) :: NumberOfNodesInXi3 = 3     ! ### PAPERBRANCH SETTING
    
  INTEGER(CMISSLIntg) :: NumberOfElementsFE
  INTEGER(CMISSIntg) :: NumberOfNodesM
  INTEGER(CMISSIntg) :: NumberOfElementsM
  INTEGER(CMISSIntg) :: NumberOfFibres
  INTEGER(CMISSIntg) :: NumberOfNodesPerLongFibre   ! fibre that touches right boundary has one additional electricity node
  INTEGER(CMISSIntg) :: NumberOfNodesPerShortFibre  ! the number of nodes on ordinary fibres not lying on the rightt boundary
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomX
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomY
  INTEGER(CMISSIntg) :: NumberOfElementsInAtomZ
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainX = 0  ! if zero it will be computed from size of atoms
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainY = 0
  INTEGER(CMISSIntg) :: NumberOfElementsInSubdomainZ = 0
  INTEGER(CMISSIntg) :: NumberOfElementsMInXi1
  INTEGER(CMISSIntg) :: NumberGlobalYFibres
  INTEGER(CMISSINTg) :: NumberOfFibreLinesPerGlobalElement
  INTEGER(CMISSIntg) :: NumberOfGlobalElementLines
  INTEGER(CMISSIntg) :: NumberOfFibreLinesTotal
  INTEGER(CMISSIntg) :: NumberOfElementsMPerFibre ! this is the total number of (0D sub-)cells (nodes) per global fibre.
  INTEGER(CMISSIntg) :: NumberOfElementsMPerFibreLine
  INTEGER(CMISSIntg) :: NumberOfNodesMPerFibreLine
  INTEGER(CMISSIntg) :: NumberOfSubdomainsX = 0, NumberOfSubdomainsY = 0, NumberOfSubdomainsZ = 0    ! if NumberOfElementsInSubdomain is set, NumberOfSubdomains will also be used. Otherwise is will be computed.
  INTEGER(CMISSIntg) :: NumberOfAtomsPerSubdomainX, NumberOfAtomsPerSubdomainY, NumberOfAtomsPerSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfAtomsLastSubdomainX, NumberOfAtomsLastSubdomainY, NumberOfAtomsLastSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfElementsLastAtomX, NumberOfElementsLastAtomY, NumberOfElementsLastAtomZ
  INTEGER(CMISSIntg) :: NumberStimulatedNodesPerFibre
  INTEGER(CMISSIntg) :: PretendedNumberOfDomainsForDomainDecomposition = 0
  
  INTEGER(CMISSIntg) :: Stat
  CHARACTER(len=256) :: CellMLModelFilename = "standard" ! standard will be replaced by the standard model file
  CHARACTER(len=1024) :: inputDirectory = "input/"
  CHARACTER(len=1024) :: FiringTimesFile = "MU_firing_times_10s.txt"
  CHARACTER(len=1024) :: InnervationZoneFile = "innervation_zone_18.txt"
  CHARACTER(len=1024) :: FibreDistributionFile = "MU_fibre_distribution_4050.txt"
  CHARACTER(len=256) :: MemoryConsumption1StTimeStep = "", MemoryConsumptionBeforeSim, Temp
  CHARACTER(len=10000) :: WorkingDirectory
  CHARACTER(len=1024) :: ScenarioName = ""
  
  LOGICAL :: CustomProfilingEnabled !< If custom profiling is compiled in
  LOGICAL :: TauProfilingEnabled !< If TAU profiling is compiled in

  INTEGER(CMISSIntg) :: Ftype,fibre_nr,NearestGP,InElement

  LOGICAL :: less_info,fast_twitch

  REAL(CMISSRP) :: TimeStart, TimeInitFinshed, TimeStretchSimFinished, TimeMainSimulationStart, TimeMainSimulationFinished
  REAL(CMISSSP), DIMENSION(2) :: DurationSystemUser     ! For receiving user and system time
  REAL(CMISSSP) :: DurationTotal
  REAL(CMISSRP) :: StimValuePerNode
  
  INTEGER(CMISSIntg) :: CustomSolverConvergenceReasonParabolic = 0
  INTEGER(CMISSIntg) :: CustomSolverConvergenceReasonNewton = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolic = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolicMin = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolicMax = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewton = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewtonMin = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewtonMax = 0
  INTEGER(CMISSIntg) :: MonodomainSolverId = 2
  INTEGER(CMISSIntg) :: MonodomainPreconditionerId = 1

  INTEGER(CMISSIntg) :: ODESolverId = 1

  INTEGER(CMISSIntg), DIMENSION(10000,100) :: MotorUnitFiringTimes
  INTEGER(CMISSIntg), ALLOCATABLE :: InnervationZoneOffset(:)
  INTEGER(CMISSIntg) :: JunctionNodeNo
  INTEGER(CMISSIntg) :: MotorUnitFires, MotorUnitRank

  INTEGER(CMISSIntg), ALLOCATABLE :: MuDistribution(:)
  REAL(CMISSRP) :: CustomTimingOdeSolver, CustomTimingParabolicSolver, CustomTimingFESolver, CustomTimingFileOutputUser, &
    & CustomTimingFileOutputSystem
  REAL(CMISSRP) :: CustomTimingOdeSolverPreLoad, CustomTimingParabolicSolverPreLoad, CustomTimingFESolverPreLoad, &
    & CustomTimingFileOutputUserPreLoad, CustomTimingFileOutputSystemPreLoad
  REAL(CMISSRP) :: TimingExportEMGUser = 0_CMISSRP
  REAL(CMISSRP) :: TimingExportEMGSystem = 0_CMISSRP

  !--------------------------------------------------------------------------------------------------------------------------------
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensionsFE=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesFE=2

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponentsFE=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MonodomainMeshComponentNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=4
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberFE=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=6
  
  
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=NumberOfSpatialCoordinates

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesM=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsM=3 !Am, Cm, Conductiity   !(scalar, since 1D)

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesM=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsFE=NumberOfSpatialCoordinates+1


  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberM=9
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM1=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM2=6 !5

  INTEGER(CMISSIntg), PARAMETER :: FieldEquationsSetUserNumberM=10

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: SolverDAEIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverParabolicIndex=2
  INTEGER(CMISSIntg), PARAMETER :: SolverFEIndex=1

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopMonodomainNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopElasticityNumber=2

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndexM,EquationsSetIndexFE
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,ComponentNumber,ElementDomain
  INTEGER(CMISSIntg) :: FibreNo, k
  INTEGER(CMISSIntg) :: NumberFiringFibres

  INTEGER(CMISSIntg), ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: shortenModelIndex!,shortenModelIndex2
  INTEGER(CMISSIntg) :: StimComponent

!  REAL(CMISSRP) :: YVALUE
  REAL(CMISSRP) :: VALUE

  INTEGER(CMISSIntg) :: Err, NodeIdx


  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis,LinearBasis,LinearBasisM
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoopMain
  TYPE(cmfe_ControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystemFE,CoordinateSystemM,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: DecompositionFE
  TYPE(cmfe_DecompositionType) :: DecompositionM
  TYPE(cmfe_EquationsType) :: EquationsFE
  TYPE(cmfe_EquationsType) :: EquationsM
  TYPE(cmfe_EquationsSetType) :: EquationsSetFE
  TYPE(cmfe_EquationsSetType) :: EquationsSetM
  TYPE(cmfe_FieldType) :: EquationsSetFieldFE
  TYPE(cmfe_FieldType) :: GeometricFieldFE
  TYPE(cmfe_FieldType) :: DependentFieldFE
  TYPE(cmfe_FieldType) :: IndependentFieldFE
  TYPE(cmfe_FieldType) :: MaterialFieldFE
  TYPE(cmfe_FieldType) :: EquationsSetFieldM
  TYPE(cmfe_FieldType) :: DependentFieldM
  TYPE(cmfe_FieldType) :: GeometricFieldM
  TYPE(cmfe_FieldType) :: IndependentFieldM
  TYPE(cmfe_FieldType) :: MaterialFieldM
  TYPE(cmfe_FieldType) :: FibreField
  TYPE(cmfe_FieldType) :: CellMLModelsField
  TYPE(cmfe_FieldType) :: CellMLStateField
  TYPE(cmfe_FieldType) :: CellMLIntermediateField
  TYPE(cmfe_FieldType) :: CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMeshFE
  TYPE(cmfe_MeshType) :: MeshFE
  TYPE(cmfe_MeshType) :: MeshM
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: RegionFE,RegionM,WorldRegion
  TYPE(cmfe_SolverType) :: SolverDAE,SolverParabolic
  TYPE(cmfe_SolverType) :: SolverFE,LinearSolverFE
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsM,SolverEquationsFE
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements
  TYPE(cmfe_MeshElementsType) :: ElementsM


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

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

!##################################################################################################################################
!##################################################################################################################################
!##################################################################################################################################

  WRITE(*,'(A,A)') TRIM(GetTimeStamp()), " Program started."

  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeStart = DurationSystemUser(2)
    
  CALL ParseParameters()
    
  !================================================================================================================================
  !  G E N E R A L   F E A T U R E S
  !================================================================================================================================
  CALL CreateRegionMesh()
  CALL CreateDecompositionFiniteElasticity()
  CALL CreateDecompositionMonodomain()

  !================================================================================================================================
  !  F I N I T E   E L A S T I C I T Y
  !================================================================================================================================
  CALL CreateFieldFiniteElasticity()

  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================

  CALL CreateFieldMonodomain()

  !================================================================================================================================
  !  EQUATIONS SET
  CALL CreateEquationsSet()
  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  CALL InitializeFieldMonodomain()

  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  CALL InitializeFieldFiniteElasticity()

  CALL CreateEquations()
  CALL InitializeCellML()
  CALL CreateProblem()

  CALL CreateControlLoops()
  CALL CreateSolvers()
  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions
  CALL SetBoundaryConditions()
  
  ! Output the data structure Problem
  !IF (DEBUGGING_PROBLEM_OUTPUT .AND. ComputationalNodeNumber == 0) THEN
  !  PRINT*, ""
  !  PRINT*, ""
  !  CALL cmfe_PrintProblem(Problem,6,30,Err)
  !  PRINT*, ""
  !  PRINT*, ""
    !PRINT*, "End the program after output of problem datastructure"
    !STOP
  !ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Calculate the bioelectrics geometric field
  CALL CalculateBioelectrics()
  
 !PRINT*, "Size IndependentFieldFE:",cmfe_getFieldSize(IndependentFieldFE, Err),"Bytes"
  
  !PRINT*, "Abort program in FortranExample.f90:450"
  CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
  !STOP
  CALL ExportEMG()
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem

  !Solve the problem -- bring to new length before applying the stimulus
  IF (InitialStretch == 1.0_CMISSRP) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, "1.) Pre-stretch disabled, initial stretch =", InitialStretch
  ELSE
    IF (ComputationalNodeNumber == 0) PRINT *, "1.) Pre-stretch simulation, initial stretch =", InitialStretch
  ENDIF
  Temp = GetMemoryConsumption()
  CALL cmfe_CustomSolverInfoReset(Err)
  IF (DEBUGGING_PARALLEL_BARRIER) CALL gdbParallelDebuggingBarrier()

  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeInitFinshed = DurationSystemUser(2)
  
  IF (InitialStretch /= 1.0_CMISSRP) THEN
    CALL cmfe_Problem_Solve(Problem,Err) ! ### PAPER SETTING: not necessary
  ENDIF
    
  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeStretchSimFinished = DurationSystemUser(2)
  
  CALL HandleSolverInfo(-1.0_CMISSRP)

  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ParameterSetUpdateConstant(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,PMax, Err)
  ENDIF


  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)

! no change for BCs -- fix at this length!!!

  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  IF (ComputationalNodeNumber == 0) PRINT*, "2.) Simulate with stimulation"
  MemoryConsumptionBeforeSim = GetMemoryConsumption()

  CALL cmfe_CustomTimingGet(CustomTimingOdeSolverPreLoad, CustomTimingParabolicSolverPreLoad, &
    & CustomTimingFESolverPreLoad, CustomTimingFileOutputUserPreLoad, CustomTimingFileOutputSystemPreLoad, Err)
  CALL cmfe_CustomTimingReset(Err)
  CALL cmfe_CustomProfilingReset(Err)   ! reset and discard all custom profiling data that was recorded so far
  
  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeMainSimulationStart = DurationSystemUser(2)

  time = 0.0_CMISSRP
  VALUE = 0.0_CMISSRP
  k = 1       ! row in firing_times input (time)
  ! main time loop
  DO WHILE(time < TimeStop-1e-10)
  
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)

    IF (ComputationalNodeNumber == 0) PRINT "(A,F0.5,A)","t = ",time," s"
    !-------------------------------------------------------------------------------------------------------------------------------
    !Set the Stimulus for monodomain at the middle of the fibres
    IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics

      !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
      !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
      !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
      !! - is the model ID also needed?
      !! - because the CellML fields should all be set up to allow direct use in the CellML code, the component number matches the index of the given variable in its associated array in the CellML generated code.
      !SUBROUTINE CELLML_FIELD_COMPONENT_GET_C(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,COMPONENT_USER_NUMBER,ERR,ERROR,*)
       !Argument variables
       !TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
       !INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map from.
       !INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
       !CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the model variable which needs to be located in the provided field.
       !INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given ID.
      CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
        & "wal_environment/I_HH",StimComponent,Err)
    
    ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
      CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
        & "Aliev_Panfilov/I_HH",StimComponent,Err)
        
    ELSEIF (ModelType == 2) THEN ! 4, "Titin"
      CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
        & "Aliev_Panfilov/I_HH",StimComponent,Err)
    ENDIF

  !  VALUE = VALUE-ABS(Vmax)/20.0_CMISSRP*StimDuration
  !!  VALUE = VALUE+ABS(Vmax)/10.0_CMISSRP*StimDuration
  !  CALL cmfe_ControlLoop_BoundaryConditionUpdate(ControlLoopFE,1,1,VALUE,Err)
  !  DO NodeIdx=1,SIZE(LeftSurfaceNodes,1)
  !    NodeNumber=LeftSurfaceNodes(NodeIdx)
  !    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
  !    IF(NodeDomain==ComputationalNodeNumber) THEN
  !      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
  !        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
  !    ENDIF
  !  ENDDO

    CALL SetStimulationAtNodes(StimValuePerNode)

    IF (ComputationalNodeNumber == 0) THEN
      PRINT "(A,I4,A,I6)", "   Number of stimulated fibres on process 0: ", NumberFiringFibres," of",NumberOfFibres
    ENDIF

    !-------------------------------------------------------------------------------------------------------------------------------
    !Solve the problem for the stimulation time
    IF (ComputationalNodeNumber == 0) Print*, "  Solve with stimulation,    time span: ", time, " to ",time+StimDuration
    CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time,time+StimDuration,ElasticityTimeStep,Err)

    CALL cmfe_CustomSolverInfoReset(Err)
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)
    
    CALL cmfe_Problem_Solve(Problem,Err)
    
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)
    CALL HandleSolverInfo(time)
    Temp = GetMemoryConsumption()
    IF (DebuggingOnlyRunShortPartOfSimulation) EXIT

    !-------------------------------------------------------------------------------------------------------------------------------
    !Now turn the stimulus off
    CALL SetStimulationAtNodes(0.0_CMISSRP)
    
    !-------------------------------------------------------------------------------------------------------------------------------
    !Solve the problem for the rest of the period
    IF (ComputationalNodeNumber == 0) PRINT*, "  Solve without stimulation, time span: ", time+StimDuration, " to ",time+StimPeriod
    CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time+StimDuration,time+StimPeriod,ElasticityTimeStep,Err)

    CALL cmfe_CustomSolverInfoReset(Err)    
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)
    
    CALL cmfe_Problem_Solve(Problem,Err)
    
    CALL cmfe_CustomProfilingStart("level 0: stimulation handling",Err)
    CALL HandleSolverInfo(time+StimDuration)
    !-------------------------------------------------------------------------------------------------------------------------------
    time = time + StimPeriod
    k = k+1

    IF (k == 2) THEN
      MemoryConsumption1StTimeStep = GetMemoryConsumption()
    ELSE
      Temp = GetMemoryConsumption()
    ENDIF
    CALL cmfe_CustomProfilingStop("level 0: stimulation handling",Err)

  ENDDO

  ! store duration
  CALL ETIME(DurationSystemUser, DurationTotal)
  TimeMainSimulationFinished = DurationSystemUser(2)
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------------------------------------------------
  CALL ExportEMG()
  CALL cmfe_TimingSummaryOutput(Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  CALL cmfe_CustomTimingGet(CustomTimingOdeSolver, CustomTimingParabolicSolver, CustomTimingFESolver, CustomTimingFileOutputUser, &
    & CustomTimingFileOutputSystem, Err)

  CALL WriteTimingFile()

  CALL WriteCustomProfilingFile()

  IF (ComputationalNodeNumber == 0 .AND. CustomProfilingEnabled) THEN
    PRINT*, TRIM(cmfe_CustomProfilingGetInfo(Err))
  ENDIF

  PRINT*, ""
  PRINT*, "--------------------------------------------------"
  PRINT*, "Process ", ComputationalNodeNumber
  PRINT*, "Timing (user time):"
  PRINT*, "   Ode Solver:       preload: ", CustomTimingOdeSolverPreLoad, " s, main: ", CustomTimingOdeSolver, " s"
  PRINT*, "   Parabolic Solver: preload: ", CustomTimingParabolicSolverPreLoad, " s, main: ", CustomTimingParabolicSolver, " s"
  PRINT*, "   3D FE Solver:     preload: ", CustomTimingFESolverPreLoad, " s, main: ", CustomTimingFESolver, " s"
  PRINT*, "   Node File Output: preload: ", CustomTimingFileOutputUserPreLoad, " s, main: ", CustomTimingFileOutputUser, " s"
  PRINT*, "           (system): preload: ", CustomTimingFileOutputSystemPreLoad, " s, main: ", CustomTimingFileOutputSystem, " s"
  PRINT*, "   Total Simulation: preload: ", (TimeStretchSimFinished - TimeInitFinshed), " s, main: ", &
    & (TimeMainSimulationFinished-TimeMainSimulationStart), " s"
  PRINT*, "   EMG Output: user: ", TimingExportEMGUser, " s, system: ", TimingExportEMGSystem, " s"
  PRINT*, ""

  WRITE(*,'(A,A)') TRIM(GetTimeStamp()), " Program successfully completed."
  STOP
  
CONTAINS

! Test whether parametrization and number of processes matches and is valid
FUNCTION CheckGeometry()
  LOGICAL :: CheckGeometry
  CheckGeometry = .TRUE.

  IF (NumberGlobalXElements <= 0 .OR. NumberGlobalYElements <= 0 .OR. NumberGlobalZElements <= 0 &
    & .OR. NumberOfInSeriesFibres <= 0) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: Number of elements (", NumberGlobalXElements, ", ",NumberGlobalYElements, ", ", &
        & NumberGlobalZElements, ", ", NumberOfInSeriesFibres, ") is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomX <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomX=", NumberOfElementsInAtomX, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomY <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomY=", NumberOfElementsInAtomY, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomZ <= 0 ) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsInAtomZ=", NumberOfElementsInAtomZ, " is invalid!"
    ENDIF
    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsMPerFibre * NumberOfInSeriesFibres /= NumberOfElementsMPerFibreLine) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Error: NumberOfElementsMPerFibreLine=",NumberOfElementsMPerFibreLine, &
        & " is not an integer multiple of F=",NumberOfInSeriesFibres,"!"
    ENDIF

    CheckGeometry = .FALSE.
  ENDIF

  IF (NumberOfElementsInAtomX > NumberGlobalXElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomX=",NumberOfElementsInAtomX, &
        & " to value of NumberGlobalXElements=",NumberGlobalXElements,"."
    ENDIF
    NumberOfElementsInAtomX = NumberGlobalXElements
  ENDIF
  
  IF (NumberOfElementsInAtomY > NumberGlobalYElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomY=",NumberOfElementsInAtomY, &
        & " to value of NumberGlobalYElements=",NumberGlobalYElements,"."
    ENDIF
    NumberOfElementsInAtomY = NumberGlobalYElements
  ENDIF
  
  IF (NumberOfElementsInAtomZ > NumberGlobalZElements) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT*, "Warning: Decreasing NumberOfElementsInAtomZ=",NumberOfElementsInAtomZ, &
        & " to value of NumberGlobalZElements=",NumberGlobalZElements,"."
    ENDIF
    NumberOfElementsInAtomZ = NumberGlobalZElements
  ENDIF

END FUNCTION CheckGeometry

SUBROUTINE ReadInputFiles()
  LOGICAL :: FileExists
  INTEGER :: Status
  INTEGER :: I
  CHARACTER(len=100000) :: Buffer
  INTEGER :: NumberOfEntries

  ! Read in firing times file
  FiringTimesFile = TRIM(InputDirectory) // TRIM(FiringTimesFile)
  INQUIRE(file=FiringTimesFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(FiringTimesFile) // """ does not exist!"
    STOP
  ENDIF

  IF (ComputationalNodeNumber == 0) PRINT*,  "Open file """ // TRIM(FiringTimesFile) // """."

  OPEN(unit=5, file=FiringTimesFile, action="read", iostat=Status)

  ! loop over maximum 10000 lines
  DO I = 1, 10000
    READ(5,*,iostat=Status) MotorUnitFiringTimes(I,:)
    IF (Status /= 0) THEN
      IF (ComputationalNodeNumber == 0) PRINT*,  "File """ // TRIM(FiringTimesFile) // """ contains ", I, " firing times."
      EXIT
    ENDIF
  ENDDO
  CLOSE(unit=5)


  !--------------------------------------------------------------------------------------------------------------------------------
  ! Read in InnervationZoneFile
  ! Gaussian distribution with mean 0 and std 2 (MATLAB: k = 2*randn(len,1), kk = int64(k);)
  InnervationZoneFile = TRIM(InputDirectory) // TRIM(InnervationZoneFile)
  INQUIRE(file=InnervationZoneFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(InnervationZoneFile) // """ does not exist!"
    STOP
  ENDIF

  OPEN(unit=5, file=InnervationZoneFile, action="read", iostat=Status)
  READ(unit=5, fmt='(A)') Buffer

  ! Determine the number of entries in the first line
  NumberOfEntries = 1
  DO I = 0, LEN_TRIM(Buffer)
    IF (Buffer(I:I) == ',') THEN
      NumberOfEntries = NumberOfEntries + 1
    ENDIF
  ENDDO

  REWIND(UNIT=5)
  IF (ComputationalNodeNumber == 0) THEN
    PRINT*, "File  """ // TRIM(InnervationZoneFile) // """ contains ", NumberOfEntries, " Entries, NumberOfFibres=", &
      & NumberOfFibres, "."
  ENDIF

  ALLOCATE(InnervationZoneOffset(NumberOfFibres), Stat=stat)
  
  IF (stat /= 0) THEN
    PRINT *, "Could not allocate ", NumberOfFibres, " objects for InnervationZoneOffset!"
    STOP
  ENDIF

  IF (NumberOfFibres <= NumberOfEntries) THEN
    READ(unit=5, fmt=*, iostat=Status) InnervationZoneOffset(1:NumberOfFibres)
  ELSE
    READ(unit=5, fmt=*, iostat=Status) InnervationZoneOffset(1:NumberOfEntries)

    ! fill with already read data
    DO I = NumberOfEntries,NumberOfFibres
      InnervationZoneOffset(I) = InnervationZoneOffset(MOD(I, NumberOfEntries)+1)
      !PRINT*, "fill I=",I," with entry ", InnervationZoneOffset(MOD(I, NumberOfEntries)+1), " from index ", &
      !  & MOD(I, NumberOfEntries)+1
    ENDDO
  ENDIF
  CLOSE(unit=5)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Read in motor unit fibre distribution
  !the MU (=motor unit) number the fiber belongs to
  FibreDistributionFile = TRIM(InputDirectory) // TRIM(FibreDistributionFile)
  INQUIRE(file=FibreDistributionFile, exist=FileExists)

  IF (.NOT. FileExists) THEN
    PRINT*, "Error: File """ // TRIM(FibreDistributionFile) // """ does not exist!"
    STOP
  ENDIF

  OPEN(unit=5, file=FibreDistributionFile, action="read", iostat=Status)
  READ(unit=5, fmt='(A)') Buffer

  ! Determine the number of entries in the first line
  NumberOfEntries = 1
  DO I = 0, LEN_TRIM(Buffer)
    IF (Buffer(I:I) == ',') THEN
      NumberOfEntries = NumberOfEntries + 1
    ENDIF
  ENDDO

  REWIND(UNIT=5)
  IF (ComputationalNodeNumber == 0) THEN
    PRINT*, "File  """ // TRIM(FibreDistributionFile) // """ contains ", NumberOfEntries, " Entries, NumberOfFibres=", &
      & NumberOfFibres, "."
  ENDIF

  ALLOCATE(MUDistribution(NumberOfFibres))

  IF (NumberOfFibres <= NumberOfEntries) THEN
    READ(unit=5, fmt=*, iostat=Status) MUDistribution(1:NumberOfFibres)
  ELSE
    READ(unit=5, fmt=*, iostat=Status) MUDistribution(1:NumberOfEntries)

    ! fill with already read data
    DO I = NumberOfEntries,NumberOfFibres
      MUDistribution(I) = MUDistribution(MOD(I, NumberOfEntries)+1)
      !PRINT*, "fill I=",I," with entry ", MUDistribution(MOD(I, NumberOfEntries)+1), " from index ", &
      !  & MOD(I, NumberOfEntries)+1
    ENDDO
  ENDIF
  CLOSE(unit=5)
END SUBROUTINE ReadInputFiles

SUBROUTINE ToLower(Str)
  CHARACTER(*), INTENT(INOUT) :: Str
  INTEGER(CMISSIntg) :: I
  DO I = 1, LEN(Str)
    SELECT CASE(Str(I:I))
      CASE("A":"Z")
        STR(I:I) = ACHAR(IACHAR(Str(I:I))+32)
    END SELECT
  ENDDO
END SUBROUTINE ToLower

SUBROUTINE ParseAssignment(Line, LineNumber, ScenarioInputFile)
  CHARACTER(LEN=*), INTENT(IN) :: Line
  INTEGER(CMISSIntg), INTENT(IN) :: LineNumber
  CHARACTER(LEN=*), INTENT(IN) :: ScenarioInputFile
  CHARACTER(LEN=1024) :: VariableName, StrValue
  INTEGER(CMISSIntg) :: StringLength
  LOGICAL :: OldTomoMechanics   

  !PRINT *, "Read line """//TRIM(Line)//"""."
  
  ! Extract variable name and value
  IF (INDEX(Line, "=") /= 0) THEN
    VariableName = TRIM(ADJUSTL(Line(1:INDEX(Line, "=")-1)))
    CALL ToLower(VariableName)
    StrValue = TRIM(ADJUSTL(Line(INDEX(Line, "=")+1:)))
    
    !PRINT *, "VariableName=["//TRIM(VariableName)//"], Value=["//TRIM(StrValue)//"]."
  ELSE
    RETURN
  ENDIF
  
  ! Parse variables
  SELECT CASE(VariableName)
    CASE ("x","numberglobalxelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalXElements
    CASE ("y","numberglobalyelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalYElements
    CASE ("z","numberglobalzelements")
      READ(StrValue, *, IOSTAT=Stat) NumberGlobalZElements
    CASE ("f","numberofinseriesfibres")
      READ(StrValue, *, IOSTAT=Stat) NumberOfInSeriesFibres
    CASE ("a","numberofelementsinatomicportionperdomain","numberofelementsinatomx","ax")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomX
    CASE ("numberofelementsinatomy","ay")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomY
    CASE ("numberofelementsinatomz","az")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInAtomZ
    CASE ("numberofelementsinsubdomainx","ex")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainX
    CASE ("numberofelementsinsubdomainy","ey")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainY
    CASE ("numberofelementsinsubdomainz","ez")
      READ(StrValue, *, IOSTAT=Stat) NumberOfElementsInSubdomainZ
    CASE ("numberofsubdomainsx","px")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsX
    CASE ("numberofsubdomainsy","py")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsY
    CASE ("numberofsubdomainsz","pz")
      READ(StrValue, *, IOSTAT=Stat) NumberOfSubdomainsZ
    CASE ("odesolverid")
      READ(StrValue, *, IOSTAT=Stat) ODESolverId
    CASE ("monodomainsolverid")
      READ(StrValue, *, IOSTAT=Stat) MonodomainSolverId
    CASE ("monodomainpreconditionerid")
      READ(StrValue, *, IOSTAT=Stat) MonodomainPreconditionerId
    CASE ("t","tend","timestop")
      READ(StrValue, *, IOSTAT=Stat) TimeStop
    CASE ("odetimestep")
      READ(StrValue, *, IOSTAT=Stat) ODETimeStep
    CASE ("pdetimestep")
      READ(StrValue, *, IOSTAT=Stat) PDETimeStep
    CASE ("elasticitytimestep")
      READ(StrValue, *, IOSTAT=Stat) ElasticityTimeStep
    CASE ("stimvalue")
      READ(StrValue, *, IOSTAT=Stat) StimValue
    CASE ("stimduration", "stimstop")
      READ(StrValue, *, IOSTAT=Stat) StimDuration
    CASE ("stimperiod", "periodd")
      READ(StrValue, *, IOSTAT=Stat) StimPeriod
    CASE ("outputtimestepstride")
      READ(StrValue, *, IOSTAT=Stat) OutputTimeStepStride
    CASE ("xi1", "numberofnodesinxi1")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi1
    CASE ("xi2", "numberofnodesinxi2")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi2
    CASE ("xi3", "numberofnodesinxi3")
      READ(StrValue, *, IOSTAT=Stat) NumberOfNodesInXi3
    CASE ("pretendednumberofdomainsfordomaindecomposition", "pretendednumberofdomains")
      READ(StrValue, *, IOSTAT=Stat) PretendedNumberOfDomainsForDomainDecomposition
    CASE ("newtonmaximumnumberofiterations")
      READ(StrValue, *, IOSTAT=Stat) NewtonMaximumNumberOfIterations
    CASE ("daerelativetolerance")
      READ(StrValue, *, IOSTAT=Stat) DAERelativeTolerance
    CASE ("daeabsolutetolerance")
      READ(StrValue, *, IOSTAT=Stat) DAEAbsoluteTolerance
    CASE ("newtontolerance")
      READ(StrValue, *, IOSTAT=Stat) NewtonTolerance
    CASE ("elasticityloopmaximumnumberofiterations")
      READ(StrValue, *, IOSTAT=Stat) ElasticityLoopMaximumNumberOfIterations
    CASE ("enableexportemg")
      READ(StrValue, *, IOSTAT=Stat) EnableExportEMG
    CASE ("oldtomomechanics")                   ! deprecated variable, use ModelType instead
      READ(StrValue, *, IOSTAT=Stat) OldTomoMechanics
      IF (OldTomoMechanics) THEN
        ModelType = 0
      ELSE
        ! If OldTomoMechanics is false, it means that ModelType=0 should not be used. 
        ! But is then ModelType=1 or ModelType=2 intended? Not clear (assume ModelType=1), thus deprecated.
        
        ModelType = 1
        PRINT *, "Warning: Parameter ""OldTomoMechanics=F"" was given, which is deprecated. Use ModelType=1 or ModelType=2 instead!"
      ENDIF
    CASE ("modeltype")
      READ(StrValue, *, IOSTAT=Stat) ModelType
    CASE ("splittingtype")
      READ(StrValue, *, IOSTAT=Stat) SplittingType
    CASE ("tklinparam")
      READ(StrValue, *, IOSTAT=Stat) TkLinParam
    CASE ("debuggingoutput")
      READ(StrValue, *, IOSTAT=Stat) DebuggingOutput
    CASE ("debuggingonlyrunshortpartofsimulation")
      READ(StrValue, *, IOSTAT=Stat) DebuggingOnlyRunShortPartOfSimulation
    CASE ("elasticitydisabled")
      READ(StrValue, *, IOSTAT=Stat) ElasticityDisabled
    CASE ("physicallength")
      READ(StrValue, *, IOSTAT=Stat) PhysicalLength
    CASE ("physicalwidth")
      READ(StrValue, *, IOSTAT=Stat) PhysicalWidth
    CASE ("physicalheight")
      READ(StrValue, *, IOSTAT=Stat) PhysicalHeight
    CASE ("physicalstimulationlength", "physicalstimulationwidth")
      READ(StrValue, *, IOSTAT=Stat) PhysicalStimulationLength
    CASE ("pmax")
      READ(StrValue, *, IOSTAT=Stat) PMax
    CASE ("conductivity")
      READ(StrValue, *, IOSTAT=Stat) Conductivity
    CASE ("am")
      READ(StrValue, *, IOSTAT=Stat) Am
    CASE ("cmfast")
      READ(StrValue, *, IOSTAT=Stat) CmFast
    CASE ("cmslow")
      READ(StrValue, *, IOSTAT=Stat) CmSlow
    CASE ("vmax")
      READ(StrValue, *, IOSTAT=Stat) Vmax
    CASE ("initialstretch")
      READ(StrValue, *, IOSTAT=Stat) InitialStretch
    CASE ("odensteps","nodesteps","nstepsode","nsteps")
      READ(StrValue, *, IOSTAT=Stat) OdeNSteps
    CASE ("pdensteps","npdesteps","nstepspde")
      READ(StrValue, *, IOSTAT=Stat) PdeNSteps
    CASE ("inputdirectory")
      InputDirectory = TRIM(ADJUSTL(StrValue))
      
      ! Append slash to input directory if necessary
      StringLength = LEN_TRIM(InputDirectory)

      IF (.NOT. InputDirectory(StringLength:StringLength) == "/") THEN
        InputDirectory(StringLength+1:StringLength+1) = "/"
      ENDIF
      InputDirectory = TRIM(InputDirectory)
    CASE ("firingtimesfile")
      FiringTimesFile = TRIM(ADJUSTL(StrValue))
    CASE ("innervationzonefile")
      InnervationZoneFile = TRIM(ADJUSTL(StrValue))
    CASE ("fibredistributionfile")
      FibreDistributionFile = TRIM(ADJUSTL(StrValue))
    CASE ("cellmlmodelfilename")
      CellMLModelFilename = TRIM(ADJUSTL(StrValue))
    CASE ("scenarioname")
      ScenarioName = TRIM(ADJUSTL(StrValue))
    CASE DEFAULT
      Read(LineNumber, *, IOSTAT=Stat) StrValue
      WRITE(*,'(A)') TRIM(ScenarioInputFile) // ":" // TRIM(StrValue) // ": Unrecognized variable """ // TRIM(VariableName) // """."
  END SELECT
  IF (Stat /= 0) THEN
    Read(LineNumber, *, IOSTAT=Stat) StrValue
    WRITE(*,'(A)') TRIM(ScenarioInputFile) // ": " // &
     &  "Could not parse value """ // TRIM(StrValue) // """ for variable """//TRIM(VariableName)//"""."            
  ENDIF
  
END SUBROUTINE ParseAssignment

SUBROUTINE ParseScenarioInputFile(ScenarioInputFile)
  CHARACTER(LEN=1024), INTENT(IN) :: ScenarioInputFile
  CHARACTER(LEN=10000) :: Line
  INTEGER(CMISSIntg) :: LineNumber, StringLength
  LineNumber = 0  
  
  OPEN(UNIT=100, FILE=ScenarioInputFile, ACTION="read", IOSTAT=Stat)
  
  IF (Stat /= 0) THEN
    PRINT*, "Could not read from Input file """ // TRIM(ScenarioInputFile) // """."
    RETURN
  ENDIF
  
  DO
    READ(UNIT=100, FMT='(A)', IOSTAT=Stat) Line
    LineNumber = LineNumber + 1
    
    IF (Stat < 0) EXIT  ! end of file reached
    
    ! Remove comments starting with # or !
    IF (INDEX(Line, "#") /= 0) THEN
      Line = TRIM(Line(1:INDEX(Line, "#")-1))
    ENDIF
    IF (INDEX(Line, "!") /= 0) THEN
      Line = Line(1:INDEX(Line, "!")-1)
    ENDIF
    
    CALL ParseAssignment(Line, LineNumber, ScenarioInputFile)
    
  ENDDO
  
END SUBROUTINE ParseScenarioInputFile

! Parse command line parameters and set numbers of elements
SUBROUTINE ParseParameters()

  INTEGER(CMISSLINTg) :: Factor, NumberArguments, ValueArgumentCount
  INTEGER(CMISSINTg) :: StringLength
  CHARACTER(LEN=10000) :: Arg
  CHARACTER(LEN=1024) :: ScenarioInputFile
  LOGICAL :: GeometryIsValid, FileExists
  INTEGER(CMISSIntg) :: I

  ! Default values
  NumberGlobalXElements = 2 ! ### PAPER SETTING
  NumberGlobalYElements = 2 ! ### PAPER SETTING
  NumberGlobalZElements = 2 ! ### PAPER SETTING
  NumberOfInSeriesFibres = 1 !1
  NumberOfElementsInAtomX = 1
  NumberOfElementsInAtomY = 1
  NumberOfElementsInAtomZ = 1

  ! loop over arguments
  NumberArguments = IARGC()
  ValueArgumentCount = 1
  DO I = 1, NumberArguments
    CALL GETARG(I, Arg)   ! get argument
    
    ! if argument ends in .sce, consider this as a scenario input file
    IF (INDEX(Arg, ".sce", .TRUE.) /= 0 .AND. INDEX(Arg, ".sce", .TRUE.) == LEN_TRIM(Arg)-3) THEN
      ScenarioInputFile = TRIM(Arg)
      CALL ParseScenarioInputFile(ScenarioInputFile)
      
    ! if argument contains '=' consider it as assignment
    ElSEIF (INDEX(Arg, "=") /= 0) THEN
      CALL ParseAssignment(Arg, I, "command line argument")
    ELSE
      SELECT CASE (ValueArgumentCount)     ! Input directory
      CASE(1)
        InputDirectory = Arg
      
        ! Append slash to input directory if necessary
        StringLength = LEN_TRIM(InputDirectory)

        IF (.NOT. InputDirectory(StringLength:StringLength) == "/") THEN
          InputDirectory(StringLength+1:StringLength+1) = "/"
        ENDIF
        InputDirectory = TRIM(InputDirectory)
      
      CASE(2)
        READ(Arg,*,Iostat=Stat)  NumberGlobalXElements
      CASE(3)
        READ(Arg,*,Iostat=Stat)  NumberGlobalYElements
      CASE(4)
        READ(Arg,*,Iostat=Stat)  NumberGlobalZElements
      CASE(5)
        READ(Arg,*,Iostat=Stat)  NumberOfInSeriesFibres
      CASE(6)
        READ(Arg,*,Iostat=Stat)  NumberOfElementsInAtomX
      CASE(7)
        READ(Arg,*,Iostat=Stat)  ODESolverId
      CASE(8)
        READ(Arg,*,Iostat=Stat)  MonodomainSolverId
      CASE(9)
        READ(Arg,*,Iostat=Stat)  MonodomainPreconditionerId
      CASE(10)
        READ(Arg,*,Iostat=Stat)  OdeNSteps
        
      ENDSELECT
      
      ValueArgumentCount = ValueArgumentCount + 1
    
    ENDIF
  ENDDO
  
  IF (NumberArguments == 0) THEN
    PRINT*, "Using default values. " // NEW_LINE('A') &
     & // "Usage: " // NEW_LINE('A') // &
     & "1)   ./cuboid [<input folder> [<X> <Y> <Z> [<F> [<NumberOfElementsInAtomX> " // &
     & "[<ODEsolver> [<Msolver> [<Mprecond> ]]]]]]] " // NEW_LINE('A') // &
     & "2)   ./cuboid [<variable>=<value>] [<input.sce>] [<variable>=<value>] " // NEW_LINE('A') // &
     & "     See the example scenario file for file format and variable names. Variables will be set in order of the arguments."
  ENDIF
!#################### the following is only necessary to correct user output, i think. ############################################
  IF (PdeNSteps/=-1) THEN
    PDETimeStep = ElasticityTimeStep/PdeNSteps
  ENDIF
  IF (OdeNSteps/=-1) THEN
    ODETimeStep = PDETimeStep/OdeNSteps
  END IF
!##################################################################################################################################

  ! direction of fibres is in Xi1=Global X direction
  NumberOfElementsFE = NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements

  NumberGlobalYFibres = NumberOfNodesInXi2 * NumberGlobalYElements
  NumberOfFibreLinesPerGlobalElement = NumberOfNodesInXi2 * NumberOfNodesInXi3
  NumberOfGlobalElementLines = NumberGlobalYElements * NumberGlobalZElements
  NumberOfFibreLinesTotal = NumberOfFibreLinesPerGlobalElement * NumberOfGlobalElementLines
  NumberOfFibres = NumberOfFibreLinesTotal * NumberOfInSeriesFibres
  NumberOfElementsMInXi1 = NumberOfNodesInXi1 - 1
  NumberOfElementsMPerFibreLine = NumberOfElementsMInXi1 * NumberGlobalXElements
  NumberOfElementsMPerFibre = NumberOfElementsMPerFibreLine / NumberOfInSeriesFibres
  NumberOfNodesPerShortFibre = NumberOfElementsMPerFibre
  NumberOfNodesPerLongFibre = NumberOfElementsMPerFibre + 1
  ! the end point nodes of fibres can not be shared between different fibres, therefore some fibres have 1 more node (the ones at the left)

  ! total number of 1D nodes
  NumberOfNodesMPerFibreLine = NumberOfNodesPerShortFibre * (NumberOfInSeriesFibres-1) + NumberOfNodesPerLongFibre
  NumberOfNodesM = NumberOfNodesMPerFibreLine * NumberOfFibreLinesTotal
  NumberOfElementsM = NumberOfElementsMPerFibre * NumberOfFibres

  ! compute number of bioelectric nodes that will be stimulated
  NumberStimulatedNodesPerFibre = MAX(1, NINT(DBLE(PhysicalStimulationLength) * (NumberOfElementsMPerFibre / PhysicalWidth)))
  
  ! previous implementation changed:: StimValuePerNode = StimValue / NumberStimulatedNodesPerFibre to::
  StimValuePerNode = StimValue
  ! ^ This was done because we need to keep the stimulation in a cell constant, no matter, how many cells there are. the dynamic of a cell is independent of the spacial discretization. Checked by Nehzat, Thomas, Aaron
!##################################################################################################################################
!  fast_twitch=.true.
!  if(fast_twitch) then
!  pathname="/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/"
!  CellMLModelFilename=trim(pathname)//"fast_2014_03_25_no_Fl_no_Fv.xml" !FAST
 
  IF (TRIM(CellMLModelFilename) == "standard") THEN
    IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
      CellMLModelFilename = TRIM(inputDirectory) // "new_slow_TK_2014_12_08.cellml"
      !StimValue = 20000.0_CMISSRP !700.0_CMISSRP!700.0_CMISSRP

    ELSEIF (ModelType == 1) THEN ! 3, , "MultiPhysStrain", numerically more stable
      CellMLModelFilename = TRIM(inputDirectory) // "Aliev_Panfilov_Razumova_2016_08_22.cellml"
      !StimValue = 90.0_CMISSRP !90.0_CMISSRP

    ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
      CellMLModelFilename = TRIM(inputDirectory) // "Aliev_Panfilov_Razumova_Titin_2016_10_10_noFv.cellml"
      
    ENDIF
  ENDIF

  ! check if file exists
  INQUIRE(file=CellMLModelFilename, exist=FileExists)
  IF (.NOT. FileExists) THEN
    PRINT*, "Error: CellML file """ // TRIM(CellMLModelFilename) // """ does not exist!"
    STOP
  ENDIF


!   &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/shorten_mod_2011_07_04.xml"
!    pathname="/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles"
!    CellMLModelFilename=trim(pathname)//"/fast_shortening_0.1vmax.xml"
!    StimValue=700.0_CMISSRP!2000.0_CMISSRP!700.0_CMISSRP
!  else !slow twitch
!    filename2= &
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_stim_2012_07_23.xml"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_2012_07_23.xml_0.401"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_twitch_2012_01_27.xml"
!    StimValue=2000.0_CMISSRP
!  endif
!##################################################################################################################################

 !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  ! set diagnostics
!  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"dignostics",["FIELD_MAPPINGS_CALCULATE"],err)

  !                          type                level list,  output file  routine list
  !CALL cmfe_DiagnosticsSetOn(cmfe_FROM_DIAG_TYPE, [3], "", ["cmfe_Problem_Solve", "PROBLEM_SOLVE     "], Err)

  !CALL cmfe_DiagnosticsSetOn(cmfe_CALL_STACK_DIAG_TYPE, [5], "diagnostics", ["FIELD_MAPPINGS_CALCULATE"], Err)

  !CALL cmfe_OutputSetOn("output.txt", Err)
  
  !CMFE_ALL_DIAG_TYPE    !<Type for setting diagnostic output in all routines \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !CMFE_IN_DIAG_TYPE     !<Type for setting diagnostic output in one routine \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !CMFE_FROM_DIAG_TYPE   !<Type for setting diagnostic output in one routine downwards \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !                          in which routine,   levelList(:), diagFilename
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE, [1],          "",&
  ! routine list
  !  & ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"], Err)
  !  & ["FIELD_MAPPINGS_CALCULATE"], Err)
  
  
  !                     type                  not output directly, filename
  ! cmfe_IN_TIMING_TYPE, cmfe_FROM_TIMING_TYPE, cmfe_ALL_TIMING_TYPE
  !CALL cmfe_TimingSetOn(cmfe_ALL_TIMING_TYPE, .TRUE.,             "", ["cmfe_Problem_Solve", "PROBLEM_SOLVE     "], Err)

  CALL cmfe_OutputSetOn("EMG",Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains = NumberOfComputationalNodes
  IF (PretendedNumberOfDomainsForDomainDecomposition == 0) THEN
    PretendedNumberOfDomainsForDomainDecomposition = NumberOfDomains
  ENDIF
  
  IF (NumberOfDomains > 0) THEN
    CALL ComputeSubdomainsWithAtoms()   ! compute domain decomposition for 3D mesh considering the atoms
  ENDIF
  GeometryIsValid = CheckGeometry()

  ! print the current directory from which the program was launched
  CALL PrintWorkingDirectory()

  ! Read in input files, stops execution if files do not exist
  CALL ReadInputFiles()

  ! output scenario information
  IF (ComputationalNodeNumber == 0) THEN
    PRINT *, "CellML file: """ // TRIM(CellMLModelFilename) // """"
    PRINT *, "Scenario name: """ // TRIM(ScenarioName) // """"
    PRINT *, ""
    PRINT *, "---------- Timing parameters -----------------------------------------------"
    PRINT *, "The time unit is 1 ms."
    PRINT "(A,F7.2,A,F5.2,A,F5.2)", "  Main loop, Δt = ", TimeStop, ", dt = ", ElasticityTimeStep
    PRINT "(A,F5.2)", "  - stimulation enabled:  Δt = ", StimDuration
    PRINT "(A,F5.2)", "  - stimulation disabled: Δt = ", (StimPeriod - StimDuration)
    PRINT *, ""

    PRINT "(A,F7.2,A,F0.5,A,I5)", "- MAIN_TIME_LOOP,         Δt =", TimeStop, ", dt = ", ElasticityTimeStep, &
      & ", N. Iter: ", CEILING(TimeStop/ElasticityTimeStep)
    
    
    IF (PdeNSteps /= -1) THEN
      PRINT "(A,F0.4,A,F0.5,2(A,I5),A)", "  - MONODOMAIN_TIME_LOOP, Δt = ", ElasticityTimeStep, ",  dt = ", PDETimeStep,&
        & ", N. Iter: ", CEILING(ElasticityTimeStep/PDETimeStep), " (set by PdeNSteps=", PdeNSteps, ")"
    ELSE
      PRINT "(A,F0.4,A,F0.5,A,I5)", "  - MONODOMAIN_TIME_LOOP, Δt = ", ElasticityTimeStep, ",  dt = ", PDETimeStep,&
        & ", N. Iter: ", CEILING(ElasticityTimeStep/PDETimeStep)
    ENDIF
      
    IF (OdeNSteps /= -1) THEN
      PRINT "(A,F0.5,A,I5,A,I5,A)", "    - SolverDAE,                       dt = ", ODETimeStep, &
        & ", N. Iter: ", CEILING(PDETimeStep/ODETimeStep), " (set by OdeNSteps=", OdeNSteps, ")"
    ELSE 
      PRINT "(A,F0.5,A,I5,A,I5,A)", "    - SolverDAE,                       dt = ", ODETimeStep, &
        & ", N. Iter: ", CEILING(PDETimeStep/ODETimeStep)
    ENDIF
    
    PRINT "(A,F0.4)", "    - SolverParabolic"
    
    IF (ElasticityDisabled) THEN
      PRINT "(A)",               "  - ELASTICITY_LOOP     (disabled) "
    ELSE
      PRINT "(A,I5)",               "  - ELASTICITY_LOOP,                                N. Iter: ",&
        & ElasticityLoopMaximumNumberOfIterations
      PRINT "(A,I4,A,E10.4)", "    - SolverFE,                 N. Iter (max): ", NewtonMaximumNumberOfIterations, &
        & ", Tol.: ",NewtonTolerance
      PRINT "(A,I4)", "      - LinearSolverFE"
      PRINT *, ""
    ENDIF
    
    IF (DebuggingOnlyRunShortPartOfSimulation) PRINT *, "Abort after first stimulation."
    PRINT *, "OutputTimeStepStride:", OutputTimeStepStride, ", EnableExportEMG:", EnableExportEMG

    ! ElasticityTimeStep should be a integer multiple of StimDuration

    ! Output problem size information
    PRINT *, ""
    PRINT *, "---------- Problem size parameters ------------------------------------------"

    PRINT "(A,3(I6,A),I12)", "# global FE-elements:      ", NumberGlobalXElements, ", ", NumberGlobalYElements, ", ", &
      & NumberGlobalZElements, &
      & ", Total: ", NumberOfElementsFE
    PRINT "(A,3(I6,A),I12)", "# local nodes per element: ", NumberOfNodesInXi1, ", ", NumberOfNodesInXi2, ", ", NumberOfNodesInXi3,&
      & ", Total: ", NumberOfNodesInXi1*NumberOfNodesInXi2*NumberOfNodesInXi3
    !PRINT "(A,I8)", "                  NumberOfInSeriesFibres: ", NumberOfInSeriesFibres
    PRINT "(A,I8)", "      NumberOfFibreLinesPerGlobalElement: ", NumberOfFibreLinesPerGlobalElement
    PRINT "(A,I8)", "              NumberOfGlobalElementLines: ", NumberOfGlobalElementLines
    !PRINT "(A,I8)", "                 NumberOfFibreLinesTotal: ", NumberOfFibreLinesTotal
    PRINT "(A,I8)", "                          NumberOfFibres: ", NumberOfFibres
    PRINT "(A,I8)", "               NumberOfElementsMPerFibre: ", NumberOfElementsMPerFibre
    PRINT "(A,I8)", "               NumberOfNodesPerLongFibre: ", NumberOfNodesPerLongFibre
    !PRINT "(A,I8)", "           NumberOfElementsMPerFibreLine: ", NumberOfElementsMPerFibreLine
    PRINT "(A,I8)", "                       NumberOfElementsM: ", NumberOfElementsM
    !PRINT "(A,I8)", "              NumberOfNodesPerShortFibre: ", NumberOfNodesPerShortFibre
    !PRINT "(A,I8)", "              NumberOfNodesMPerFibreLine: ", NumberOfNodesMPerFibreLine
    PRINT "(A,I8)", "                          NumberOfNodesM: ", NumberOfNodesM
    PRINT *,""
    PRINT "(3(A,I5),A)", "              Non-decomposable atom: ", NumberOfElementsInAtomX, &
      & "x", NumberOfElementsInAtomY, "x", NumberOfElementsInAtomZ, " elements"
    PRINT "(A,3(I5,A))", "                        Subdomain e: ", NumberOfElementsInSubdomainX, "x", &
      & NumberOfElementsInSubdomainY, "x", &
      & NumberOfElementsInSubdomainZ, " elements"
    PRINT "(A,4(I5,A))", "             Domain decomposition p: ", NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", &
      & NumberOfSubdomainsZ, &
      & " subdomains"
    PRINT "(2(A,I5))", " Number of initial domains for domain decomposition: ", PretendedNumberOfDomainsForDomainDecomposition, &
      & ", number of processes: ", NumberOfDomains
    PRINT "(A,I5)", " Number of different processes for a fibre: ", NumberOfSubdomainsX
    
    PRINT *, ""
    PRINT *, "---------- Physical parameters -----------------------------------------------"
    PRINT "(4(A,F5.2))", "       Dimensions [cm]: ",PhysicalLength,"x",PhysicalWidth,"x",PhysicalHeight,&
     & ",          InitialStretch: ", InitialStretch
    IF (NumberStimulatedNodesPerFibre == 1) THEN
      PRINT "(A,F11.2,A,F8.5,A,I3,A,F9.2,A)", " Stimulation [uA/cm^2]: ",StimValue," on ",PhysicalStimulationLength," cm, i.e. ",&
        & NumberStimulatedNodesPerFibre," node"
    ELSE
      PRINT "(A,F11.2,A,F8.5,A,I3,A,F9.2,A)", " Stimulation [uA/cm^2]: ",StimValue," on ",PhysicalStimulationLength," cm, i.e. ",&
        & NumberStimulatedNodesPerFibre," nodes, ",StimValuePerNode, " per node"
    ENDIF
    PRINT "(3(A,F5.2))", " PMax:", PMax, ",        VMax: ", VMax, ",   Conductivity: ", Conductivity
    PRINT "(3(A,F7.2))", "   Am:", Am, ", Cm (fast):", CmFast, ",     Cm (slow):", CmSlow
    
    IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
      PRINT *, "ModelType: 0 (""3a"", MultiPhysStrain, Old mechanics formulation that works in parallel.)"
    ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
      PRINT *, "ModelType: 1 (""3"", MultiPhysStrain)"
    ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
      PRINT *, "ModelType: 2 (""4"", Titin)"
    ENDIF

    PRINT *, ""
    PRINT *, "---------- Solvers -----------------------------------------------------------"
  
  
    SELECT CASE(SplittingType)
      CASE(0)
        PRINT *, "Splitting:       0 Godunov"
      CASE(1)
        PRINT *, "Splitting:       1 Strang"
      CASE DEFAULT
        PRINT *, "Splitting:       0 Godunov (default)"
    END SELECT
  
    SELECT CASE(ODESolverId)
      CASE(1)
        PRINT *, "(0D) ODE:        1 Explicit Euler"
      CASE(2)
        PRINT *, "(0D) ODE:        2 BDF"
        PRINT *, "                   Absolute Tolerance: ",DAEAbsoluteTolerance
        PRINT *, "                   Relative Tolerance: ",DAERelativeTolerance
      CASE(3)
        PRINT *, "(0D) ODE:        3 General Linear (GL)"
      CASE(4)
        PRINT *, "(0D) ODE:        4 Crank-Nicolson"
      CASE(5)
        PRINT *, "(0D) ODE:        5 Improved Euler (Heun's method)"
      CASE DEFAULT
        PRINT *, "(0D) ODE:        Explicit Euler (default)" 
    END SELECT

    SELECT CASE(MonodomainSolverId)
      CASE(1)
        PRINT *, "(1D) Monodomain: 1 SOLVER_DIRECT_LU"
      CASE(2)
        PRINT *, "(1D) Monodomain: 2 SOLVER_ITERATIVE_GMRES"
      CASE(3)
        PRINT *, "(1D) Monodomain: 3 SOLVER_ITERATIVE_CONJUGATE_GRADIENT"
      CASE(4)
        PRINT *, "(1D) Monodomain: 4 SOLVER_ITERATIVE_CONJGRAD_SQUARED"
      CASE DEFAULT
        PRINT *, "(1D) Monodomain: SOLVER_ITERATIVE_GMRES (default)" 
    END SELECT
    
    SELECT CASE(MonodomainPreconditionerId)
      CASE(1)
        PRINT *, "                 1 NO_PRECONDITIONER"
      CASE(2)
        PRINT *, "                 2 JACOBI_PRECONDITIONER"
      CASE(3)
        PRINT *, "                 3 BLOCK_JACOBI_PRECONDITIONER"
      CASE(4)
        PRINT *, "                 4 SOR_PRECONDITIONER"
      CASE(5)
        PRINT *, "                 5 INCOMPLETE_CHOLESKY_PRECONDITIONER"
      CASE(6)
        PRINT *, "                 6 INCOMPLETE_LU_PRECONDITIONER"
      CASE(7)
        PRINT *, "                 7 ADDITIVE_SCHWARZ_PRECONDITIONER"
      CASE DEFAULT
        PRINT *, "                 NO_PRECONDITIONER (default)"
    END SELECT
    PRINT *, "------------------------------------------------------------------------------"
    PRINT *, ""

    ! Output some static (compile-time) settings
    CALL cmfe_CustomProfilingGetEnabled(CustomProfilingEnabled, TAUProfilingEnabled, Err)

    IF (CustomProfilingEnabled) THEN
      PRINT*, "Custom Profiling is enabled."
    ELSE
      PRINT*, "Custom Profiling is disabled. (Enable with -DUSE_CUSTOM_PROFILING)"
    ENDIF

    IF (TAUProfilingEnabled) THEN
      PRINT*, "TAU Profiling is enabled."
    ELSE
      PRINT*, "TAU Profiling is disabled."
    ENDIF

    IF (.NOT. GeometryIsValid) THEN
      PRINT*, "Parametrization is invalid. Aborting."
    ENDIF
  ENDIF

  IF (.NOT. GeometryIsValid) THEN
    STOP
  ENDIF

  IF (ODESolverId/=5 .AND. SplittingType==1) THEN
    PRINT *, "Strang-Splitting must be used with Improved Euler method! &
     & Use ODESolverId=5 instead."
    STOP
  ENDIF
  
END SUBROUTINE ParseParameters

SUBROUTINE PrintWorkingDirectory()

  IF (ComputationalNodeNumber == 0) THEN
    CALL SYSTEM("pwd > pwd.txt", Stat)
    IF (Stat == 0) THEN
      OPEN(UNIT=100, FILE="pwd.txt", ACTION="read", IOSTAT=Stat)
      IF (Stat == 0) THEN
        READ(100,"(A)", IOSTAT=Stat) WorkingDirectory
        CLOSE(UNIT=100)
        IF (Stat == 0) THEN
          PRINT*, "Working Directory: """ // TRIM(WorkingDirectory) // """."
        ELSE
          PRINT*, "Error reading pwd.txt"
        ENDIF
      ELSE
        PRINT*, "Error opening pwd.txt"
      ENDIF
    ELSE
      PRINT*, "Error calling 'pwd'"
    ENDIF
  ENDIF
END SUBROUTINE PrintWorkingDirectory
  
SUBROUTINE CreateRegionMesh()
  INTEGER(CMISSIntg) :: ElementInFibreNo
  INTEGER(CMISSIntg) :: FibreNo
  INTEGER(CMISSIntg) :: NodeUserNumber
  INTEGER(CMISSIntg) :: ElementMGlobalNumber
  !-------------------------------------------------------------------------------------------------------------------------------
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  
  IF (ComputationalNodeNumber == 0) THEN
    IF (NumberOfComputationalNodes == 1) THEN
      PRINT "(A)", "Running with 1 process."
    ELSE 
      PRINT "(A,I6,A)", "Running with ",NumberOfComputationalNodes," processes."
    ENDIF
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemFE,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberFE,CoordinateSystemFE,Err)
  !Set the coordinate system to be 3D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemFE,3,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemFE,Err)


  ! CREATE A 1D COORDINATE SYSTEM
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemM,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberM,CoordinateSystemM,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemM,1,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL cmfe_Region_Initialise(RegionFE,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberFE,WorldRegion,RegionFE,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionFE,CoordinateSystemFE,Err)
  CALL cmfe_Region_LabelSet(RegionFE,"Region3D",Err)
  CALL cmfe_Region_CreateFinish(RegionFE,Err)


  ! CREATE A SECOND REGION
  CALL cmfe_Region_Initialise(RegionM,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberM,WorldRegion,RegionM,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionM,CoordinateSystemM,Err)
  CALL cmfe_Region_LabelSet(RegionM,"Region1D",Err)
  CALL cmfe_Region_CreateFinish(RegionM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  ! CREATE A SECOND LINEAR BASIS FOR THE 1D GRID
  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasisM,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumberM,LinearBasisM,Err)
  CALL cmfe_Basis_TypeSet(LinearBasisM,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasisM,1,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasisM,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasisM,[2],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasisM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a mesh with 8 three-dimensional elements
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMeshFE,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,RegionFE,GeneratedMeshFE,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMeshFE,[QuadraticBasis,LinearBasis],Err)
  !Define the mesh on the region

  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMeshFE,[PhysicalLength,PhysicalWidth,PhysicalHeight],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMeshFE,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(MeshFE,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMeshFE,MeshUserNumberFE,MeshFE,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! create the monodomain mesh
  !Create a mesh in the region
    
  CALL cmfe_Mesh_Initialise(MeshM,Err)
  
  !                          user_number,     region,  n_dim, mesh(out)
  CALL cmfe_Mesh_CreateStart(MeshUserNumberM, RegionM, 1,     MeshM, Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(MeshM,1,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(MeshM,NumberOfElementsM,Err)

  CALL cmfe_MeshElements_Initialise(ElementsM,Err)   
  !                                  mesh,  meshComponentNumber,           basis,        meshElements
  CALL cmfe_MeshElements_CreateStart(MeshM, MonodomainMeshComponentNumber, LinearBasisM, ElementsM, Err)
  
  ! Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(RegionM, NumberOfNodesM, Nodes, Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  
  ! Set adjacent nodes for each element
  NodeUserNumber = 1
  ElementMGlobalNumber = 1
  DO FibreNo = 1, NumberOfFibres
    !PRINT *, "Fibre ", FibreNo
    DO ElementInFibreNo = 1,NumberOfElementsMPerFibre
      
      !                               meshElements, globalElementNumber,  elementUserNodes (user numbers)
      CALL cmfe_MeshElements_NodesSet(ElementsM,    ElementMGlobalNumber, [NodeUserNumber,NodeUserNumber+1], Err)

      !PRINT "(A,I3.3,A,I3.3,A,I5.5,A,I7.7,A,I7.7)", &
      !  & "a ", ComputationalNodeNumber, ": fibre ", FibreNo, " 1D el. no. ", ElementInFibreNo, " has nodes ", &
      !  & NodeUserNumber, ", ", NodeUserNumber+1

      ! If at the end of a fibre line, increment node to starting node of next fibre line
      IF (MOD(ElementInFibreNo, NumberOfElementsMPerFibreLine) == 0) THEN
        NodeUserNumber = NodeUserNumber+1
      ENDIF
     
      NodeUserNumber = NodeUserNumber+1
      ElementMGlobalNumber = ElementMGlobalNumber+1
    ENDDO
  ENDDO
  !write(*,*) "Finished setting up 1D elements"

  CALL cmfe_MeshElements_CreateFinish(ElementsM, Err)
  CALL cmfe_Mesh_CreateFinish(MeshM, Err)
                 
  !PRINT*, "internal elements ", ELEMENTS_MAPPING%INTERNAL_START,"to",ELEMENTS_MAPPING%INTERNAL_FINISH
  !PRINT*, "element parameters:"
  
  !PRINT*, "index        element_no      interpolation_parameters"
  !DO element_idx=ELEMENTS_MAPPING%INTERNAL_START, ELEMENTS_MAPPING%INTERNAL_FINISH
  
  !  ne = ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
    
    ! get interpolation parameters of element
    ! version which is used with real preallocated variable names:
  !  CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,&
  !    & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

  !  INTERPOLATION_PARAMETERS=> &
  !    & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%INTERPOLATION_PARAMETERS                
    
    ! direct version
  !  CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    
  !  DO component_idx = 1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
  !    PRINT*, element_idx, ne, INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx)
  !  ENDDO
  !ENDDO
  


END SUBROUTINE CreateRegionMesh

! return ceil(X/Y) for integers
FUNCTION CeilDiv(X,Y)
  INTEGER(CMISSIntg), INTENT(IN) :: X, Y
  REAL(CMISSRP) :: CeilDiv
  IF (MOD(X,Y) == 0) THEN
    CeilDiv = X/Y
  ELSE
    CeilDiv = X/Y + 1
  ENDIF
END FUNCTION CeilDiv

FUNCTION GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, &
  & ParameterNumberOfSubdomainsX, ParameterNumberOfSubdomainsY, ParameterNumberOfSubdomainsZ)
  INTEGER(CMISSIntg), INTENT(IN) :: NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ
  INTEGER(CMISSIntg), INTENT(IN) :: ParameterNumberOfSubdomainsX, ParameterNumberOfSubdomainsY, ParameterNumberOfSubdomainsZ
  INTEGER(CMISSIntg) :: NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ
  INTEGER(CMISSIntg) :: nEmptySubdomainsX, nEmptySubdomainsY, nEmptySubdomainsZ, nEmptySubdomains
  INTEGER(CMISSIntg) :: nSubdomains
  INTEGER(CMISSIntg) :: GetNumberOfUsedSubdomains, NumberOfAtomsPerSubdomain
  
  IF (ParameterNumberOfSubdomainsX <= 0 .OR. ParameterNumberOfSubdomainsY <= 0 .OR. ParameterNumberOfSubdomainsZ <= 0) THEN
    GetNumberOfUsedSubdomains = HUGE(GetNumberOfUsedSubdomains)
    RETURN
  ENDIF
  
  NumberOfAtomsPerSubdomainX = CEILING(DBLE(NumberOfAtomsX) / ParameterNumberOfSubdomainsX)
  NumberOfAtomsPerSubdomainY = CEILING(DBLE(NumberOfAtomsY) / ParameterNumberOfSubdomainsY)
  NumberOfAtomsPerSubdomainZ = CEILING(DBLE(NumberOfAtomsZ) / ParameterNumberOfSubdomainsZ)
  NumberOfAtomsPerSubdomain = NumberOfAtomsPerSubdomainX*NumberOfAtomsPerSubdomainY*NumberOfAtomsPerSubdomainX

  ! decrease number of subdomains to exclude now empty subdomains
  nEmptySubdomainsX = FLOOR(DBLE(NumberOfAtomsPerSubdomainX*ParameterNumberOfSubdomainsX - NumberOfAtomsX) &
    & / NumberOfAtomsPerSubdomainX)
  nEmptySubdomainsY = FLOOR(DBLE(NumberOfAtomsPerSubdomainY*ParameterNumberOfSubdomainsY - NumberOfAtomsY) &
    & / NumberOfAtomsPerSubdomainY)
  nEmptySubdomainsZ = FLOOR(DBLE(NumberOfAtomsPerSubdomainZ*ParameterNumberOfSubdomainsZ - NumberOfAtomsZ) &
    & / NumberOfAtomsPerSubdomainZ)

  nEmptySubdomains = ParameterNumberOfSubdomainsX*ParameterNumberOfSubdomainsY*ParameterNumberOfSubdomainsZ &
    & - (ParameterNumberOfSubdomainsX-nEmptySubdomainsX)*(ParameterNumberOfSubdomainsY-nEmptySubdomainsY) &
    & * (ParameterNumberOfSubdomainsZ-nEmptySubdomainsZ)

  NumberOfSubdomainsX = ParameterNumberOfSubdomainsX - nEmptySubdomainsX
  NumberOfSubdomainsY = ParameterNumberOfSubdomainsY - nEmptySubdomainsY
  NumberOfSubdomainsZ = ParameterNumberOfSubdomainsZ - nEmptySubdomainsZ
  nSubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  GetNumberOfUsedSubdomains = nSubdomains
END FUNCTION GetNumberOfUsedSubdomains
  
SUBROUTINE ComputeSubdomainsWithAtoms()
  REAL(CMISSDP) :: OptimalSideLength
  INTEGER(CMISSIntg) :: NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfAtoms
  REAL(CMISSDP) :: NumberOfSubdomainsXFloat, NumberOfSubdomainsYFloat, NumberOfSubdomainsZFloat
  REAL(CMISSDP) :: DiffX, DiffY, DiffZ
  INTEGER(CMISSIntg) :: NumberOfAtomsPerSubdomain, nUnusedSubdomains
  INTEGER(CMISSIntg) :: nEmptySubdomainsX, nEmptySubdomainsY, nEmptySubdomainsZ, nEmptySubdomains
  INTEGER(CMISSIntg) :: nSubdomains
  INTEGER(CMISSIntg) :: NormalNumberOfElements, actualNumberOfElements
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXDecreased, DiffNumberOfDomainsYDecreased, DiffNumberOfDomainsZDecreased
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXYDecreased, DiffNumberOfDomainsXZDecreased, DiffNumberOfDomainsYZDecreased
  INTEGER(CMISSIntg) :: DiffNumberOfDomainsXYZDecreased, MinDiffNumberOfDomains
  INTEGER(CMISSIntg) :: PretendedNumberOfDomains     !< this is a value for the number of domains to be used in domain decomposition. It can be higher than the actual number of processes, because sometimes domain decomposition output produces less domains than requested.
  LOGICAL :: DEBUGGING = .FALSE.
  
  IF (DEBUGGING) DEBUGGING = (ComputationalNodeNumber == 0)  
  
  ! Test if a valid domain decomposition was specified via parameters NumberOfElementsInSubdomainXYZ and NumberOfSubdomainsXYZ
  ! if yes, use this domain decomposition and return, if not, control proceeds to computing the domain decomposition from paraameters of "atom" sizes
  IF (NumberOfElementsInSubdomainX /= 0 .AND. NumberOfElementsInSubdomainY /= 0 .AND. NumberOfElementsInSubdomainZ /= 0) THEN
    
    ! if specified number of subdomains does not match number of processes
    IF (NumberOfSubdomainsX * NumberOfSubdomainsY * NumberOfSubdomainsZ /= NumberOfComputationalNodes) THEN
      
      IF (ComputationalNodeNumber == 0) THEN
        PRINT *, "The specified number of subdomains ",NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", NumberOfSubdomainsZ, & 
          & " = ", NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ, " does not match the number of processes (", &
          & NumberOfComputationalNodes,"). Now the domain decomposition is computed from atom sizes."
      ENDIF 
    ELSE
      ! number of subdomains does match number of processes, take the specified domain decomposition

      IF (NumberOfElementsInSubdomainX*NumberOfElementsInSubdomainY*NumberOfElementsInSubdomainZ*NumberOfComputationalNodes &
        & /= NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements) THEN

        IF (ComputationalNodeNumber == 0) THEN
          PRINT *, "The specified number of 3D elements per subdomain e=",NumberOfElementsInSubdomainX, "x", &
            & NumberOfElementsInSubdomainY, "x", NumberOfElementsInSubdomainZ, " * number of subdomains (", &
            & NumberOfComputationalNodes,"), p=",NumberOfSubdomainsX, "x", NumberOfSubdomainsY, "x", NumberOfSubdomainsZ, &
            & " does not match the total number of 3D elements (", &
            & NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements, &
            & "). Now the domain decomposition is computed from atom sizes."
        ENDIF 
      ELSE
        ! The following global variables are used later in the code to assign elements to domains and therefore need to be set here:
        ! NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ

        ! For this approach without given "atom" sizes, set the atom size to be of size 1 element

        ! NumberOfSubdomainsXYZ: already set as parameter and checked, that it is valid

        ! NumberOfAtomsPerSubdomainXYZ
        NumberOfAtomsPerSubdomainX = NumberOfElementsInSubdomainX
        NumberOfAtomsPerSubdomainY = NumberOfElementsInSubdomainY
        NumberOfAtomsPerSubdomainZ = NumberOfElementsInSubdomainZ

        ! NumberOfAtomsLastSubdomainXYZ
        NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX
        NumberOfAtomsLastSubdomainY = NumberOfAtomsPerSubdomainY
        NumberOfAtomsLastSubdomainZ = NumberOfAtomsPerSubdomainZ

        ! NumberOfElementsInAtomXYZ
        ! set atom sizes to 1
        NumberOfElementsInAtomX = 1
        NumberOfElementsInAtomY = 1
        NumberOfElementsInAtomZ = 1

        ! NumberOfElementsLastAtomXYZ
        NumberOfElementsLastAtomX = 1
        NumberOfElementsLastAtomY = 1
        NumberOfElementsLastAtomZ = 1

        IF (ComputationalNodeNumber == 0) THEN
          PRINT *, "Using given domain decomposition."
        ENDIF

        RETURN
      ENDIF
    ENDIF

  ENDIF


  ! The following global variables are used later in the code to assign elements to domains and therefore need to be set here:
  ! NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ

  ! Domain decomposition for the 3D finite elasticity mesh is computed as follows: 
  ! The user can specify the size of an 'atomic' cuboid, e.g. with ax*ay*az elements. These 'atoms' then must not be split to multiple processors.
  ! In that way it is possibile to forbid the subdivision of a 1D fibre mesh to multiple processes, or to restrict it to only 2 processes per fibre.
  ! The decomposition is then done such that the cuts are 2D plains, in a way that the subdomains are as cube-shaped as possible 
  ! and as evenly distributed as possible while fulfilling the constraint.
  ! This means that not all processes may get the same amount of elements and that some processes will be left idle.
  ! This is intended and makes sense for scaling measurements:
  ! E.g. the case of 6x6x1 elements total with 10 processes. It will give 2x2x1 portions to each of the 9 processes, not using 1 process. This is better than using all 10 processes.

  PretendedNumberOfDomains = PretendedNumberOfDomainsForDomainDecomposition
  
  ! ----------  
  ! An atom is a cuboid of NumberOfElementsInAtomX x NumberOfElementsInAtomY x NumberOfElementsInAtomZ 3D finite elasticity elements that will not be distributed to multiple processes.
  ! So this is an undividable unit for domain decomposition.
  ! Determine how many atoms are possible in each direction, fractions at the end count as full atom
  NumberOfAtomsX = CeilDiv(NumberGlobalXElements, NumberOfElementsInAtomX)
  NumberOfAtomsY = CeilDiv(NumberGlobalYElements, NumberOfElementsInAtomY)
  NumberOfAtomsZ = CeilDiv(NumberGlobalZElements, NumberOfElementsInAtomZ)
  NumberOfAtoms = NumberOfAtomsX * NumberOfAtomsY * NumberOfAtomsZ

  IF (DEBUGGING) PRINT *, "PretendedNumberOfDomains=",PretendedNumberOfDomains
  IF (DEBUGGING) PRINT *, "NumberOfAtoms: ",NumberOfAtomsX,NumberOfAtomsY,NumberOfAtomsZ,"=",NumberOfAtoms

  ! subdomain = all the atoms that belong to one process
  ! determine number of subdomains in each direction, such that domains are equally sized

  ! NumberOfSubdomainsXFloat = NumberGlobalXElements / OptimalSideLength
  ! NumberOfSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
  ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
  ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat 
  !         = NumberGlobalXElements * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**3
  !         = NumberOfElementsFE / OptimalSideLength**3
  ! => OptimalSideLength = (NumberOfElementsFE / nSubdomains)**(1./3)

  OptimalSideLength = (DBLE(NumberOfElementsFE) / PretendedNumberOfDomains)**(1./3)

  IF (DEBUGGING) PRINT *, "OptimalSideLength:", OptimalSideLength,", avg. number elements z: ", &
    & DBLE(NumberGlobalZElements) / OptimalSideLength, ">", NumberOfAtomsZ,"?"
  
  IF (DBLE(NumberGlobalZElements) / OptimalSideLength > NumberOfAtomsZ - 1e-3) THEN
      
    NumberOfSubdomainsZFloat = MIN(DBLE(NumberGlobalZElements) / OptimalSideLength, DBLE(NumberOfAtomsZ))
    
    IF (DEBUGGING) PRINT *, "YES, NumberOfAtomsZ =", NumberOfAtomsZ, ", OptimalSideLength =",OptimalSideLength,", z=",&
       & NumberGlobalZElements, ", z/Opt=", DBLE(NumberGlobalZElements) / OptimalSideLength, ", NumberOfSubdomainsZFloat=", &
       & NumberOfSubdomainsZFloat

    OptimalSideLength = SQRT(NumberOfSubdomainsZFloat * NumberGlobalXElements * NumberGlobalYElements / PretendedNumberOfDomains)
    IF (DEBUGGING) PRINT *, "new OptimalSideLength=", OptimalSideLength
    
    NumberOfSubdomainsYFloat = MIN(DBLE(NumberGlobalYElements) / OptimalSideLength, DBLE(NumberOfAtomsY))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsY=", NumberOfAtomsY, ", OptimalSideLength=",OptimalSideLength,", y=",&
      & NumberGlobalYElements, ", y/Opt=", DBLE(NumberGlobalYElements) / OptimalSideLength, ", NumberOfSubdomainsYFloat=", &
      & NumberOfSubdomainsYFloat

    NumberOfSubdomainsXFloat = MIN(DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsZFloat * NumberOfSubdomainsYFloat), &
      & DBLE(NumberOfAtomsX))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsX=", NumberOfAtomsX, ", NumberOfSubdomainsXFloat=", &
       & DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsZFloat * NumberOfSubdomainsYFloat), ", final: ", &
       & NumberOfSubdomainsXFloat
    
  ! begin partioning in x direction
  ELSE
    
    NumberOfSubdomainsXFloat = MIN(DBLE(NumberGlobalXElements) / OptimalSideLength, DBLE(NumberOfAtomsX))
    ! now IF nAtomX is smaller than the optimal NumberOfSubdomainsXFloat, we must take nAtomX instead
    ! NumberOfSubdomainsXFloat is given
    ! NumberOfSubdomainsYFloat = NumberGlobalYElements / OptimalSideLength
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
    ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat 
    !         = NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / (OptimalSideLength)**2
    ! => OptimalSideLength = (NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / nSubdomains)**(1./2)
    ! 

    IF (DEBUGGING) PRINT *, "NO, NumberOfAtomsX =", NumberOfAtomsX, ", OptimalSideLength =",OptimalSideLength,", x=",&
       & NumberGlobalXElements, ", x/Opt=", DBLE(NumberGlobalXElements) / OptimalSideLength, ", NumberOfSubdomainsXFloat=", &
       & NumberOfSubdomainsXFloat

    OptimalSideLength = SQRT(NumberOfSubdomainsXFloat * NumberGlobalYElements * NumberGlobalZElements / PretendedNumberOfDomains)
    IF (DEBUGGING) PRINT *, "new OptimalSideLength=", OptimalSideLength, "=sqrt(",NumberOfSubdomainsXFloat," * ",&
      & NumberGlobalYElements," * ",NumberGlobalZElements,")/",PretendedNumberOfDomains

    NumberOfSubdomainsYFloat = MIN(DBLE(NumberGlobalYElements) / OptimalSideLength, DBLE(NumberOfAtomsY))
    IF (DEBUGGING) PRINT *, "NumberOfAtomsY=", NumberOfAtomsY, ", OptimalSideLength=",OptimalSideLength,", y=",&
      & NumberGlobalYElements, ", y/Opt=", DBLE(NumberGlobalYElements) / OptimalSideLength, ", NumberOfSubdomainsYFloat=", &
      & NumberOfSubdomainsYFloat

    ! now IF nAtomY is smaller than the optimal NumberOfSubdomainsYFloat, take NumberOfAtomsY
    ! NumberOfSubdomainsXFloat .AND. NumberOfSubdomainsYFloat are given
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / OptimalSideLength
    ! nSubdomains = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberOfSubdomainsZFloat
    !         = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / OptimalSideLength
    ! => OptimalSideLength = NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / nSubdomains
    ! NumberOfSubdomainsZFloat = NumberGlobalZElements / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements / nSubdomains)
    !               = (NumberGlobalZElements * nSubdomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat * NumberGlobalZElements)
    !               = nSubdomains / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat)
    !
     
    NumberOfSubdomainsZFloat = MIN(DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat), &
      & DBLE(NumberOfAtomsZ))

    IF (DEBUGGING) PRINT *, "NumberOfAtomsZ=", NumberOfAtomsZ, ", NumberOfSubdomainsZFloat=", &
       & DBLE(PretendedNumberOfDomains) / (NumberOfSubdomainsXFloat * NumberOfSubdomainsYFloat), ", final: ", &
       & NumberOfSubdomainsZFloat
  ENDIF
       
  IF (DEBUGGING) PRINT *, "nSubdomainsFloat=", NumberOfSubdomainsXFloat,NumberOfSubdomainsYFloat,NumberOfSubdomainsZFloat,"=",&
    & NumberOfSubdomainsXFloat*NumberOfSubdomainsYFloat*NumberOfSubdomainsZFloat

  NumberOfSubdomainsX = MAX(1, CEILING(NumberOfSubdomainsXFloat))     ! try to round up to use the maximum number of processes possible, if value gets to high, it will be decreased anyway
  NumberOfSubdomainsY = MAX(1, CEILING(NumberOfSubdomainsYFloat))
  NumberOfSubdomainsZ = MAX(1, CEILING(NumberOfSubdomainsZFloat))

  IF (DEBUGGING) PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
    & NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  ! adjust number of subdomains such that total number is <= number of domains (ideally '=')
  
  IF (GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY, &
     & NumberOfSubdomainsZ) > PretendedNumberOfDomains) THEN
    DiffNumberOfDomainsXDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsYDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsXYDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ)
    DiffNumberOfDomainsXZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsYZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ-1)
    DiffNumberOfDomainsXYZDecreased = PretendedNumberOfDomains - &
      & GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX-1, NumberOfSubdomainsY-1, &
      & NumberOfSubdomainsZ-1)

    IF (DiffNumberOfDomainsXDecreased < 0) DiffNumberOfDomainsXDecreased = HUGE(DiffNumberOfDomainsXDecreased)
    IF (DiffNumberOfDomainsYDecreased < 0) DiffNumberOfDomainsYDecreased = HUGE(DiffNumberOfDomainsYDecreased)
    IF (DiffNumberOfDomainsZDecreased < 0) DiffNumberOfDomainsZDecreased = HUGE(DiffNumberOfDomainsZDecreased)
    IF (DiffNumberOfDomainsXYDecreased < 0) DiffNumberOfDomainsXYDecreased = HUGE(DiffNumberOfDomainsXYDecreased)
    IF (DiffNumberOfDomainsXZDecreased < 0) DiffNumberOfDomainsXZDecreased = HUGE(DiffNumberOfDomainsXZDecreased)
    IF (DiffNumberOfDomainsYZDecreased < 0) DiffNumberOfDomainsYZDecreased = HUGE(DiffNumberOfDomainsYZDecreased)
    IF (DiffNumberOfDomainsXYZDecreased < 0) DiffNumberOfDomainsXYZDecreased = HUGE(DiffNumberOfDomainsXYZDecreased)
    
    MinDiffNumberOfDomains = MIN(DiffNumberOfDomainsXDecreased, MIN(DiffNumberOfDomainsYDecreased, &
      & MIN(DiffNumberOfDomainsZDecreased, MIN(DiffNumberOfDomainsXYDecreased, MIN(DiffNumberOfDomainsXZDecreased, &
      & MIN(DiffNumberOfDomainsYZDecreased, DiffNumberOfDomainsXYZDecreased))))))
      
    IF (DEBUGGING) PRINT *, "diffs: x:",DiffNumberOfDomainsXDecreased, ", y:",DiffNumberOfDomainsYDecreased, &
      & ", z:",DiffNumberOfDomainsZDecreased, ", xy:",DiffNumberOfDomainsXYDecreased, ", xz:",DiffNumberOfDomainsXZDecreased, &
      & ", yz:", DiffNumberOfDomainsYZDecreased, "xyz:",DiffNumberOfDomainsXYZDecreased

      
    IF (MinDiffNumberOfDomains == DiffNumberOfDomainsXDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) THEN
      IF (DEBUGGING) PRINT *, "best to decrease X by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsYDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Y by 1"
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Z by 1"
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXYDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease X and Y by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease X and Z by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsYZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains)) &
    & THEN
      IF (DEBUGGING) PRINT *, "best to decrease Y and Z by 1"
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSEIF (MinDiffNumberOfDomains == DiffNumberOfDomainsXYZDecreased .AND. MinDiffNumberOfDomains /= HUGE(MinDiffNumberOfDomains & 
    )) THEN
      IF (DEBUGGING) PRINT *, "best to decrease X, Y and Z by 1"
      NumberOfSubdomainsX = NumberOfSubdomainsX-1
      NumberOfSubdomainsY = NumberOfSubdomainsY-1
      NumberOfSubdomainsZ = NumberOfSubdomainsZ-1
    ELSE
      IF (DEBUGGING) PRINT *, "it does not help to decrease X,Y or Z by 1, start iterative procedure"
      DO WHILE (GetNumberOfUsedSubdomains(NumberOfAtomsX, NumberOfAtomsY, NumberOfAtomsZ, NumberOfSubdomainsX, &
        & NumberOfSubdomainsY, NumberOfSubdomainsZ) > PretendedNumberOfDomains)
        DiffX = NumberOfSubdomainsX - NumberOfSubdomainsXFloat
        DiffY = NumberOfSubdomainsY - NumberOfSubdomainsYFloat
        DiffZ = NumberOfSubdomainsZ - NumberOfSubdomainsZFloat
        
        IF (DiffX >= DiffY .AND. DiffX >= DiffZ) THEN
          IF (NumberOfSubdomainsX /= 1) THEN
            NumberOfSubdomainsX = NumberOfSubdomainsX - 1
          ELSEIF (DiffY >= DiffZ) THEN
            IF (NumberOfSubdomainsY /= 1) THEN
              NumberOfSubdomainsY = NumberOfSubdomainsY - 1
            ELSE
              NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
            ENDIF
          ELSE
            IF (NumberOfSubdomainsZ /= 1) THEN
              NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
            ELSE
              NumberOfSubdomainsY = NumberOfSubdomainsY - 1
            ENDIF
          ENDIF
              
        ELSEIF (DiffY >= DiffZ) THEN    ! DiffY >= DiffZ, DiffY > DiffX
          IF (NumberOfSubdomainsY /= 1) THEN
            NumberOfSubdomainsY = NumberOfSubdomainsY - 1
          ELSE
            IF (DiffX >= DiffZ) THEN
              IF (NumberOfSubdomainsX /= 1) THEN
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ELSE
                NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
              ENDIF
            ELSE
              IF (NumberOfSubdomainsZ /= 1) THEN
                NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
              ELSE
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ENDIF
            ENDIF
          ENDIF
      
        ELSE       ! DiffZ > DiffY, DiffZ >= DiffX
          IF (NumberOfSubdomainsZ /= 1) THEN
            NumberOfSubdomainsZ = NumberOfSubdomainsZ - 1
          ELSE
            IF (DiffX >= DiffY) THEN
              IF (NumberOfSubdomainsX /= 1) THEN
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ELSE
                NumberOfSubdomainsY = NumberOfSubdomainsY - 1
              ENDIF
            ELSE
              IF (NumberOfSubdomainsY /= 1) THEN
                NumberOfSubdomainsY = NumberOfSubdomainsY - 1
              ELSE
                NumberOfSubdomainsX = NumberOfSubdomainsX - 1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
          
        IF (DEBUGGING) PRINT *, "Diff: ", DiffX,DiffY,DiffZ, ", nSubdomains:",NumberOfSubdomainsX, NumberOfSubdomainsY, &
          & NumberOfSubdomainsZ, "=", NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ
      ENDDO
    ENDIF  
  ENDIF
  
  
  IF (DEBUGGING) PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
    &  NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  ! determine shape of subdomain, i.e. NumberOfAtomsPerSubdomain
  ! NumberOfAtomsPerSubdomainX * NumberOfSubdomainsX = NumberOfAtomsX  =>  NumberOfAtomsPerSubdomainX = NumberOfAtomsX / NumberOfSubdomainsX

  NumberOfAtomsPerSubdomainX = CEILING(DBLE(NumberOfAtomsX) / NumberOfSubdomainsX)
  NumberOfAtomsPerSubdomainY = CEILING(DBLE(NumberOfAtomsY) / NumberOfSubdomainsY)
  NumberOfAtomsPerSubdomainZ = CEILING(DBLE(NumberOfAtomsZ) / NumberOfSubdomainsZ)
  NumberOfAtomsPerSubdomain = NumberOfAtomsPerSubdomainX*NumberOfAtomsPerSubdomainY*NumberOfAtomsPerSubdomainZ

  ! decrease number of subdomains to exclude now empty subdomains
  nEmptySubdomainsX = FLOOR(DBLE(NumberOfAtomsPerSubdomainX*NumberOfSubdomainsX - NumberOfAtomsX) / NumberOfAtomsPerSubdomainX)
  nEmptySubdomainsY = FLOOR(DBLE(NumberOfAtomsPerSubdomainY*NumberOfSubdomainsY - NumberOfAtomsY) / NumberOfAtomsPerSubdomainY)
  nEmptySubdomainsZ = FLOOR(DBLE(NumberOfAtomsPerSubdomainZ*NumberOfSubdomainsZ - NumberOfAtomsZ) / NumberOfAtomsPerSubdomainZ)

  IF (DEBUGGING) PRINT *, "NumberOfAtomsPerSubdomain: ",NumberOfAtomsPerSubdomainX,NumberOfAtomsPerSubdomainY, &
    & NumberOfAtomsPerSubdomainZ
  
  nEmptySubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ &
    & - (NumberOfSubdomainsX-nEmptySubdomainsX)*(NumberOfSubdomainsY-nEmptySubdomainsY)*(NumberOfSubdomainsZ-nEmptySubdomainsZ)

  IF (DEBUGGING) PRINT *, "nEmptySubdomains: ",nEmptySubdomainsX,nEmptySubdomainsY,nEmptySubdomainsZ,&
    & ", total: ",nEmptySubdomains

  NumberOfSubdomainsX = NumberOfSubdomainsX - nEmptySubdomainsX
  NumberOfSubdomainsY = NumberOfSubdomainsY - nEmptySubdomainsY
  NumberOfSubdomainsZ = NumberOfSubdomainsZ - nEmptySubdomainsZ
  nSubdomains = NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ

  nUnusedSubdomains = NumberOfDomains - nSubdomains
  IF (DEBUGGING) PRINT *, "PretendedNumberOfDomains: ", PretendedNumberOfDomains,", nSubdomains:", nSubdomains, &
    & ", nUnusedSubdomains:",nUnusedSubdomains
  
  IF (nUnusedSubdomains /= 0) THEN
    IF (ComputationalNodeNumber == 0) THEN
      PRINT *, "The computed decomposition contains ", nUnusedSubdomains, " subdomains less than given processes. " // &
        & NEW_LINE('A') // " This is intended and no error. Restart program with ", nSubdomains, &
        & " processes instead of ", NumberOfDomains, " or use PretendedNumberOfDomains."
    ENDIF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
    CALL MPI_FINALIZE(Err)
    STOP
  ENDIF

  NumberOfElementsInSubdomainX = NumberOfElementsInAtomX*NumberOfAtomsPerSubdomainX
  NumberOfElementsInSubdomainY = NumberOfElementsInAtomY*NumberOfAtomsPerSubdomainY
  NumberOfElementsInSubdomainZ = NumberOfElementsInAtomZ*NumberOfAtomsPerSubdomainZ

  IF (DEBUGGING) THEN
    PRINT *, "NumberOfAtomsPerSubdomain: ",NumberOfAtomsPerSubdomainX, NumberOfAtomsPerSubdomainY, NumberOfAtomsPerSubdomainZ, &
      & "=", NumberOfAtomsPerSubdomain
    PRINT *, "nSubdomains: ",NumberOfSubdomainsX, NumberOfSubdomainsY, NumberOfSubdomainsZ, "=", &
      & NumberOfSubdomainsX*NumberOfSubdomainsY*NumberOfSubdomainsZ, &
      & ", NumberOfAtomsPerSubdomain: ", NumberOfAtomsPerSubdomainX,NumberOfAtomsPerSubdomainY,NumberOfAtomsPerSubdomainZ, "=", &
      & NumberOfAtomsPerSubdomain, ", NumberOfElementsInAtom: ", NumberOfElementsInAtomX, NumberOfElementsInAtomY, &
      & NumberOfElementsInAtomZ, "=", NumberOfElementsInAtomX*NumberOfElementsInAtomY*NumberOfElementsInAtomX, & 
      ", NumberGlobalElements ", NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements,"=",&
      &  NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements

    PRINT *, "subdomain size: ",NumberOfElementsInSubdomainX,"x",&
       & NumberOfElementsInSubdomainY, "x", NumberOfElementsInSubdomainZ
    PRINT *, "total size: ", NumberOfSubdomainsX*NumberOfElementsInSubdomainX, "x",&
       & NumberOfSubdomainsY*NumberOfElementsInSubdomainY, "x", &
       & NumberOfSubdomainsZ*NumberOfElementsInSubdomainZ, & 
       " >= ", NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements
  ENDIF

  ! determine last subdomain sizes
  ! NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX - (NumberOfAtomsPerSubdomainX*NumberOfSubdomainsX - NumberOfAtomsX)
  NumberOfAtomsLastSubdomainX = NumberOfAtomsPerSubdomainX*(1-NumberOfSubdomainsX) + NumberOfAtomsX
  NumberOfAtomsLastSubdomainY = NumberOfAtomsPerSubdomainY*(1-NumberOfSubdomainsY) + NumberOfAtomsY
  NumberOfAtomsLastSubdomainZ = NumberOfAtomsPerSubdomainZ*(1-NumberOfSubdomainsZ) + NumberOfAtomsZ

  ! NumberOfElementsLastAtomX = NumberOfElementsInAtomX - (NumberOfAtomsX*NumberOfElementsInAtomX - NumberGlobalXElements)
  ! NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1 - NumberOfAtomsX) + NumberGlobalXElements
  NumberOfElementsLastAtomX = NumberOfElementsInAtomX*(1-NumberOfAtomsX) + NumberGlobalXElements
  NumberOfElementsLastAtomY = NumberOfElementsInAtomY*(1-NumberOfAtomsY) + NumberGlobalYElements
  NumberOfElementsLastAtomZ = NumberOfElementsInAtomZ*(1-NumberOfAtomsZ) + NumberGlobalZElements


  IF (DEBUGGING) PRINT *, "NumberOfAtomsLastSubdomain: ", NumberOfAtomsLastSubdomainX, NumberOfAtomsLastSubdomainY, &
    & NumberOfAtomsLastSubdomainZ, ", NumberOfElementsLastAtom: ",NumberOfElementsLastAtomX, NumberOfElementsLastAtomY, &
    & NumberOfElementsLastAtomZ

  ! for special subdomains (at x+,y+,z+ border) output the number of missing elements
  
  NormalNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ

  IF (ComputationalNodeNumber == 0) PRINT *, nSubdomains, " process(es), normal number of elements per process: ", &
    & normalNumberOfElements

  ! x+ face
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY&
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF ((NumberOfSubdomainsY-1)*(NumberOfSubdomainsZ-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsY-1)*(NumberOfSubdomainsZ-1), &
      & " process(es) on    x+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
        
  ! y+ face
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF ((NumberOfSubdomainsX-1)*(NumberOfSubdomainsZ-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsX-1)*(NumberOfSubdomainsZ-1), &
      & " process(es) on    y+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
    
  ! z+ face
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF ((NumberOfSubdomainsX-1)*(NumberOfSubdomainsY-1) > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, (NumberOfSubdomainsX-1)*(NumberOfSubdomainsY-1), &
      & " process(es) on    z+ boundary have ", normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ y+ edge
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ
        
  IF (NumberOfSubdomainsZ-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsZ-1, " process(es) on x+,y+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ z+ edge
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (NumberOfSubdomainsY-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsY-1, " process(es) on x+,z+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! y+ z+ edge
  actualNumberOfElements = &
    NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (NumberOfSubdomainsX-1 > 0 .AND. actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, NumberOfSubdomainsX-1, " process(es) on y+,z+ boundary have ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
  ! x+ y+ z+ subdomain
  actualNumberOfElements = &
    ((NumberOfAtomsLastSubdomainX-1)*NumberOfElementsInAtomX + NumberOfElementsLastAtomX) &
    *((NumberOfAtomsLastSubdomainY-1)*NumberOfElementsInAtomY + NumberOfElementsLastAtomY) &
    *((NumberOfAtomsLastSubdomainZ-1)*NumberOfElementsInAtomZ + NumberOfElementsLastAtomZ)
        
  IF (actualNumberOfElements < normalNumberOfElements) THEN
    IF (ComputationalNodeNumber == 0) PRINT *, "          1  process on  x+,y+,z+ boundary has  ", &
      & normalNumberOfElements-actualNumberOfElements," elements less (only ", &
      & 100.*actualNumberOfElements/normalNumberOfElements,"%)"
  ENDIF
  
END SUBROUTINE ComputeSubdomainsWithAtoms

SUBROUTINE CreateDecompositionFiniteElasticity()
  INTEGER(CMISSIntg) :: DomainNo
  INTEGER(CMISSIntg) :: ElementFENo
  INTEGER(CMISSIntg) :: SubdomainZ, SubdomainY, SubdomainX
  INTEGER(CMISSIntg) :: AtomX, AtomY, AtomZ
  INTEGER(CMISSIntg) :: X, Y, Z, XIndex, YIndex, ZIndex
  INTEGER(CMISSIntg) :: NumberOfAtomsCurrentSubdomainX, NumberOfAtomsCurrentSubdomainY, NumberOfAtomsCurrentSubdomainZ
  INTEGER(CMISSIntg) :: NumberOfElementsCurrentAtomX, NumberOfElementsCurrentAtomY, NumberOfElementsCurrentAtomZ
  LOGICAL :: DEBUGGING = .FALSE.
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Create a decomposition for global FE elements
  ! The following global variables are read: NumberOfSubdomainsXYZ, NumberOfAtomsPerSubdomainXYZ, NumberOfAtomsLastSubdomainXYZ, NumberOfElementsInAtomXYZ, NumberOfElementsLastAtomXYZ
  CALL cmfe_Decomposition_Initialise(DecompositionFE,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberFE,MeshFE,DecompositionFE,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(DecompositionFE,.TRUE.,Err)

  IF (DEBUGGING) DEBUGGING = (ComputationalNodeNumber == 2)
  
  IF(NumberOfDomains > 1) THEN
    CALL cmfe_Decomposition_TypeSet(DecompositionFE, CMFE_DECOMPOSITION_USER_DEFINED_TYPE, Err)
    
    ! assign domain values
        
    ! iterate over subdomains
    DO SubdomainZ = 1,NumberOfSubdomainsZ
      DO SubdomainY = 1,NumberOfSubdomainsY
        DO SubdomainX = 1,NumberOfSubdomainsX
        
          DomainNo = (SubdomainZ-1)*NumberOfSubdomainsY*NumberOfSubdomainsX + (SubdomainY-1)*NumberOfSubdomainsX + SubdomainX-1
          
          NumberOfAtomsCurrentSubdomainX = NumberOfAtomsPerSubdomainX
          IF (SubdomainX == NumberOfSubdomainsX) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainX = NumberOfAtomsLastSubdomainX
          ENDIF
            
          NumberOfAtomsCurrentSubdomainY = NumberOfAtomsPerSubdomainY          
          IF (SubdomainY == NumberOfSubdomainsY) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainY = NumberOfAtomsLastSubdomainY
          ENDIF
            
          NumberOfAtomsCurrentSubdomainZ = NumberOfAtomsPerSubdomainZ
          IF (SubdomainZ == NumberOfSubdomainsZ) THEN  ! last subdomain
            NumberOfAtomsCurrentSubdomainZ = NumberOfAtomsLastSubdomainZ
          ENDIF
            
          IF (DEBUGGING) THEN
            PRINT *, "subdomain (",SubdomainX,",",SubdomainY,",",SubdomainZ,") has ", NumberOfAtomsCurrentSubdomainX,",", &
              & NumberOfAtomsCurrentSubdomainY,",",NumberOfAtomsCurrentSubdomainZ," atoms"
          ENDIF
          
          DO AtomZ = 1,NumberOfAtomsCurrentSubdomainZ
            DO AtomY = 1,NumberOfAtomsCurrentSubdomainY
              DO AtomX = 1,NumberOfAtomsCurrentSubdomainX
                  
                NumberOfElementsCurrentAtomX = NumberOfElementsInAtomX
                IF (SubdomainX == NumberOfSubdomainsX .AND. AtomX == NumberOfAtomsCurrentSubdomainX) THEN  ! last atom
                  NumberOfElementsCurrentAtomX = NumberOfElementsLastAtomX
                ENDIF
                  
                NumberOfElementsCurrentAtomY = NumberOfElementsInAtomY
                IF (SubdomainY == NumberOfSubdomainsY .AND. AtomY == NumberOfAtomsCurrentSubdomainY) THEN  ! last atom
                  NumberOfElementsCurrentAtomY = NumberOfElementsLastAtomY
                ENDIF
                  
                NumberOfElementsCurrentAtomZ = NumberOfElementsInAtomZ
                IF (SubdomainZ == NumberOfSubdomainsZ .AND. AtomZ == NumberOfAtomsCurrentSubdomainZ) THEN  ! last atom
                  NumberOfElementsCurrentAtomZ = NumberOfElementsLastAtomZ
                ENDIF
                  
                IF (DEBUGGING) PRINT *, "  atom (",AtomX,",",AtomY,",",AtomZ,") has ",NumberOfElementsCurrentAtomX,",",&
                  & NumberOfElementsCurrentAtomY,",",NumberOfElementsCurrentAtomZ," elements"
                  
                DO Z = 1,NumberOfElementsCurrentAtomZ
                  DO Y = 1,NumberOfElementsCurrentAtomY
                    DO X = 1,NumberOfElementsCurrentAtomX
                      
                      ZIndex = (SubdomainZ-1)*NumberOfAtomsPerSubdomainZ*NumberOfElementsInAtomZ &
                        & + (AtomZ-1)*NumberOfElementsInAtomZ + Z
                      YIndex = (SubdomainY-1)*NumberOfAtomsPerSubdomainY*NumberOfElementsInAtomY &
                        & + (AtomY-1)*NumberOfElementsInAtomY + Y
                      XIndex = (SubdomainX-1)*NumberOfAtomsPerSubdomainX*NumberOfElementsInAtomX &
                        & + (AtomX-1)*NumberOfElementsInAtomX + X
                      
                      ElementFENo = (ZIndex-1)*NumberGlobalXElements*NumberGlobalYElements + (YIndex-1)*NumberGlobalXElements &
                        & + XIndex
                      
                      IF (DEBUGGING) PRINT *,"      (x,y,z)=(",XIndex,",",YIndex,",",ZIndex,") i=",ElementFENo," domain ",DomainNo
     
                      !                                        DECOMPOSITION,   GLOBAL_ELEMENT_NUMBER, DOMAIN_NUMBER
                      CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE, ElementFENo,           DomainNo, Err)
                    
                    ENDDO ! X
                  ENDDO   ! Y
                ENDDO     ! Z
              ENDDO ! AtomX
            ENDDO   ! AtomY
          ENDDO     ! AtomZ
        ENDDO ! SubdomainX
      ENDDO   ! SubdomainY
    ENDDO     ! SubdomainZ
                    
        
    IF (DEBUGGING) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
      CALL MPI_ABORT(MPI_COMM_WORLD, Err, Err)
      PRINT *, "stop in FortranExample.f90:1733"
      STOP
    ENDIF
    
  ELSE
    ! single process
    CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  ENDIF

  ! finish decomposition
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(DecompositionFE,Err)
END SUBROUTINE CreateDecompositionFiniteElasticity

SUBROUTINE CreateDecompositionMonodomain()

  INTEGER(CMISSIntg) :: DomainNo
  INTEGER(CMISSIntg) :: ElementMGlobalNumber, ElementMInFibreNo
  INTEGER(CMISSIntg) :: FibreNo
  INTEGER(CMISSIntg) :: FEElementGlobalNumber
  INTEGER(CMISSIntg) :: FEElementXIdx, FEElementYIdx, FEElementZIdx
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Create a decompositions for monodomain elements on fibres
  !PRINT "(I3.3,A)", ComputationalNodeNumber, ": Create decomposition for monodomain elements"
  
  CALL cmfe_Decomposition_Initialise(DecompositionM,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberM,MeshM,DecompositionM,Err)
         
  IF(NumberOfDomains > 1) THEN
    ! set decomposition to be of user defined type
    CALL cmfe_Decomposition_TypeSet(DecompositionM, CMFE_DECOMPOSITION_USER_DEFINED_TYPE, Err)
    
    ! assign the same domains to bioelectric nodes as the containing FE elements
    ! loop over elementsM
    ElementMGlobalNumber = 1
    DO FibreNo = 1, NumberOfFibres
      !PRINT *, "Fibre ",FibreNo
    
      DO ElementMInFibreNo = 1,NumberOfElementsMPerFibre
        
        ! compute global ElementFE that contains the current elementM
        FEElementZIdx = INT(INT((FibreNo-1) / NumberGlobalYFibres) / NumberOfNodesInXi3)
        FEElementYIdx = INT(MOD((FibreNo-1), NumberGlobalYFibres)  / NumberOfNodesInXi2)
        FEElementXIdx = (ElementMInFibreNo-1) / NumberOfNodesInXi1
        
        FEElementGlobalNumber = FEElementZIdx * NumberGlobalYElements * NumberGlobalXElements &
         & + FEElementYIdx * NumberGlobalXElements + FEElementXIdx + 1
        
        !PRINT "(I1.1,2(A,I2.2),4(A,I3.3))", ComputationalNodeNumber,": Fibre", FibreNo, ", ElementM ", ElementMGlobalNumber, &
        !  & ", Element (",FEElementXIdx,",",FEElementYIdx,",",FEElementZIdx, &
        !  & ")=", FEElementGlobalNumber
      
        ! get the domain of the global ElementFE
        !                                        DECOMPOSITION,   USER_ELEMENT_NUMBER,   DOMAIN_NUMBER
        CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE, FEElementGlobalNumber, DomainNo, Err)

        ! set the domain of the ElementM to the same domain as the containing global element
        !                                        DECOMPOSITION,   GLOBAL_ELEMENT_NUMBER, DOMAIN_NUMBER
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionM,  ElementMGlobalNumber,  DomainNo, Err)
        !PRINT "(I1.1,A,I5.5,A,I2,A,I2)", ComputationalNodeNumber, ": fibre ", FibreNo, ", 1D el. no. ", ElementMGlobalNumber, &
        !  & " to domain no. ", DomainNo
          
        ElementMGlobalNumber = ElementMGlobalNumber + 1
      ENDDO
    ENDDO
  ELSE
    ! single process
    CALL cmfe_Decomposition_TypeSet(DecompositionM,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  ENDIF

  ! finish decomposition
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
    
  ! In this call it is checked whether all processes have at least one node
  CALL cmfe_Decomposition_CreateFinish(DecompositionM,Err)

  IF (.FALSE.) THEN
    PRINT*, "Print elements mapping in FortranExample.f90:1565"
    PRINT*, "--------------- ElementsMapping --------------------------"
    CALL cmfe_PrintElementsMapping(DecompositionM,Err)
  ENDIF
  
  !PRINT*, "Print nodes mapping in FortranExample.f90:1549"
  !CALL cmfe_PrintNodesMapping(DecompositionM,Err)
  
END SUBROUTINE CreateDecompositionMonodomain

SUBROUTINE CreateFieldFiniteElasticity()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberFE,RegionFE,GeometricFieldFE,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldFE,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldFE,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldFE,Err)

!  CALL cmfe_Field_ParameterSetUpdateStart(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
!  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMeshFE,GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,RegionFE,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(MaterialFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberFE,RegionFE,MaterialFieldFE,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldFE,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL cmfe_Field_VariableTypesSet(MaterialFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,Err)
  ELSEIF (ModelType == 2) THEN  ! 4, "Titin"
    CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
  ENDIF

  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

  IF (ModelType == 1 .OR. ModelType == 2) THEN !
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF

  IF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF

  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldFE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !Set Mooney-Rivlin constants c10 and c01.
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
     & Err)
     
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable

    !Set Material-Parameters [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,C(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,C(5),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,C(6),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,C(7),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,C(8),Err)
 
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    !Set Mooney-Rivlin constants c10 and c01.
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
     & Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,TkLinParam,&
     & Err)

  ENDIF



  !Create the dependent field for FE with 2 variables and * components
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL cmfe_Field_Initialise(DependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberFE,RegionFE,DependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(DependentFieldFE,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldFE,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)

  CALL cmfe_Field_VariableTypesSet(DependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
!  CALL cmfe_Field_ScalingTypeSet(DependentFieldFE,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the mechanics mesh
  CALL cmfe_Field_Initialise(IndependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldIndependentUserNumberFE,RegionFE,IndependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(IndependentFieldFE,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(IndependentFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(IndependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(IndependentFieldFE,CMFE_FIELD_INDEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldFE,2,Err)
  CALL cmfe_Field_VariableTypesSet(IndependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_DimensionSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_DimensionSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
  ENDIF

  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,5,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
       & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5, &
      & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !lambda_a
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (unbound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (bound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (unbound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (bound)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! activation for titin
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  ENDIF

  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI1
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,2, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI2
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,3, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! number of nodes in XI3
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,4, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) ! if fibre starts in current FE element
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,5, &
    & CMFE_FIELD_CONSTANT_INTERPOLATION,Err) ! number of in series fibres
  CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_FE",Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_FE",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_FE",Err)
  ENDIF

  CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"subgrid_info",Err)
  CALL cmfe_Field_CreateFinish(IndependentFieldFE,Err)

END SUBROUTINE CreateFieldFiniteElasticity

SUBROUTINE CreateFieldMonodomain()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldM,Err)
  ! FIELD_CREATE_START(fieldUserNumber,REGION,FIELD)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberM,RegionM,GeometricFieldM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldM,DecompositionM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldM,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldM,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL cmfe_Field_Initialise(MaterialFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberM,RegionM,MaterialFieldM,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldM,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsM,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Am,Err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,CmFast,Err)
  !Set Conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & Conductivity,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 3 variables. [U: 1 component, DelUdelN: 1 component, V: 3]
  !DependentFieldM: FIELD_U_VARIABLE_TYPE: 1) Vm, FIELD_DELUDELN_VARIABLE_TYPE: 1)dVm/dn, FIELD_V_VARIABLE_TYPE: 1),2),3): GeometryM3D, 3D-position of geometry
  CALL cmfe_Field_Initialise(DependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberM,RegionM,DependentFieldM,Err)
  CALL cmfe_Field_TypeSet(DependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldM,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL cmfe_Field_VariableTypesSet(DependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
   & CMFE_FIELD_V_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
  !additional the V_Var_Type with 3 components for the 3-d position of the geometry
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM, & 
    & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
!  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"GeometryM3D",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the bioelectrics mesh
  CALL cmfe_Field_Initialise(IndependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldIndependentUserNumberM,RegionM,IndependentFieldM,Err)
  CALL cmfe_Field_TypeSet(IndependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(IndependentFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(IndependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(IndependentFieldM,CMFE_FIELD_INDEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldM,4,Err)
  CALL cmfe_Field_VariableTypesSet(IndependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE, &
   & CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_U2_VARIABLE_TYPE],Err)

  !first variable:   CMFE_FIELD_U_VARIABLE_TYPE
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !first variable:   CMFE_FIELD_U_VARIABLE_TYPE -- 1) active stress
    CALL cmfe_Field_DimensionSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_M",Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_M",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_DimensionSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !normalised sarcomere-based active stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the unbound titin stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the bound titin stress
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the unbound titin stress in the XF-direction
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,5, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the bound titin stress in the XF-direction
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the titin activation (stress without force-length relation)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Active_Stress_M",Err)
  ENDIF

  !second variable:   CMFE_FIELD_V_VARIABLE_TYPE -- 1) motor unit number   2) fibre type   3) fibre number   4) nearest Gauss point   5) containing element local number 6) in-fibre contiguous node number
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,FieldIndependentNumberOfComponentsM2,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,2, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,4, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,5, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,6, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"fibre_info",Err)

  !third variable:   FIELD_U1_VARIABLE_TYPE -- 1) half-sarcomere length   2) inital half-sarcomere length   3) initial node distance
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)

IF (ModelType == 0 .OR. ModelType == 1) THEN ! 3,3a MultiPhysStrain
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3,Err)
ELSEIF (ModelType == 2) THEN ! 4, "Titin"
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,4,Err)
ENDIF
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,2, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

  IF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF

  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,"half-sarcomere_length",Err)

  !fourth variable:   FIELD_U2_VARIABLE_TYPE -- 1) old node distance   2) maximum contraction velocity   3) relative contraction velocity   4) velocity before 1 time step   5) velocity before 2 time step   6) velocity before 3 time steps   7) node distance to right node (only computed on some nodes)
  CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,7,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,1, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,2, &
   & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,4, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,5, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,6, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,7, &
   & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,"contraction_velocity",Err)

  CALL cmfe_Field_CreateFinish(IndependentFieldM,Err)

END SUBROUTINE CreateFieldMonodomain

SUBROUTINE CreateEquationsSet()
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity

  CALL cmfe_Field_Initialise(EquationsSetFieldFE,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetFE,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], &
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], &
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
     & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], & 
     & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  ENDIF
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumberFE,IndependentFieldFE,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL cmfe_Field_Initialise(EquationsSetFieldM,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM, &
     & [CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
     & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], &
     & FieldEquationsSetUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  ENDIF
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetM,Err)

  !Create the equations set dependent field variables for monodomain
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetM,FieldIndependentUserNumberM,IndependentFieldM,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetM,Err)

  !Create the equations set materials field variables for monodomain
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetM,Err)
    
END SUBROUTINE CreateEquationsSet

SUBROUTINE InitializeFieldMonodomain()
  INTEGER(CMISSIntg) :: FibreLineNo
  INTEGER(CMISSIntg) :: FibreInLineNo
  INTEGER(CMISSIntg) :: FibreNo, NodeNo, NodeIdx
  INTEGER(CMISSIntg) :: MotorUnitRank
  INTEGER(CMISSIntg) :: FEElementGlobalNumber, FEElementLocalNumber
  INTEGER(CMISSIntg) :: FEElementXIdx, FEElementYIdx, FEElementZIdx
  INTEGER(CMISSIntg) :: NodeMUserNumber
  INTEGER(CMISSIntg) :: MElementUserNumber, ElementIdx
  INTEGER(CMISSIntg) :: ElementDomain, NodeDomain
  INTEGER(CMISSIntg), DIMENSION(2) :: ElementUserNodes
  REAL(CMISSDP) :: InitialBioelectricNodeDistance

  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  ! OldTomoMechanics (ModelType 0)
  !first variable (U)
  !  components:
  !    1) active stress
  !
  ! .NOT. OldTomoMechanics (ModelType 1)
  !first variable (U)
  !  components:
  !    1) A_1
  !    2) A_2
  !    3) x_1
  !    4) x_2
  !
  !second variable (V)
  !  components:
  !    1) motor unit number
  !    2) fibre type
  !    3) fibre number
  !    4) nearest Gauss point
  !    5) FE element local number that contains bioelectric node
  !
  !init the motor unit number to 101
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
   & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,101,Err) !

  ! initialize values
  !init the fibre type to 1
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,1,Err) !Ftype=1
  !init the fibre number, the nearest Gauss point info and the inElem info to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, &
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM, & 
    & CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0,Err)  ! local element number
  !third variable (U1):
  !  components:
  !    1) half-sarcomere length
  !    2) initial half-sarcomere length
  !    3) initial node distance
  InitialBioelectricNodeDistance = PhysicalLength / (NumberOfElementsMPerFibreLine)
  !PRINT *, "InitialBioelectricNodeDistance = ", InitialBioelectricNodeDistance
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & 1.0_CMISSRP,Err) ! lengths in the cell model are in /micro/meters!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & InitialBioelectricNodeDistance, Err)  !lengths in the cell model are in /micro/meters!!!
  !fourth variable (U2):
  !  components:
  !    1) old node distance
  !    2) maximum contraction velocity
  !    3) relative contraction velocity
  !    4) node distance to previous node before 1 time step   
  !    5) node distance to previous node before 2 time steps
  !    6) node distance to previous node before 3 time steps
  !    7) node distance to right node (only computed on some nodes)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & InitialBioelectricNodeDistance,Err)  !lengths in the cell model are in /micro/meters!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & Vmax / (NumberOfElementsMPerFibreLine),Err)    !velocity in the cell model is in micro/meters/millisecond!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)    !velocity in the cell model is in /micro/meters/per/millisecond!!!
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6, &
   & InitialBioelectricNodeDistance,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7, &
   & 0.0_CMISSRP,Err)
     
  ! Store containing FE element local numbers for bioelectric nodes in IndependentFieldM V comp. 5
  ! loop over fibres and bioelectric nodes
  DO FibreNo = 1, NumberOfFibres
    
    ! get rank of fibre
    MotorUnitRank = MUDistribution(MOD(FibreNo-1, 4050)+1)
    
    ! loop over elements of fibre, for each element, the right node is considered, except for the first element (on ElementIdx=0), then the left node of the first element is considered
    DO ElementIdx = 0, NumberOfElementsMPerFibre
      
      ! compute the bioelectric element user number
      IF (ElementIdx == 0) THEN   ! index 0 and 1 are both for the first element, but first for the left node, then for the right node
        MElementUserNumber = (FibreNo-1) * NumberOfElementsMPerFibre + 1
      ELSE
        MElementUserNumber = (FibreNo-1) * NumberOfElementsMPerFibre + ElementIdx
      ENDIF
      
      ! get the nodes of the current element
      CALL cmfe_MeshElements_NodesGet(ElementsM, MElementUserNumber, ElementUserNodes, Err)
      ! ElementUserNodes(1:2): node user number
          
      IF (ElementIdx == 0) THEN
        NodeMUserNumber = ElementUserNodes(1)   ! take left node of first element in first iteration per fibre
      ELSE
        NodeMUserNumber = ElementUserNodes(2)
      ENDIF
      
      !PRINT "(I1.1,3(A,I3))", ComputationalNodeNumber,": Fibre", FibreNo, ", Element in fibre ", ElementIdx, &
      !  & ", NodeMUserNumber:", NodeMUserNumber        
      
      ! get domain number of element and node
      !                                        decomposition,  elementUserNumber, domain
      CALL cmfe_Decomposition_ElementDomainGet(DecompositionM, MElementUserNumber, ElementDomain, Err)
      
      !                                     decomposition,  nodeUserNumber,  meshComponentNumber, domain
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM, NodeMUserNumber, 1,                   NodeDomain, Err)
      
      IF (NodeDomain == ComputationalNodeNumber) THEN
      
        ! compute number of containing FE element
        FEElementZIdx = INT(INT((FibreNo-1) / NumberGlobalYFibres) / NumberOfNodesInXi3)
        FEElementYIdx = INT(MOD((FibreNo-1), NumberGlobalYFibres)  / NumberOfNodesInXi2)
        FEElementXIdx = (ElementIdx-1) / NumberOfNodesInXi1
        IF (ElementIdx == 0) THEN
          FEElementXIdx = 0
        ENDIF
      
        FEElementGlobalNumber = FEElementZIdx * NumberGlobalYElements * NumberGlobalXElements &
         & + FEElementYIdx * NumberGlobalXElements + FEElementXIdx + 1
      
        ! get local FEElementLocalNumber from FEElementGlobalNumber
        CALL cmfe_BioelectricFiniteElasticity_GetLocalElementNumber(GeometricFieldFE, FEElementGlobalNumber, FEElementLocalNumber, &
          & Err)
      
        !PRINT "(I1.1,2(A,I2.2),6(A,I3.3))", ComputationalNodeNumber,": Fibre", FibreNo, ", Element in fibre ", ElementIdx, &
        !  & ", Element (",FEElementXIdx,",",FEElementYIdx,",",FEElementZIdx, &
        !  & ")=", FEElementGlobalNumber, ", MU ",MotorUnitRank, ", local FE ", FEElementLocalNumber
      
        ! set number of containing FE element
        !                                     FIELD,              VARIABLE_TYPE,              
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
        !   FIELD_SET_TYPE,             VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE, ERR
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  5,                FEElementLocalNumber, Err)
          
      
        ! set motor unit number
        !                                      FIELD,             VARIABLE_TYPE               
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
          ! FIELD_SET_TYPE              VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  1,                MotorUnitRank, Err)
        
        ! set fibre number
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM, CMFE_FIELD_V_VARIABLE_TYPE, &
          ! FIELD_SET_TYPE              VERSION_NO, DERIVATIVE_NO, USER_NODE_NUMBER, COMPONENT_NUMBER, VALUE
          & CMFE_FIELD_VALUES_SET_TYPE, 1,          1,             NodeMUserNumber,  3,                FibreNo, Err)
        
      ENDIF
    ENDDO
  ENDDO
  
END SUBROUTINE InitializeFieldMonodomain

SUBROUTINE InitializeFieldFiniteElasticity()
  INTEGER(CMISSIntg) :: ElementUserNumber
  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  !second variable of IndependentFieldFE
  !  components:
  !    1) number of nodes in Xi(1) direction per element
  !    2) number of nodes in Xi(2) direction per element
  !    3) number of nodes in Xi(3) direction per element
  !    4) beginning of fibres in this FE element? 1=yes, 0=no
  !    5) number of in series fibres
  !
  !initialise as if the fibres would not start in any element, and adjust below
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & NumberOfNodesInXi1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & NumberOfNodesInXi2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & NumberOfNodesInXi3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & 0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
   & NumberOfInSeriesFibres,Err)

  ! Set V component 4) to 1 if fibre starts in element
  DO ElementUserNumber=1, NumberOfElementsFE, NumberGlobalXElements/NumberOfInSeriesFibres

    ! only if element is assigned to own domain
    !                                        DECOMPOSITION,   USER_ELEMENT_NUMBER, DOMAIN_NUMBER
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE, ElementUserNumber,   ElementDomain,Err)
    IF (ElementDomain == ComputationalNodeNumber) THEN
      !fibres begin in this element
      !                                        FIELD,              VARIABLE_TYPE,             FIELD_SET_TYPE,
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      !  USER_ELEMENT_NUMBER, COMPONENT_NUMBER, VALUE
       & ElementUserNumber,   4,                1,     Err)
      !PRINT "(A,I3.3,A,I5.5,A)", "x ", ComputationalNodeNumber, ": 3D el. ", ElementUserNumber, " has beginning fibre "
    ENDIF
  ENDDO

END SUBROUTINE InitializeFieldFiniteElasticity

SUBROUTINE CreateEquations()
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  
  CALL cmfe_Equations_Initialise(EquationsM,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsM,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsM,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(EquationsM,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetM,Err)
    
  !Create the equations set equations for Finite Elasticity
  CALL cmfe_Equations_Initialise(EquationsFE,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  !CALL cmfe_Equations_SparsityTypeSet(EquationsFE,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsFE,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(EquationsFE,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  IF (DebuggingOutput) THEN
    CALL cmfe_Equations_OutputTypeSet(EquationsFE,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  ENDIF
  
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetFE,Err)

END SUBROUTINE CreateEquations

SUBROUTINE InitializeCellML()
  REAL(CMISSRP) :: InitialPressure
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,RegionM,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
  
  !                            CellML env.    ,URI,                modelIndex (out)
  CALL cmfe_CellML_ModelImport(CellML,CellMLModelFilename,shortenModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !,- to set from this side

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"razumova/L_S",Err) !(= l_{hs})
   ! muss theoretisch auch übergeben werden, ist aber noch nicht im subcell model enthalten:
   ! CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"razumova/velo",Err) ( = d{l_{hs}} / d{t})
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Aliev_Panfilov/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/l_hs",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/velo",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Aliev_Panfilov/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/l_hs",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/rel_velo",Err)
  ENDIF
!  CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"razumova/rel_velo",Err)
!
!  CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex2,"wal_environment/I_HH",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex2,"razumova/L_S",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex2,"razumova/rel_velo",Err)
  !,- to get from the CellML side
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !,  --> set "razumova/stress" as wanted!
  !,  --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted!
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
  ! wird in CellML als Stress berechnet, dann aber als Rate ausgegeben. passt das so?: Achtung, Modell nicht implementiert in ModelType==0.
    CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"razumova/stress",Err)
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"Razumova/ActiveStress",Err)
    CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"Razumova/Activation",Err)
  ENDIF
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex2,"razumova/stress",Err)
  !,- and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  
  IF (ODESolverId==2) THEN      ! BDF solver
   ! To initialize the BDF solver, OpenCMISS requires to have set CELLML%MAX_NUMBER_OF_INTERMEDIATE.
   ! For now, we can manipulate it here:
    CALL cmfe_CellML_IntermediateMaxNumberSet(CellML,1,Err) ! todo: exact number required or just something > 0?
  ENDIF
  
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  
  !Map the half-sarcomere length L_S
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !Map the transmembrane voltage V_m
    !                                       CellML env., "from"-field
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,      DependentFieldM, & 
    !  variable type,          comp. no., "from"-field param. set
     & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    !  "to"-cellml model user no., "to"-variable,      "to"-cellml variable parameter set
     & shortenModelIndex,          "wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
     
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS", &
     & CMFE_FIELD_VALUES_SET_TYPE,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/stress", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM, &
     & CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/l_hs",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the sarcomere relative contraction velocity
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM, & 
     & CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM, &
     & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Aliev_Panfilov/V_m", &
     & CMFE_FIELD_VALUES_SET_TYPE,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/A_1", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/A_2", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/x_1", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/x_2", &
     & CMFE_FIELD_VALUES_SET_TYPE,IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    !Map the sarcomere half length L_S
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/l_hs",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the sarcomere relative contraction velocity
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/rel_velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
      & shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE, &
     & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/ActiveStress",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/Activation",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_VALUES_SET_TYPE,Err)

  ENDIF

!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/L_S",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the sarcomere relative contraction velocity
!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/rel_velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the transmembrane voltage V_m
!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
!   & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
!  !Map the active stress
!  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"razumova/stress",CMFE_FIELD_VALUES_SET_TYPE, &
!   & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & -79.974_CMISSRP,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & 0.0_CMISSRP,Err)

  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & 0.0_CMISSRP,Err)

  ENDIF

  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,1,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,2,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,3,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)

  IF (ModelType == 0 .OR. ModelType == 1) THEN    ! 3,3a, "MultiPhysStrain"
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
     & 0.0_CMISSRP,Err)

  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    InitialPressure = -2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
     & InitialPressure,Err)

  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !Set up the models field
  CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & shortenModelIndex,Err)

  !Create the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field
  IF (ModelType==0 .OR. ModelType==2 .OR. ODESolverId==2) THEN ! ModelType==0: 3a, "MultiPhysStrain", ModelType==2: 4, "Titin", ODESolverId== 2: BDF
    CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber, & 
     & CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  ENDIF

  !Create the CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

END SUBROUTINE InitializeCellML

SUBROUTINE CreateProblem()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE],Problem,Err)
    ENDIF
  ELSEIF (ModelType == 1) THEN ! 3, , "MultiPhysStrain", numerically more stable
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE],Problem,Err)
      ENDIF
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
       & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
       & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE],Problem,Err)
    ENDIF
  ENDIF

  CALL cmfe_Problem_CreateFinish(Problem,Err)

END SUBROUTINE CreateProblem

SUBROUTINE CreateControlLoops()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for StimDuration with the Stimulus applied.
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,ElasticityTimeStep,ElasticityTimeStep,Err)

  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OutputTimeStepStride,Err)
  IF (DebuggingOutput) THEN
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_TIMING_OUTPUT,Err)
  ELSE
    ! output types:
    ! CONTROL_LOOP_NO_OUTPUT = 0 !<No output from the control loop (no output of MainTime_* files)
    ! CONTROL_LOOP_PROGRESS_OUTPUT = 1 !<Progress output from control loop (also output MainTime_* files)
    ! CONTROL_LOOP_TIMING_OUTPUT = 2 !<Timing output from the control loop (also output MainTime_* files)
    ! CONTROL_LOOP_FILE_OUTPUT = -1 <Only MainTime_* files output
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_TIMING_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,cmfe_CONTROL_LOOP_FILE_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,-1,Err)
  ENDIF


  !set the monodomain loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoopM,0.0_CMISSRP,ElasticityTimeStep,PDETimeStep,Err)
  ! question from Aaron: is the following necessary, or is the TimesSet-CALL before sufficient?
  IF (PdeNSteps /= -1) THEN
    CALL cmfe_ControlLoop_NumberOfIterationsSet(ControlLoopM,PdeNSteps,Err)
  ENDIF
  
  IF (DebuggingOutput) THEN
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
  ELSE
    CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
    !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_FILE_OUTPUT,Err)
  ENDIF

  !set the finite elasticity loop (simple type)
  IF (.NOT. elasticityDisabled) THEN
    CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
    CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
    CALL cmfe_ControlLoop_TypeSet(ControlLoopFE,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
    CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,ElasticityLoopMaximumNumberOfIterations,Err)
    CALL cmfe_ControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

    IF (DebuggingOutput) THEN
      CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopFE,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
    ELSE
      CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopFE,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)
    ENDIF
  ENDIF
    
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

END SUBROUTINE CreateControlLoops

SUBROUTINE CreateSolvers()
  TYPE(cmfe_SolverType) :: linearSolver

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  ! create parabolic and ode solver, create subsolver structure
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_DAETimeStepSet(SolverDAE,ODETimeStep,Err)  ! #timestepset
  
  ! might be handy some day:
  ! CALL cmfe_Solver_DAETimesSet(SolverDAE,?0.0_CMISSRP?,?0.001_CMISSRP?,Err)
  
  SELECT CASE(ODESolverId) !use bdf instead of default-explicitEuler
    CASE(2) ! BDF
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_BDF,Err)
      CALL cmfe_Solver_DAEbdfSetTolerance(SolverDAE,0.0000001_CMISSRP,0.0000001_CMISSRP,err) !default values were both: 1.0E-7
    CASE(3) ! GL, not stable yet
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_GL,Err)
    CASE(4) ! Crank-Nicolson, not stable yet
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_CRANK_NICOLSON,Err)
    CASE(5) ! improved Euler (Heun)
      CALL cmfe_Solver_DAEEulerSolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EULER_IMPROVED,Err)
    CASE DEFAULT
      IF (ComputationalNodeNumber == 0) THEN
        PRINT *, ""
        PRINT *, "Warning: For the DAE Problem (0D) the standard explicit Euler method is used, which is slow. " // NEW_LINE('A') &
          & // "          Consider to use another DAE solver via the command line option 'ODESolverId'."
        PRINT *, ""
      ENDIF
  END SELECT
  !-------------------------------------------------------------------------------------------  
  !Set the Number of ODE time steps. CARE: This makes cmfe_Solver_DAETimeStepSet() obsolete!
  IF(ODESolverId==1 .AND. OdeNSteps/=-1) THEN
    CALL cmfe_Solver_DAEEulerForwardSetNSteps(SolverDAE,OdeNSteps,Err)
  ELSEIF(ODESolverId==5 .AND. OdeNSteps/=-1) THEN
    CALL cmfe_Solver_DAEEulerImprovedSetNSteps(SolverDAE,OdeNSteps,Err)
  END IF
  
  !> \todo or not-todo... solve the CellML equations on the GPU for efficiency (later)
  !CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EXTERNAL,Err)

  IF (DebuggingOutput) THEN
    CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_TIMING_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  !-----------------------------------------------------------------------------------------------------
  !Create the parabolic solver
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  

  IF (SplittingType == 1) THEN    ! Strang splitting      
    CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,Err)
  ELSE    ! Godunov splitting
    CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  ENDIF
  
  ! data structure is as follows:
  ! SolverParabolic
  !   LinearSolver
  !      IterativeSolver
  !   LinkedSolver(1) -> LinearSolver
  ! Retrieve linear solver
  NULLIFY(linearSolver%solver)
  CALL cmfe_Solver_DynamicLinearSolverGet(SolverParabolic, linearSolver, Err)

  ! direct:
  ! CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  ! CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
  ! CMFE_SOLVER_DIRECT_LU
  ! CMFE_SOLVER_DIRECT_CHOLESKY
  ! CMFE_SOLVER_DIRECT_SVD
  !
  ! iterative:
  ! CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE, Err)
  ! CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
  ! CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER, Err)
  ! solver:
  ! CMFE_SOLVER_ITERATIVE_GMRES
  ! CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT
  ! CMFE_SOLVER_ITERATIVE_CONJGRAD_SQUARED
  ! preconditioner:
  ! CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_JACOBI_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_SOR_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER
  
  !MonodomainSolverId, MonodomainPreconditionerId, ODESolverId
  IF (MonodomainSolverId <= 1) THEN    ! direct solver
    
    CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  
    SELECT CASE(MonodomainSolverId)
      CASE(1)
        CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_DIRECT_LU, Err)
      CASE DEFAULT
        PRINT*, "Warning: Wrong MonodomainSolverId ",MonodomainSolverId
        CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_DIRECT_LU, Err)
    END SELECT
  ELSE ! iterative solver
    
    CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE, Err)
    
    SELECT CASE(MonodomainSolverId)
      CASE(2)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_GMRES, Err)
      CASE(3)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
      CASE(4)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJGRAD_SQUARED, Err)
      CASE DEFAULT
        PRINT*, "Warning: Wrong MonodomainSolverId ",MonodomainSolverId
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_GMRES, Err)
    END SELECT
    
    SELECT CASE(MonodomainPreconditionerId)
      CASE(1)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER, Err)
      CASE(2)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_JACOBI_PRECONDITIONER, Err)
      CASE(3)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, Err)
      CASE(4)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_SOR_PRECONDITIONER, Err)
      CASE(5)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, &
         & CMFE_SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, Err)
      CASE(6)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER, Err)
      CASE(7)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, & 
         & CMFE_SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER, Err)
      CASE DEFAULT
        PRINT *, "Warning: Wrong MonodomainPreconditionerId, ", MonodomainPreconditionerId
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER, Err)
    END SELECT
    
  ENDIF
  
  IF (DEBUGGING_PROBLEM_OUTPUT) THEN
    PRINT*, ""
    PRINT*, "before cmfe_Solver_LinearTypeSet"
    CALL cmfe_PrintSolver(SolverParabolic, 5, 10, Err)
  ENDIF
  
  ! Recreate linearSolver subtype as direct solver (instead of iterative solver)
  !CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  
  IF (DEBUGGING_PROBLEM_OUTPUT) THEN
    PRINT*, ""
    PRINT*, "========================================================================="
    PRINT*, "After cmfe_Solver_LinearTypeSet"
    CALL cmfe_PrintSolver(SolverParabolic, 5, 10, Err)
  ENDIF
  
  !SOLVER_LINEAR_DIRECT_TYPE_SET
  
  !CALL cmfe_SOLVER_DYNAMICLINEARITYTYPESET(SolverParabolic, CMFE_SOLVER_DYNAMIC_LINEAR, Err)
  !CALL cmfe_Solver_NonlinearTypeSet(SolverParabolic,CMFE_SOLVER_NONLINEAR_NEWTON,Err)

  IF (DebuggingOutput) THEN
    CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_TIMING_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  !Create the Finte Elasticity solver
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_Solver_Initialise(LinearSolverFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverFEIndex,SolverFE,Err)

  ! only enable output for specified ComputeNode
  IF (DebuggingOutput) THEN
    IF (ComputationalNodeNumber == 0) THEN
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_TIMING_OUTPUT,Err)
      CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
    ELSE
      CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
    ENDIF
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
!  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err) ! 721 steps
  CALL cmfe_Solver_NewtonMaximumIterationsSet(SolverFE,NewtonMaximumNumberOfIterations,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)

  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)
  
  

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  
  !                           problem, control loop indetifiers (in),                    solverindex (in), solver (out)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  
  !                                   in              in     out
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
    
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsM,Err)
  
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
   
   
  CALL cmfe_Solver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  
  
  !PRINT*, ""
  !PRINT*, "After cmfe_Solver_SolverEquationsGet"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  
  
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_SPARSE_MATRICES,Err)
  
  
  !PRINT*, ""
  !PRINT*, "After cmfe_SolverEquations_SparsityTypeSet"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_FULL_MATRICES,Err)
  
  
  ! TYPE(cmfe_SolverEquationsType), INTENT(IN) :: solverEquations !<The solver equations to add the equations set for.
  ! TYPE(cmfe_EquationsSetType), INTENT(IN) :: equationsSet !<The equations set to add.
  ! INTEGER(INTG), INTENT(OUT) :: equationsSetIndex !<On return, the index of the added equations set in the solver equations
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)
  
  !PRINT*, ""
  !PRINT*, "cmfe_PrintSolverEquationsType:"
  !CALL cmfe_PrintSolverEquations(SolverEquationsM, 5, 10, Err)
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, Err)
  !STOP

  !Create the problem solver Finite Elasticity equations
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],SolverFEIndex,SolverFE,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_FULL_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  IF (DEBUGGING_PROBLEM_OUTPUT) THEN
    PRINT*, "==========================="
    PRINT*, "SolverFE"
    CALL cmfe_PrintSolver(SolverFE, 6, 3, Err)
  ENDIF
  
  !SolverFE%solver%NONLINEAR_SOLVER%NEWTON_SOLVER%LINEAR_SOLVER%OUTPUT_TYPE = SOLVER_MATRIX_OUTPUT
  
  !STOP
  
END SUBROUTINE CreateSolvers

SUBROUTINE SetBoundaryConditions()
  INTEGER(CMISSIntg) :: NodeMUserNumber
  INTEGER(CMISSIntg) :: NodeIdx

  !Prescribe boundary conditions for monodomain
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  
  !PRINT*, "Print SolverEquationsM Solver"
  
  !CALL cmfe_PrintSolver(SolverEquationsM, 5, 10, Err)
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  
  !SolverEquationsM
  !SolverEquationsM%solverEquations%SOLVER%SOLVERS%SOLVERS(2)% &
  !  & PTR%LINKED_SOLVERS(1)%PTR%LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE = 1e-5 


  !PRINT*, "After changing tolerance"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  ! Finalize solver, set PETSC parameters, such as tolerances, solver type etc.
  ! The parameters are already set in the data structure SolverEquations
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x.
  DO NodeIdx = 1, SIZE(LeftSurfaceNodes, 1)
    NodeNumber = LeftSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE, NodeNumber, 1, NodeDomain, Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      !                                    BOUNDARY_CONDITIONS,  FIELD,            VARIABLE_TYPE,
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE, DependentFieldFE, CMFE_FIELD_U_VARIABLE_TYPE, &
      !   VERSION_NUMBER, DERIVATIVE_NUMBER, USER_NODE_NUMBER, COMPONENT_NUMBER, CONDITION,
        & 1,              1,                 NodeNumber,       1,                CMFE_BOUNDARY_CONDITION_FIXED, &   ! CMFE_BOUNDARY_CONDITION_FIXED_USER_CONTROLLED
      ! VALUE
        & 0.0_CMISSRP, Err)
    ENDIF
  ENDDO

 !Set x=PhysicalWidth nodes to InitialStretch x displacement
  DO NodeIdx = 1, SIZE(RightSurfaceNodes, 1)
    NodeNumber = RightSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE, NodeNumber, 1, NodeDomain, Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, PhysicalLength*InitialStretch, Err)
        
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber,&
      !  & 1,CMFE_BOUNDARY_CONDITION_FIXED, 0.0_CMISSRP, Err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO NodeIdx = 1, SIZE(FrontSurfaceNodes, 1)
    NodeNumber = FrontSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO NodeIdx = 1, SIZE(BottomSurfaceNodes, 1)
    NodeNumber = BottomSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

END SUBROUTINE SetBoundaryConditions

SUBROUTINE CalculateBioelectrics()
  !Calculate the bioelectrics geometric field
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)

  !PRINT*, "before cmfe_BioelectricsFiniteElasticity_UpdateGeometricField"
  !CALL cmfe_OutputInterpolationParameters(Problem,DependentFieldM, SolverParabolic,Err)
  
  CALL cmfe_BioelectricsFiniteElasticity_UpdateGeometricField(ControlLoopM,.TRUE.,Err)

  CALL SLEEP(2)
  
  !PRINT*, "after cmfe_BioelectricsFiniteElasticity_UpdateGeometricField"
  !CALL cmfe_OutputInterpolationParameters(Problem,DependentFieldM, SolverParabolic,Err)
  
  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)
END SUBROUTINE CalculateBioelectrics

SUBROUTINE ExportEMG()
  REAL(CMISSSP), DIMENSION(2) :: TimeStart, TimeStop
  REAL(CMISSSP) :: Total
 
  IF (.NOT. EnableExportEMG) RETURN
   
  IF (ComputationalNodeNumber == 0) WRITE(*,'(A)',advance='no') "Output EMG Data ... "
  CALL ETIME(TimeStart, Total)
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EMGExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EMGExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EMGExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EMGExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  CALL ETIME(TimeStop, Total)
  TimingExportEMGUser = TimingExportEMGUser + (TimeStop(1) - TimeStart(1))
  TimingExportEMGSystem = TimingExportEMGSystem + (TimeStop(2) - TimeStart(2))
  
  IF (ComputationalNodeNumber == 0) PRINT*, " done"

END SUBROUTINE ExportEMG

SUBROUTINE SetStimulationAtNodes(StimValuePerNode)
  REAL(CMISSRP), INTENT(IN) :: StimValuePerNode
  INTEGER(CMISSIntg) :: StimulatedNodeBegin, StimulatedNodeEnd, StimulatedNodeNo
  LOGICAL :: CurrentFibreFires
  
  NumberFiringFibres = 0
  !loop over all neuromuscular junctions (middle point of the fibres)
  DO FibreNo = 1, NumberOfFibres
    CurrentFibreFires = .FALSE.
    
    ! get middle point of fibre
    JunctionNodeNo = (FibreNo-1) * NumberOfNodesPerLongFibre + (NumberOfNodesPerLongFibre+1)/2
    
    ! add innervation zone offset
    JunctionNodeNo = JunctionNodeNo + InnervationZoneOffset(FibreNo)
    
    ! compute first node to stimulation
    StimulatedNodeBegin = JunctionNodeNo - NumberStimulatedNodesPerFibre/2
    
    ! assert limits
    StimulatedNodeBegin = MIN(FibreNo*NumberOfNodesPerLongFibre-NumberStimulatedNodesPerFibre, &
      & MAX((FibreNo-1)*NumberOfNodesPerLongFibre, StimulatedNodeBegin))
    
    StimulatedNodeEnd = StimulatedNodeBegin + NumberStimulatedNodesPerFibre-1
    
    ! loop over nodes of current fibre to be stimulated
    DO StimulatedNodeNo = StimulatedNodeBegin, StimulatedNodeEnd
      
      !                                     decomposition,  nodeUserNumber, meshComponentNumber, domain
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM, StimulatedNodeNo, 1,                   NodeDomain, Err)
      IF (NodeDomain == ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & StimulatedNodeNo,1,MotorUnitRank,Err)

        IF ((MotorUnitRank <= 0) .OR. (MotorUnitRank >= 101)) THEN
          PRINT*, ComputationalNodeNumber,": Warning! MotorUnitRank=",MotorUnitRank,", set to 100. Fibre ", FibreNo, &
             & ", Node", StimulatedNodeNo
          MotorUnitRank=100
        ELSE
          !PRINT*, "MotorUnitFiringTimes row k=", k, ": MU rank=", MotorUnitRank, ", StimComponent=",StimComponent
          MotorUnitFires = MotorUnitFiringTimes(k, MotorUnitRank)   ! determine if mu fires
          IF (MotorUnitFires == 1) THEN
            !PRINT*, "Fibre ",FibreNo,": MU ", MotorUnitRank, " fires, StimulatedNodeNo=",StimulatedNodeNo, &
            !  & ", StimComponent=",StimComponent,", StimValue=", StimValue

            CurrentFibreFires = .TRUE.
            CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,StimulatedNodeNo,StimComponent,StimValuePerNode,Err)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
    ! If there was a stimulated node in the current fibre, increase number of stimulated fibres counter for output
    IF (CurrentFibreFires) NumberFiringFibres = NumberFiringFibres + 1
  ENDDO

END SUBROUTINE SetStimulationAtNodes
  
FUNCTION GetTimeStamp()
  INTEGER(CMISSIntg), DIMENSION(3) :: Today, Now
  CHARACTER(LEN=100) :: TimeString
  CHARACTER(LEN=100) :: GetTimeStamp

  CALL IDATE(Today)  ! day, month, year
  CALL ITIME(Now)    ! h, min, s

  WRITE (TimeString, '(I2.2, A, I2.2, A, I4.4, A, I2.2, A, I2.2, A, I2.2)') &
  & Today(1), '.', Today(2), '.', Today(3), ' ', Now(1), ':', Now(2), ':', Now(3)

  GetTimeStamp = TimeString
END FUNCTION GetTimeStamp

FUNCTION GetMemoryConsumption()
  CHARACTER(LEN=100) :: GetMemoryConsumption
  CHARACTER(LEN=100) ::  &
   & pid, comm, state, ppid, pgrp, session, tty_nr, &
   & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
   & utime, stime, cutime, cstime, priority, nice, &
   & O, itrealvalue, starttime, Description, Limit
  CHARACTER(LEN=10000) :: Debug
  INTEGER(CMISSIntg) :: I
  INTEGER(CMISSLintg) :: MemoryConsumption, &
   & VmSize, VmRSS, Shared, Text, Lib, Data, Dt, RssLimBytes, RssAnon, Pagesize, VSizeBytes

  ! Critical Section
  DO I=0,NumberOfComputationalNodes
    !CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
    IF (I == ComputationalNodeNumber) THEN

      ! read memory page size
      Stat = 0
      IF (I==0) CALL SYSTEM("getconf PAGESIZE > pagesize", Stat)
      IF (Stat == 0) THEN
        OPEN(UNIT=10, FILE="pagesize", ACTION="read", IOSTAT=Stat)
        IF (Stat == 0) THEN
          READ(10,*, IOSTAT=Stat) Pagesize
          CLOSE(UNIT=10)
        ELSE
          PRINT*, "Error opening pagesize."
        ENDIF
      ELSE
        PRINT*, "Error calling 'getconf PAGESIZE'"
        Pagesize = 4096
      ENDIF

      ! read from /proc/self/stat
      ! see http://man7.org/linux/man-pages/man5/proc.5.html for reference
      OPEN(UNIT=10, FILE="/proc/self/stat", ACTION="read", IOSTAT=stat)
      IF (STAT /= 0) THEN
        PRINT*, "Could not read memory consumption from /proc/self/stat."
      ELSE
        !READ(10,"(A)",IOSTAT=stat,advance='no') Debug
        !PRINT*, "proc/self/stat: "//TRIM(Debug)
        READ(10,*, IOSTAT=stat) pid, comm, state, ppid, pgrp, session, tty_nr, &
           & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
           & utime, stime, cutime, cstime, priority, nice, &
           & O, itrealvalue, starttime, VSizeBytes, VmRSS, RssLimBytes
        CLOSE(UNIT=10)

        MemoryConsumption = VSizeBytes

        IF (ComputationalNodeNumber == 0) THEN
          WRITE(*, "(A,F7.3,A,I11,A,F7.3,A)") "     MemoryConsumption:", (MemoryConsumption/1e9), " GB (", &
            & MemoryConsumption, " B), Resident: ", (VmRSS*PageSize/1e9), " GB"
        ENDIF

        WRITE(GetMemoryConsumption, *) MemoryConsumption
      ENDIF

      ! read from /proc/self/limits
      IF (.TRUE.) THEN
        OPEN(UNIT=10, FILE="/proc/self/limits", ACTION="read", IOSTAT=stat)
        IF (STAT /= 0) THEN
          PRINT*, "Could not read limits from /proc/self/limits."
        ELSE
          DO
            READ(10, "(A26,A21)", IOSTAT=Stat) Description, Limit
            !PRINT*, "Description:["//TRIM(Description)//"], Limit:["//TRIM(Limit)//"]"
            IF (Stat /= 0) EXIT
            !IF (INDEX(Description, "Max resident set") /= 0) THEN
              !IF (TRIM(Limit) == "unlimited") THEN
              !  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident has no soft limit)"
              !ELSE
              !  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident is limited to ", Limit,")"
              !ENDIF
            !ENDIF
          ENDDO
          CLOSE(UNIT=10)
        ENDIF
      ENDIF

      ! read from /proc/self/statm
      OPEN(UNIT=10, FILE="/proc/self/statm", ACTION="read", IOSTAT=stat)
      IF (STAT /= 0) THEN
        PRINT*, "Could not read memory consumption from /proc/self/statm."
      ELSE
        READ(10,*, IOSTAT=stat) VmSize, VmRSS, Shared, Text, Lib, Data, Dt
        CLOSE(UNIT=10)

        ! VmRSS = RssAnon + Shared (all as number of pages)
        ! RssAnon is the resident set size of anonymous memory (real memory in RAM, not laid out in files)
        RssAnon = VmRSS - Shared

        IF (ComputationalNodeNumber == 0) THEN
          !PRINT*, "VmSize: ", VmSize, ", VmRSS: ", VmRSS, ", Shared: ", Shared, ", Text: ", Text, ", Lib:", Lib, ", Data:", Data
          !PRINT*, "RssAnon: ", RSSAnon, ", RssLimBytes: ", RssLimBytes, ", Pagesize: ", Pagesize

          ! Output Percentage
          !WRITE(*, "(3(A,F7.3),A,F5.1,A)") "     VmSize:", (VmSize*Pagesize/1.e9), " GB, RssAnon:", (RssAnon*Pagesize/1.e9), &
          !  & " GB, RssLimBytes: ", (RssLimBytes/1.e9), " GB (", (REAL(RssAnon*Pagesize) / RssLimBytes * 100.0), "%)"
        ENDIF
      ENDIF

      !PRINT*, TRIM(cmfe_CustomProfilingGetInfo(Err))
    ENDIF
  ENDDO

END FUNCTION GetMemoryConsumption

SUBROUTINE WriteTimingFile()

  CHARACTER(len=256) :: Filename = "duration."
  CHARACTER(len=20)  :: ComputationalNodeNumberStr
  CHARACTER(len=100) :: Hostname
  LOGICAL :: FileExists
  REAL(CMISSSP), DIMENSION(2) :: Elapsed     ! For receiving user and system time
  REAL(CMISSSP) :: DurationTotal
  REAL(CMISSRP) :: DurationInit, DurationStretchSim, DurationIntInit, DurationMainSim
  CHARACTER(LEN=100) :: TimeStampStr, MemoryConsumptionEnd

  ! create filename
  WRITE(ComputationalNodeNumberStr, '(I0.5)') ComputationalNodeNumber     ! convert ComputationalNodeNumber to string
  Filename = TRIM(Filename) // TRIM(ComputationalNodeNumberStr) // ".csv"

  ! Check if file exists
  INQUIRE(file=Filename, exist=FileExists)

  ! Write Comment in first line if file does not yet exist
  IF(.NOT. FileExists) THEN
    OPEN(unit=123, file=Filename, iostat=stat)
    IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'
    WRITE(123,'(A)') '# Stamp; Host; NProc; X; Y; Z; F; Total FE; Total M; End Time; ' // &
      & 'Dur. Init; Stretch Sim; Int. Init; Main Sim; Total; Total (User); Total (System); ' // &
      & 'ODE; Parabolic; FE; FE before Main Sim; Mem. Consumption after 1st timestep; Memory Consumption At End; ' // &
      & 'Parabolic reason; Newton reason; parabolic n. iter; min; max; newton n. iter; min; max; ' // &
      & 'level 0: stimulation handling; level 0: problem solve; level 1: MAIN_TIME_LOOP overhead; ' // &
      & 'level 1: MONODOMAIN_TIME_LOOP overhead; level 1: ELASTICITY_LOOP overhead; level 1: SolverDAE solve; ' // &
      & 'level 1: SolverParabolic solve; level 1: SolverFE solve; level 1: interpolate 1D->3D; level 1: interpolate 3D->1D; ' // &
      & 'level 1: file output; level 2: solver overhead; level 2: 0D solve; level 2: 1D solve; level 2: 3D solve; ' // &
      & 'level 3: 1D assembly; level 3: 1D solve; level 3: 1D other; level 3: 3D assembly; level 3: 3D solve; ' // &
      & 'level 3: 3D other;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;' // &
      & 'duration FESolverPreLoad; duration OdeSolverPreLoad; duration ParabolicSolverPreLoad; ' // &
      & 'duration FileOutputPreLoad (user); duration export EMG user; duration export EMG system; duration FileOutput (user); ' // &
      & 'duration FileOutput (system); duration FileOutputPreload (system); MonodomainSolverId; MonodomainPreconditionerId; ' // &
      & 'ODESolverId; NumberOfElementsInAtomX; NumberOfElementsInAtomY; NumberOfElementsInAtomZ; NumberOfSubdomainsX; ' // & 
      & 'NumberOfSubdomainsY; NumberOfSubdomainsZ; ModelType; ScenarioName'
      
    CLOSE(unit=123)
  ENDIF

  ! Write line in file
  OPEN(unit=123, file=Filename, iostat=stat, access='append')
  IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

  DurationInit = TimeInitFinshed - TimeStart
  DurationStretchSim = TimeStretchSimFinished - TimeInitFinshed
  DurationIntInit = TimeMainSimulationStart - TimeStretchSimFinished
  DurationMainSim = TimeMainSimulationFinished - TimeMainSimulationStart

  CALL ETIME(Elapsed, DurationTotal)
  CALL HOSTNM(Hostname)
  MemoryConsumptionEnd = GetMemoryConsumption()

  TimeStampStr = GetTimeStamp()

  IF (CustomProfilingEnabled) THEN

    WRITE(123,"(4A,7(I11,A),(F8.3,A),11(F0.8,A),2(A,A),8(I7,A),21(F25.13,A),A,A,9(F8.3,A),3(I2,A),6(I8,A),I1,2A)") &
      & TRIM(TimeStampStr), ';', &
      & TRIM(Hostname(1:22)), ';', &          ! end of 4A
      & NumberOfComputationalNodes, ';', &
      & NumberGlobalXElements, ';', &
      & NumberGlobalYElements, ';', &
      & NumberGlobalZElements, ';', &
      & NumberOfInSeriesFibres, ';', &
      & NumberOfElementsFE, ';', &
      & NumberOfElementsM, ';', &             ! end of 7(I11,A)
      & TimeStop, ';', &                      ! (F8.3,A)
      & DurationInit, ';', &
      & DurationStretchSim, ';', &
      & DurationIntInit, ';', &
      & DurationMainSim,';',  &
      & DurationTotal, ';', &
      & Elapsed(1), ';', &
      & Elapsed(2), ';', &
      & CustomTimingOdeSolver, ';', &
      & CustomTimingParabolicSolver, ';', &
      & CustomTimingFESolver, ';', &
      & CustomTimingFESolverPreLoad, ';', &   ! end of 11(F0.8,A)
      & TRIM(ADJUSTL(MemoryConsumption1StTimeStep)), ';', &
      & TRIM(ADJUSTL(MemoryConsumptionEnd)), ';', &   ! end of 2(A,A)
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';', &   ! end of 8(I7,A)
      ! Custom Profiling durations
      & cmfe_CustomProfilingGetDuration("level 0: stimulation handling", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 0: problem solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: MAIN_TIME_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: MONODOMAIN_TIME_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: ELASTICITY_LOOP overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverDAE solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverParabolic solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: SolverFE solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: interpolate 1D->3D", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: interpolate 3D->1D", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 1: file output", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: solver overhead", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 0D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 1D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 2: 3D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D assembly", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 1D other", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D assembly", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D solve", Err), ';', &
      & cmfe_CustomProfilingGetDuration("level 3: 3D other", Err), ';', &                 ! end of 21(F25.13,A)
      & ';;;;;;;;;;;;;;', &             ! end of A (14 free columns)
      & ';;;;;;;;;;;;;;;;;;;;;;;;', &   ! end of A (3*8=24 free columns) (8(I17,A,I5,A,I7,A))
      ! custom profiling memory consumption
      !& cmfe_CustomProfilingGetMemory("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector cmiss DP", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector cmiss INTG", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage diag", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, offdiag", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, row ind.", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage, col. ind.", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed matrix petsc, compr. row storage (local to global mapping)", Err), ';', &
      !& cmfe_CustomProfilingGetMemory("distributed vector petsc", Err), ';', &
      !& cmfe_CustomProfilingGetSizePerElement("distributed vector petsc", Err), ';', &
      !& cmfe_CustomProfilingGetNumberObjects("distributed vector petsc", Err), ';', &       ! end of 8(I17,A,I5,A,I7,A)
      & CustomTimingFESolverPreLoad, ';', &
      & CustomTimingOdeSolverPreLoad, ';', &
      & CustomTimingParabolicSolverPreLoad, ';', &
      & CustomTimingFileOutputUserPreLoad, ';', &
      & TimingExportEMGUser, ';', &
      & TimingExportEMGSystem, ';', &
      & CustomTimingFileOutputUser, ';', &
      & CustomTimingFileOutputSystem, ';', &
      & CustomTimingFileOutputSystemPreLoad, ';', &     ! end of 9(F8.3,A)
      & MonodomainSolverId, ';', &
      & MonodomainPreconditionerId, ';', &
      & ODESolverId, ';', &                             ! end of 3(I2,A)
      & NumberOfElementsInAtomX, ';', &
      & NumberOfElementsInAtomY, ';', &
      & NumberOfElementsInAtomZ, ';', &                 
      & NumberOfSubdomainsX, ';', &
      & NumberOfSubdomainsY, ';', &
      & NumberOfSubdomainsZ, ';', &                            ! end of 6(I8,A)
      & ModelType, ';', &                                       ! I1
      & TRIM(ScenarioName)

  ELSE  ! custom profiling is disabled
    
    WRITE(123,"(4A,7(I11,A),(F8.3,A),11(F0.8,A),2(A,A),8(I7,A),9(F8.3,A),3(I2,A),6(I8,A),I1,2A)") &
      & TRIM(TimeStampStr), ';', &
      & TRIM(Hostname(1:22)), ';', &                    ! end of 4A
      & NumberOfComputationalNodes, ';', &
      & NumberGlobalXElements, ';', &
      & NumberGlobalYElements, ';', &
      & NumberGlobalZElements, ';', &
      & NumberOfInSeriesFibres, ';', &
      & NumberOfElementsFE, ';', &
      & NumberOfElementsM, ';', &                       ! end of 7(I11,A)
      & TimeStop, ';', &                                ! (F8.3,A)
      & DurationInit, ';', &
      & DurationStretchSim, ';', &
      & DurationIntInit, ';', &
      & DurationMainSim,';',  &
      & DurationTotal, ';', &
      & Elapsed(1), ';', &
      & Elapsed(2), ';', &
      & CustomTimingOdeSolver, ';', &
      & CustomTimingParabolicSolver, ';', &
      & CustomTimingFESolver, ';', &
      & CustomTimingFESolverPreLoad, ';', &             ! end of 11(F0.8,A)
      & TRIM(ADJUSTL(MemoryConsumption1StTimeStep)), ';', &
      & TRIM(ADJUSTL(MemoryConsumptionEnd)), ';', &     ! end of 2(A,A)
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;', &    ! end of 8(I7,A)
      & CustomTimingFESolverPreLoad, ';', &
      & CustomTimingOdeSolverPreLoad, ';', &
      & CustomTimingParabolicSolverPreLoad, ';', &
      & CustomTimingFileOutputUserPreLoad, ';', &
      & TimingExportEMGUser, ';', &
      & TimingExportEMGSystem, ';', &
      & CustomTimingFileOutputUser, ';', &
      & CustomTimingFileOutputSystem, ';', &
      & CustomTimingFileOutputSystemPreLoad, ';', &   ! end of 9(F8.3,A)
      & MonodomainSolverId, ';', &
      & MonodomainPreconditionerId, ';', &
      & ODESolverId, ';', &                           ! end of 3(I2,A)
      & NumberOfElementsInAtomX, ';', &
      & NumberOfElementsInAtomY, ';', &
      & NumberOfElementsInAtomZ, ';', &
      & NumberOfSubdomainsX, ';', &
      & NumberOfSubdomainsY, ';', &
      & NumberOfSubdomainsZ, ';', &                          ! end of 6(I8,A)
      & ModelType, ';', &                                     ! I1
      & TRIM(ScenarioName)
  ENDIF

  CLOSE(unit=123)
END SUBROUTINE WriteTimingFile

SUBROUTINE WriteCustomProfilingFile()
  CHARACTER(len=256) :: Filename = "profiling."
  CHARACTER(len=20)  :: ComputationalNodeNumberStr
  CHARACTER(len=100) :: Hostname

  CALL HOSTNM(Hostname)

  ! create filename
  WRITE(ComputationalNodeNumberStr, '(I0.5)') ComputationalNodeNumber     ! convert ComputationalNodeNumber to string
  Filename = TRIM(Filename) // TRIM(ComputationalNodeNumberStr) // ".txt"

  OPEN(unit=123, file=Filename, iostat=stat, access='append')
  IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

  ! Write header to file
  WRITE(123,'(3A,I5)') '# CustomProfiling at ',GetTimeStamp(), ", number of computational nodes: ", NumberOfComputationalNodes
  WRITE(123,'(A)') '# Hostname, NProc, X,Y,Z,F, nFE, nM, tEnd'
  WRITE(123,'(2A,7(I11,A),(F8.3,A))')  &
    & TRIM(Hostname(1:22)), ';', &
    & NumberOfComputationalNodes, ';', &
    & NumberGlobalXElements, ';', &
    & NumberGlobalYElements, ';', &
    & NumberGlobalZElements, ';', &
    & NumberOfInSeriesFibres, ';', &
    & NumberOfElementsFE, ';', &
    & NumberOfElementsM, ';', &
    & TimeStop, ';'

  ! write info
  WRITE(123,*) cmfe_CustomProfilingGetInfo(Err)

  CLOSE(unit=123)
END SUBROUTINE WriteCustomProfilingFile

SUBROUTINE HandleSolverInfo(TimeStep)
  REAL(CMISSRP), INTENT(IN) :: TimeStep
  CHARACTER(len=256) :: Filename = "iterations.csv"
  LOGICAL :: FileExists

  RETURN      ! disable solver info for the moment
  
  IF (ComputationalNodeNumber == 0) THEN

    CALL cmfe_CustomSolverInfoGet(CustomSolverConvergenceReasonParabolic, CustomSolverConvergenceReasonNewton, &
      & CustomSolverNumberIterationsParabolic, CustomSolverNumberIterationsParabolicMin, CustomSolverNumberIterationsParabolicMax, &
      & CustomSolverNumberIterationsNewton, CustomSolverNumberIterationsNewtonMin, CustomSolverNumberIterationsNewtonMax, Err)

    ! Check if file exists
    INQUIRE(file=Filename, exist=FileExists)

    ! Write Comment in first line if file does not yet exist
    IF(.NOT. FileExists) THEN
      OPEN(unit=123, file=Filename, iostat=stat)
      IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'
      WRITE(123,'(A)') '# Timestep; Parabolic reason; Newton reason; parabolic n. iter; min; max; newton n. iter; min; max '
      WRITE(123,'(A)') '# Parabolic reason: 1=RTOL_NORMAL, 2=RTOL, 3=ATOL, 4=ITS, 5=CG_NEG_CURVE, 6=CG_CONSTRAINED, ' // &
        & '7=STEP_LENGTH, 8=HAPPY_BREAKDOWN, 9=ATOL_NORMAL'
      WRITE(123,'(A)') '# Newton reason: 2=FNORM_ABS ||F||<atol, 3=FNORM_RELATIVE ||F|| < rtol*||F_initial||, ' // &
        & '4=SNORM_RELATIVE Newton computed step size small || delta x || < stol || x ||, ' // &
        & '5=ITS, 7=TR_DELTA'
      CLOSE(unit=123)
    ENDIF

    ! Write line in file
    OPEN(unit=123, file=Filename, iostat=stat, access='append')
    IF (stat /= 0 ) PRINT*, 'Failed to open File \"'// TRIM(Filename) // '\" for writing!.'

    WRITE(123,"(F0.8,A,8(I7,A))") &
      & TimeStep, ';', &
      & CustomSolverConvergenceReasonParabolic, ';', &
      & CustomSolverConvergenceReasonNewton, ';', &
      & CustomSolverNumberIterationsParabolic, ';', &
      & CustomSolverNumberIterationsParabolicMin,';',  &
      & CustomSolverNumberIterationsParabolicMax, ';', &
      & CustomSolverNumberIterationsNewton, ';', &
      & CustomSolverNumberIterationsNewtonMin, ';', &
      & CustomSolverNumberIterationsNewtonMax, ';'

    CLOSE(unit=123)
  ENDIF

END SUBROUTINE HandleSolverInfo

SUBROUTINE gdbParallelDebuggingBarrier()
  INTEGER(CMISSIntg) :: Gdb_Resume
  INTEGER(CMISSIntg) :: MPI_IERROR, MPI_COMM_WORLD
  INTEGER(CMISSIntg) :: ComputationalNodeNumber, NumberOfComputationalNodes
  Gdb_Resume = 0

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ComputationalNodeNumber,MPI_IERROR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NumberOfComputationalNodes,MPI_IERROR)

  IF (NumberOfComputationalNodes > 1) THEN
    PRINT*, "Node ", ComputationalNodeNumber, ", UID ",GETPID()," is waiting for Gdb_Resume=", Gdb_Resume &
      & , " to become 1 " // NEW_LINE('A') // "sudo gdb cuboid ",GETPID(), NEW_LINE('A') //"select-frame 2" // &
      & NEW_LINE('A') // "set var gdb_resume = 1" // NEW_LINE('A') // &
      & "info locals" // NEW_LINE('A') // "next"
    DO WHILE (Gdb_Resume == 0)
      CALL Sleep(1)
    ENDDO
    PRINT*, "Node ", ComputationalNodeNumber, " resumes because gdb_resume=", Gdb_Resume, "."
  ENDIF

END SUBROUTINE gdbParallelDebuggingBarrier

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE


