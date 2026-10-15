!
!////////////////////////////////////////////////////////////////////////
!
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
#include "Includes.h"
module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f
   
      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#ifdef FLOW
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(inout)       :: Q(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         use VariableConversion, only: GetGradientValues_f
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(inout)       :: U(NGRAD)
         procedure(GetGradientValues_f)     :: GetGradients
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedGradVars_f


      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)    :: x(NDIM)
         real(kind=RP), intent(in)    :: t
         real(kind=RP), intent(in)    :: nHat(NDIM)
         real(kind=RP), intent(in)    :: Q(NCONS)
         real(kind=RP), intent(in)    :: U_x(NGRAD)
         real(kind=RP), intent(in)    :: U_y(NGRAD)
         real(kind=RP), intent(in)    :: U_z(NGRAD)
         real(kind=RP), intent(inout) :: flux(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedNeumann_f

#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors & 
#ifdef FLOW
         , thermodynamics_ &
         , dimensionless_  &
         , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
         , multiphase_ &
#endif
      )
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         use fluiddata
         use physicsstorage
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
      END SUBROUTINE UserDefinedPeriodicOperation_f
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
#ifdef ACOUSTIC
,Qbase, Lambbase, Lamb &
#endif
)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NCONS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#ifdef ACOUSTIC
         real(kind=RP),             intent(in)  :: Qbase(NCONSB)
         real(kind=RP),             intent(in)  :: Lambbase(NDIM)
         real(kind=RP),             intent(in)  :: Lamb(NDIM)
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                 , thermodynamics_ &
                                                 , dimensionless_  &
                                                 , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
                                                 , multiphase_ &
#endif
                                                 , monitors, &
                                                   elapsedTime, &
                                                   CPUTime   )
         use SMConstants
         USE HexMeshClass
         use FluidData
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)                        :: mesh
         REAL(KIND=RP)                         :: time
         integer                               :: iter
         real(kind=RP)                         :: maxResidual
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         type(Monitor_t),        intent(in)    :: monitors
         real(kind=RP),             intent(in) :: elapsedTime
         real(kind=RP),             intent(in) :: CPUTime
      END SUBROUTINE UserDefinedFinalize_f

      SUBROUTINE UserDefinedTermination_f
         implicit none
      END SUBROUTINE UserDefinedTermination_f
   end interface
   
end module ProblemFileFunctions

module MeshInterpolationTestField
   use SMConstants
   implicit none

   ! Lamb vector files of the test (same format as the Lamb vector files of the solver)
   character(len=*), parameter :: LAMB_SOURCE_FILE       = "RESULTS/Source_0000000000.Lamb.hsol"
   character(len=*), parameter :: LAMB_INTERPOLATED_FILE = "RESULTS/Source_0000000000.Lamb.Interpolated.hsol"

   ! .true. in the run that sets the initial condition (step 1), .false. after a restart (step 3)
   logical :: isSourceRun = .false.

   contains
!
!     Field of the test: polynomials of degree <= 3 in each direction (5 conservative variables)
!
      pure function MeshInterpolationField(x) result(Q)
         real(kind=RP), intent(in) :: x(NDIM)
         real(kind=RP)             :: Q(5)

         Q(1) = 1.0_RP  + 0.1_RP*x(1)*x(2) + 0.01_RP*x(3)**3
         Q(2) = 0.2_RP*x(1)**3 - 0.1_RP*x(2)*x(3)
         Q(3) = 0.1_RP*x(1)*x(2)**2*x(3)
         Q(4) = 0.05_RP*x(1)**2*x(3) + 0.02_RP*x(2)**3
         Q(5) = 10.0_RP + 0.5_RP*x(1)**3*x(2)**2*x(3)

      end function MeshInterpolationField
!
!     Lamb vector field of the test: polynomials of degree <= 3 in each direction (NDIM components)
!
      pure function MeshInterpolationLambField(x) result(L)
         real(kind=RP), intent(in) :: x(NDIM)
         real(kind=RP)             :: L(NDIM)

         L(1) = 0.3_RP*x(1)*x(2)*x(3) - 0.2_RP*x(3)**2
         L(2) = 1.0_RP + 0.1_RP*x(1)**3*x(3)
         L(3) = 0.05_RP*x(2)**3 - 0.4_RP*x(1)*x(2)

      end function MeshInterpolationLambField
!
!     IO offsets of the elements (accumulated number of nodes in global element order)
!
      subroutine MeshInterpolationOffsets(mesh, offsets)
         use HexMeshClass
         class(HexMesh)      , intent(in)  :: mesh
         integer, allocatable, intent(out) :: offsets(:)

         integer, allocatable :: sizes(:)
         integer              :: eID

         allocate(sizes(mesh % no_of_allElements), offsets(mesh % no_of_elements))
         sizes = 0
         do eID = 1, mesh % no_of_elements
            sizes(mesh % elements(eID) % globID) = product(mesh % elements(eID) % Nxyz + 1)
         end do
         do eID = 1, mesh % no_of_elements
            offsets(eID) = sum(sizes(1:mesh % elements(eID) % globID - 1))
         end do

      end subroutine MeshInterpolationOffsets
!
!     Write the Lamb vector field of the test on the nodes of the mesh.
!     The file has the format of the .hsol files (see SolutionFile.f90). It is written here with
!     stream IO because the procedures of SolutionFile cannot be linked in the problem file library
!
      subroutine WriteMeshInterpolationLambFile(mesh)
         use HexMeshClass
         use SolutionFile, only: SOLUTION_FILE, POS_INIT_DATA, BEGINNING_DATA, NO_OF_SAVED_REFS, SOLFILE_STR_LEN
         class(HexMesh), intent(in) :: mesh

         integer, parameter             :: END_OF_FILE = 99
         integer                        :: fid, eID, i, j, k
         integer, allocatable           :: offsets(:)
         integer(kind=AddrInt)          :: pos, fileSize
         character(len=SOLFILE_STR_LEN) :: title
         real(kind=RP)                  :: refs(NO_OF_SAVED_REFS)
         real(kind=RP), allocatable     :: L(:,:,:,:)

         ! Header: name, file type, nodes, number of elements, iteration, time, reference values, terminator
         title = LAMB_SOURCE_FILE
         refs  = 0.0_RP
         open(newunit=fid, file=LAMB_SOURCE_FILE, action="write", status="replace", form="unformatted", access="stream")
         write(fid, pos=1) title, SOLUTION_FILE, mesh % nodeType, mesh % no_of_allElements, 0, 0.0_RP, refs, BEGINNING_DATA

         ! Data of each element: rank, dimensions and array
         call MeshInterpolationOffsets(mesh, offsets)
         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
            allocate(L(NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            do k = 0, e % Nxyz(3);  do j = 0, e % Nxyz(2);  do i = 0, e % Nxyz(1)
               L(:,i,j,k) = MeshInterpolationLambField(e % geom % x(:,i,j,k))
            end do;                 end do;                 end do
            pos = POS_INIT_DATA + (e % globID - 1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*offsets(eID)*SIZEOF_RP
            write(fid, pos=pos) 4, NDIM, e % Nxyz + 1
            write(fid) L
            deallocate(L)
            end associate
         end do

         ! End of file indicator
         inquire(unit=fid, size=fileSize)
         write(fid, pos=fileSize+1) END_OF_FILE
         close(fid)

      end subroutine WriteMeshInterpolationLambFile
!
!     Maximum error of the interpolated Lamb vector file with respect to the exact field.
!     A file that does not match the mesh (type, nodes, elements or orders) gives huge(1.0)
!
      function MeshInterpolationLambError(mesh) result(maxError)
         use HexMeshClass
         use SolutionFile, only: SOLUTION_FILE, POS_INIT_DATA, POS_FILETYPE
         class(HexMesh), intent(in) :: mesh
         real(kind=RP)              :: maxError

         integer                    :: fid, eID, i, j, k, rank, dims(4), header(3)
         integer, allocatable       :: offsets(:)
         integer(kind=AddrInt)      :: pos
         logical                    :: exists
         real(kind=RP), allocatable :: L(:,:,:,:)

         maxError = huge(1.0_RP)

         inquire(file=LAMB_INTERPOLATED_FILE, exist=exists)
         if ( .not. exists ) then
            write(STD_OUT,'(A,A)') "Missing file ", LAMB_INTERPOLATED_FILE
            return
         end if

         ! Header: file type, nodes and number of elements
         open(newunit=fid, file=LAMB_INTERPOLATED_FILE, action="read", status="old", form="unformatted", access="stream")
         read(fid, pos=POS_FILETYPE) header
         if ( any(header .ne. [SOLUTION_FILE, mesh % nodeType, mesh % no_of_allElements]) ) then
            write(STD_OUT,'(A,3I6)') "Unexpected header of the interpolated Lamb vector file: ", header
            close(fid)
            return
         end if

         call MeshInterpolationOffsets(mesh, offsets)

         maxError = 0.0_RP
         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
            pos = POS_INIT_DATA + (e % globID - 1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*offsets(eID)*SIZEOF_RP
            read(fid, pos=pos) rank
            read(fid) dims
            if ( (rank .ne. 4) .or. any(dims .ne. [NDIM, e % Nxyz + 1]) ) then
               maxError = huge(1.0_RP)
               exit
            end if
            allocate(L(NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            read(fid) L
            do k = 0, e % Nxyz(3);  do j = 0, e % Nxyz(2);  do i = 0, e % Nxyz(1)
               maxError = max(maxError, maxval(abs(L(:,i,j,k) - MeshInterpolationLambField(e % geom % x(:,i,j,k)))))
            end do;                 end do;                 end do
            deallocate(L)
            end associate
         end do
         close(fid)

      end function MeshInterpolationLambError

end module MeshInterpolationTestField

         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            use MeshInterpolationTestField
            implicit none
            class(hexmesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

!
!           ---------------------------------------
!           Navier-Stokes initial condition of the test
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refvalues_ % AOAtheta*(pi/180.0_RP)
            phi   = refvalues_ % AOAphi*(pi/180.0_RP)
      
!           Polynomial field (degree <= 3 in each direction), exactly represented with N >= 3
!           in the straight-sided elements of the source mesh. It is interpolated with
!           horses.tools and checked in UserDefinedFinalize after the restart on the other mesh
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx
                  mesh % elements(eID) % storage % q(:,i,j,k) = MeshInterpolationField(mesh % elements(eID) % geom % x(:,i,j,k))
               end do;        end do;        end do
               end associate
            end do
!
!           Lamb vector file of the source mesh, also interpolated with horses.tools
!
            isSourceRun = .true.
            call WriteMeshInterpolationLambFile(mesh)

            end associate
#endif
!
!           ------------------------------------------------------
!           Incompressible Navier-Stokes default initial condition
!           ------------------------------------------------------
!
#if defined(INCNS)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP, 1.0_RP,0.0_RP,0.0_RP,0.0_RP] 
               end do;        end do;        end do
               end associate
            end do
#endif

!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#ifdef CAHNHILLIARD
            call random_seed()
         
            do eid = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eid) % Nxyz(1), &
                          Ny => mesh % elements(eid) % Nxyz(2), &
                          Nz => mesh % elements(eid) % Nxyz(3) )
               associate(e => mesh % elements(eID) % storage)
               call random_number(e % c) 
               e % c = 2.0_RP * (e % c - 0.5_RP)
               end associate
               end associate
            end do
#endif
!
!           -----------------------------------
!           Acoustic default initial condition
!           -----------------------------------
!
#if defined(ACOUSTIC)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  ! mesh % elements(eID) % storage % q(:,i,j,k) = [0.0_RP, 0.0_RP,0.0_RP,0.0_RP,0.0_RP] 
                  mesh % elements(eID) % storage % q(:,i,j,k) = 0.0_RP 
               end do;        end do;        end do
               end associate
            end do
#endif
         end subroutine UserDefinedInitialCondition

         subroutine UserDefinedState1(x, t, nHat, Q & 
#ifdef FLOW
            , thermodynamics_ &
            , dimensionless_  &
            , refValues_ & 
#endif
#ifdef CAHNHILLIARD
            , multiphase_ &
#endif
            )
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
#ifdef FLOW
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         end subroutine UserDefinedState1

#ifdef FLOW

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors &
#ifdef FLOW
         , thermodynamics_ &
         , dimensionless_  &
         , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
         , multiphase_ &
#endif
      )
!
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            use fluiddata
            use physicsstorage
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
#endif
#ifdef ACOUSTIC
, Qbase, Lambbase, Lamb &
#endif
)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#ifdef ACOUSTIC
            real(kind=RP),             intent(in)  :: Qbase(NCONSB)
            real(kind=RP),             intent(in)  :: Lambbase(NDIM)
            real(kind=RP),             intent(in)  :: Lamb(NDIM)
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            use MeshInterpolationTestField
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
#if defined(NAVIERSTOKES)
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName = "Mesh interpolation (tools)"
            REAL(KIND=RP)                      :: maxError, maxErrorLamb
            INTEGER                            :: eID, i, j, k
            TYPE(FTAssertionsManager), POINTER :: sharedManager
!
!           The solution must be the polynomial field of the initial condition. After the
!           restart from the interpolated file (different mesh, polynomial order and nodes)
!           this checks the interpolation of horses.tools ("mesh interpolation")
!
            maxError = 0.0_RP
            do eID = 1, mesh % no_of_elements
               associate( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3);  do j = 0, e % Nxyz(2);  do i = 0, e % Nxyz(1)
                  maxError = max(maxError, maxval(abs(e % storage % Q(:,i,j,k) - MeshInterpolationField(e % geom % x(:,i,j,k)))))
               end do;                 end do;                 end do
               end associate
            end do

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = 1.0_RP, &
                               actualValue   = maxError + 1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Maximum error of the interpolated solution")
!
!           After the restart (step 3), check also the interpolated Lamb vector file
!
            maxErrorLamb = 0.0_RP
            if ( .not. isSourceRun ) then
               maxErrorLamb = MeshInterpolationLambError(mesh)
               CALL FTAssertEqual(expectedValue = 1.0_RP, &
                                  actualValue   = maxErrorLamb + 1.0_RP, &
                                  tol           = 1.d-11, &
                                  msg           = "Maximum error of the interpolated Lamb vector")
            end if

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)

            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "Maximum error (solution):    ", maxError
               WRITE(6,*) "Maximum error (Lamb vector): ", maxErrorLamb
                error stop 99
            END IF
            WRITE(6,*)

            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
