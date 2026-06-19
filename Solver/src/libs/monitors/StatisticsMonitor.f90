#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
module StatisticsMonitor
   use SMConstants
   use HexMeshClass
   use StorageClass
   use Utilities, only: GreatestCommonDivisor
#ifdef _HAS_MPI_
   use mpi
#endif

   private
   public     StatisticsMonitor_t, U, V, W, UU, VV, WW, UV, UW, VW
   public     NO_OF_VARIABLES_Sij, NO_OF_VARIABLES
!
!  Commands for the parameter file
!  -------------------------------
#define START_COMMAND "@start"
#define PAUSE_COMMAND "@pause"
#define STOP_COMMAND "@stop"
#define RESET_COMMAND "@reset"
#define DUMP_COMMAND "@dump"
!
!  Statistics monitor states
!  -------------------------
#define OFF 0
#define WAITING_STATE 1
#define STOP_STATE 2
#define ACTIVE_STATE 3
!
!  Statistics monitor content
!  --------------------------
   integer            :: NO_OF_VARIABLES
   integer, parameter :: NO_OF_VARIABLES_Sij = 9
   integer, parameter ::  U  = 1
   integer, parameter ::  V  = 2
   integer, parameter ::  W  = 3
   integer, parameter ::  UU = 4
   integer, parameter ::  VV = 5
   integer, parameter ::  WW = 6
   integer, parameter ::  UV = 7
   integer, parameter ::  UW = 8
   integer, parameter ::  VW = 9

   type StatisticsMonitor_t
      integer        :: state
      integer        :: sampling_interval
      integer        :: dump_interval
      integer        :: reset_interval
      integer        :: starting_iteration
      real(kind=RP)  :: starting_time
      integer        :: no_of_samples
      logical        :: saveGradients
      logical        :: saveLambVector
      logical        :: saveSoundVelocitySquared
      logical        :: saveGradSoundVelocitySquared
      logical        :: savePressureSquared
      contains
         procedure   :: Construct    => StatisticsMonitor_Construct
         procedure   :: Update       => StatisticsMonitor_Update
         procedure   :: UpdateValues => StatisticsMonitor_UpdateValues
         procedure   :: GetState     => StatisticsMonitor_GetState
         procedure   :: WriteLabel   => StatisticsMonitor_WriteLabel
         procedure   :: WriteValue   => StatisticsMonitor_WriteValue
         procedure   :: WriteFile    => StatisticsMonitor_WriteFile
   end type StatisticsMonitor_t
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine StatisticsMonitor_Construct(self, mesh, saveGradients, saveLambVector, saveSoundVelocitySquared, saveGradSoundVelocitySquared, savePressureSquared)
         use ParamfileRegions
         use PhysicsStorage, only: NCONS, NGRAD
         implicit none
         class(StatisticsMonitor_t)    :: self
         class(HexMesh)                :: mesh
         logical, intent(in)           :: saveGradients, saveLambVector, saveSoundVelocitySquared, saveGradSoundVelocitySquared
         logical, intent(in)           :: savePressureSquared
         integer, allocatable          :: Nsample, i0, Ndump, Nreset
         real(kind=RP), allocatable    :: t0
         integer                       :: eID
         character(len=LINE_LENGTH)    :: paramFile

         NO_OF_VARIABLES = NO_OF_VARIABLES_Sij + NCONS
         self % saveGradients = saveGradients
         self % saveLambVector = saveLambVector
         self % saveSoundVelocitySquared = saveSoundVelocitySquared
         self % saveGradSoundVelocitySquared = saveGradSoundVelocitySquared
         self % savePressureSquared = savePressureSquared
         if (saveGradients) NO_OF_VARIABLES = NO_OF_VARIABLES + NGRAD * NDIM
         if (saveLambVector) NO_OF_VARIABLES = NO_OF_VARIABLES + NDIM
         if (saveSoundVelocitySquared) NO_OF_VARIABLES = NO_OF_VARIABLES + 1
         if (saveGradSoundVelocitySquared) NO_OF_VARIABLES = NO_OF_VARIABLES + NDIM
         if (savePressureSquared) NO_OF_VARIABLES = NO_OF_VARIABLES + 1
!
!        Search for the parameters in the case file
!        ------------------------------------------
         call get_command_argument(1, paramFile)
!
!        Search if the #define statistics section is defined
!        ---------------------------------------------------
         if ( findIfStatsAreActive(paramFile) ) then
            self % state = WAITING_STATE

         else
            self % state = OFF
            return

         end if

         call readValueInRegion ( trim ( paramFile ) , "sampling interval"  , Nsample , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "starting iteration" , i0      , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "starting time"      , t0      , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "dump interval"      , Ndump   , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "reset interval"     , Nreset  , "#define statistics" , "#end" )
         
!
!        Check the input data
!        --------------------
         if ( allocated(i0) ) self % starting_iteration = i0      ! Both iteration and time can start the statistics calculation
         if ( allocated(t0) ) self % starting_time = t0           !

         if ( allocated(i0) .and. (.not. allocated(t0)) ) self % starting_time = huge(1.0_RP) ! Statistics controlled by iterations
         if ( allocated(t0) .and. (.not. allocated(i0)) ) self % starting_iteration = huge(1) ! Statistics controlled by time

         if ( (.not. allocated(i0)) .and. (.not. allocated(t0)) ) then  ! 
            self % starting_time = 0.0_RP                               ! Start since the beginning by default
            self % starting_iteration = 0                               !
         end if                                                         !

         if ( allocated(Nsample) ) then            !
            self % sampling_interval = Nsample     !  Set 1 by default
         else                                      !
            self % sampling_interval = 1           ! 
         end if     
         
         if ( allocated(Ndump) ) then              !
            self % dump_interval = Ndump           !  Set huge by default
         else                                      !
            self % dump_interval = huge(1)         ! 
         end if                                    !
         
         if ( allocated(Nreset) ) then                                        !
            if (allocated(Ndump) ) then                                       !
               self % dump_interval  = GreatestCommonDivisor(Nreset,Ndump)    !  Set huge by default
            else                                                              !
               self % dump_interval  = Nreset                                 !
            end if                                                            !
            self % reset_interval = Nreset                                    !
         else                                                                 !
            self % reset_interval = huge(1)                                   ! 
         end if                                                               !
!
!        Initial state: waiting
!        ----------------------
         self % state = WAITING_STATE
         self % no_of_samples = 0
!
!        Construct the data
!        ------------------
         do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
            call e % storage % stats % Construct(NO_OF_VARIABLES, e % Nxyz)
            end associate
         end do

      end subroutine StatisticsMonitor_Construct

      subroutine StatisticsMonitor_WriteFile(self, mesh, iter, t, solution_file)
         use PhysicsStorage
         implicit none
         class(StatisticsMonitor_t) :: self
         class(HexMesh)             :: mesh
         integer, intent(in)        :: iter
         real(kind=RP), intent(in)  :: t
         character(len=*), intent(in)  :: solution_file
         character(len=LINE_LENGTH)    :: fileName
         
         integer :: lblimits, ublimits

         write(fileName,'(A,A)') trim(solution_file),'.stats.hsol'
         if ( self % state .ne. OFF) then
            call mesh % SaveStatistics(iter, t, trim(fileName), self % saveGradients)
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
            ! Save Lamb vector
            lblimits = NO_OF_VARIABLES_Sij + NCONS
            if (self % saveGradients) lblimits = lblimits + NDIM * NGRAD
            lblimits = lblimits + 1
            ublimits = lblimits + NDIM - 1
            write(fileName,'(A,A)') trim(solution_file),'.Lamb.stats.hsol'
            if (self % saveLambVector) call mesh % SaveLambVectorStatistics(iter, t, trim(fileName), lblimits, ublimits)
#endif
#ifdef NAVIERSTOKES
            ! Save sound velocity
            lblimits = NO_OF_VARIABLES_Sij + NCONS
            if (self % saveGradients) lblimits = lblimits + NDIM * NGRAD
            if (self % saveLambVector) lblimits = lblimits + NDIM
            lblimits = lblimits + 1
            ublimits = lblimits + 1 - 1
            write(fileName,'(A,A)') trim(solution_file),'.SoundVelocitySquared.stats.hsol'
            if (self % saveSoundVelocitySquared) call mesh % saveSoundVelocitySquaredStatistics(iter, t, trim(filename), lblimits, ublimits)
            ! Save gradient of sound velocity squared
            lblimits = NO_OF_VARIABLES_Sij + NCONS
            if (self % saveGradients) lblimits = lblimits + NDIM * NGRAD
            if (self % saveLambVector) lblimits = lblimits + NDIM
            if (self % saveSoundVelocitySquared) lblimits = lblimits + 1
            lblimits = lblimits + 1
            ublimits = lblimits + NDIM - 1
            write(fileName,'(A,A)') trim(solution_file),'.GradientSoundVelocitySquared.stats.hsol'
            if (self % saveGradSoundVelocitySquared) call mesh % saveGradientSoundVelocitySquaredStatistics(iter, t, trim(filename), lblimits, ublimits)
#endif
#ifdef ACOUSTIC
            ! Save pressure squared
            lblimits = NO_OF_VARIABLES_Sij + NCONS
            if (self % saveGradients) lblimits = lblimits + NDIM * NGRAD
            lblimits = lblimits + 1
            ublimits = lblimits
            write(fileName,'(A,A)') trim(solution_file),'.PressureSquared.stats.hsol'
            if (self % savePressureSquared) call mesh % SavePressureSquaredStatistics(iter, t, trim(fileName), lblimits, ublimits)
#endif
         end if

      end subroutine StatisticsMonitor_WriteFile

      subroutine StatisticsMonitor_Update(self, mesh, iter, t, solution_file)
         implicit none
         class(StatisticsMonitor_t) :: self
         class(HexMesh)             :: mesh
         integer, intent(in)        :: iter
         real(kind=RP), intent(in)  :: t
         character(len=*), intent(in)  :: solution_file
!
!        ---------------
!        Local variables
!        ---------------
!
         logical     :: reset, dump
         character(len=LINE_LENGTH)    :: fileName
   
         if ( self % state .eq. OFF ) return
         if ( mod(iter, self % sampling_interval) .ne. 0) return
!
!        Read the parameter file to check the state
!        ------------------------------------------
         call self % GetState(reset, dump)
!
!        Dump the contents if requested
!        ------------------------------
         if ( dump .or. ( (mod(iter, self % dump_interval) == 0) .and. (iter > self % starting_iteration .or. t > self % starting_time) ) ) then
            write(fileName,'(A,A,I10.10,A)') trim(solution_file),'.stats.',iter,'.hsol'
            call mesh % SaveStatistics(iter, t, trim(fileName), self % saveGradients)
            write(STD_OUT,'(A,A,A)') '   *** Saving statistics file as "',trim(fileName),'".'
         end if
!
!        Reset the statistics if requested
!        ---------------------------------
         if ( reset .or. ( (mod(iter, self % reset_interval) == 0) .and. (iter > self % starting_iteration) ) ) then
            call mesh % ResetStatistics
            self % no_of_samples = 0
         end if
!
!        Update the state if it is waiting
!        ---------------------------------
         if ( self % state .eq. WAITING_STATE ) then
            if ( (t .gt. self % starting_time) .or. (iter .ge. self % starting_iteration) ) then
               self % state = ACTIVE_STATE
            end if
         end if
            
!        Update the values if the state is "active"
!        ------------------------------------------
         if(self % state .eq. ACTIVE_STATE) call self % UpdateValues(mesh)

      end subroutine StatisticsMonitor_Update

      subroutine StatisticsMonitor_UpdateValues(self, mesh)
         use PhysicsStorage
         use MappedGeometryClass, only: vCross, ComputeVorticity
#if (defined(NAVIERSTOKES) && !(defined(SPALARTALMARAS)) && !(defined(TRANSPORT)))
         use VariableConversion_NS, only: getVelocityGradients_STATE, Pressure
#elif (defined(NAVIERSTOKES) && (defined(SPALARTALMARAS)) )
         use VariableConversion_NSSA, only: getVelocityGradients_STATE, Pressure
#elif (defined(NAVIERSTOKES) && (defined(TRANSPORT)))
         use VariableConversion_NSTPT, only: getVelocityGradients_STATE, Pressure
#elif defined(INCNS)
         use VariableConversion_iNS, only: getVelocityGradients
#elif defined(MULTIPHASE)
         use VariableConversion_MU, only: getVelocityGradients
         use FluidData_MU, only: dimensionless
#endif
         implicit none
         class(StatisticsMonitor_t)    :: self
         class(HexMesh)              :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID
         integer  :: i, j, k
         real(RP) :: ratio, inv_nsamples_plus_1
         real(RP) :: rfactor1, rfactor2
         integer, dimension(8) :: limits
         real(kind=RP), dimension(NDIM) :: vorticity, velocity, LambVector
         real(kind=RP), dimension(NDIM) :: gradvel_x, gradvel_y, gradvel_z, uDivRho
         real(kind=RP) :: press, invRho, invRho2, soundVelocity
         real(kind=RP), dimension(NDIM) :: gradPress, gradSoundVelocity

#ifdef NAVIERSTOKES
         !  if gradients are not saved, limits(2) is equal to limits(5), the latter wont be used
            limits(1) = NO_OF_VARIABLES_Sij + IRHO
            limits(2) = NO_OF_VARIABLES_Sij + NCONS
            limits(3:7) = limits(2)
            if (self % saveGradients) then
               limits(3) = NO_OF_VARIABLES_Sij + NCONS + NGRAD ! U_x
               limits(4) = NO_OF_VARIABLES_Sij + NCONS + 2*NGRAD ! U_y
               limits(5) = NO_OF_VARIABLES_Sij + NCONS + NDIM * NGRAD ! U_z
               limits(6:8) = limits(5)
            end if
            if (self % saveLambVector) limits(6) = limits(5) + NDIM
            if (self % saveSoundVelocitySquared) limits(7) = limits(6) + 1
            if (self % saveGradSoundVelocitySquared) limits(8) = limits(7) + NDIM

            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            do eID = 1, size(mesh % elements)
               associate(e    => mesh % elements(eID), &
                         data => mesh % elements(eID) % storage % stats % data)

               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / e % storage % Q(IRHO,i,j,k)
                  rfactor2 = inv_nsamples_plus_1 / POW2( e % storage % Q(IRHO,i,j,k) )
                  data(U,i,j,k)  = data(U,i,j,k)  * ratio + e % storage % Q(IRHOU,i,j,k) * rfactor1
                  data(V,i,j,k)  = data(V,i,j,k)  * ratio + e % storage % Q(IRHOV,i,j,k) * rfactor1
                  data(W,i,j,k)  = data(W,i,j,k)  * ratio + e % storage % Q(IRHOW,i,j,k) * rfactor1
                  data(UU,i,j,k) = data(UU,i,j,k) * ratio + POW2( e % storage % Q(IRHOU,i,j,k) ) * rfactor2
                  data(VV,i,j,k) = data(VV,i,j,k) * ratio + POW2( e % storage % Q(IRHOV,i,j,k) ) * rfactor2
                  data(WW,i,j,k) = data(WW,i,j,k) * ratio + POW2( e % storage % Q(IRHOW,i,j,k) ) * rfactor2
                  data(UV,i,j,k) = data(UV,i,j,k) * ratio + e % storage % Q(IRHOU,i,j,k) * e % storage % Q(IRHOV,i,j,k) * rfactor2
                  data(UW,i,j,k) = data(UW,i,j,k) * ratio + e % storage % Q(IRHOU,i,j,k) * e % storage % Q(IRHOW,i,j,k) * rfactor2
                  data(VW,i,j,k) = data(VW,i,j,k) * ratio + e % storage % Q(IRHOV,i,j,k) * e % storage % Q(IRHOW,i,j,k) * rfactor2
                  data(limits(1):limits(2),i,j,k) = data(limits(1):limits(2),i,j,k) * ratio + e % storage % Q(:,i,j,k) * inv_nsamples_plus_1
                  if (self % saveGradients) then
                      data(limits(2)+1:limits(3),i,j,k) = data(limits(2)+1:limits(3),i,j,k) * ratio + e % storage % U_x(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(3)+1:limits(4),i,j,k) = data(limits(3)+1:limits(4),i,j,k) * ratio + e % storage % U_y(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(4)+1:limits(5),i,j,k) = data(limits(4)+1:limits(5),i,j,k) * ratio + e % storage % U_z(:,i,j,k) * inv_nsamples_plus_1
                  end if
                  ! Save Lamb vector: u x w
                  if (self % saveLambVector) then
                     call getVelocityGradients_State(e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), gradvel_x, gradvel_y, gradvel_z)
                     call ComputeVorticity(gradvel_x, gradvel_y, gradvel_z, vorticity)
                     velocity = e % storage % Q(IRHOU:IRHOW,i,j,k) / e % storage % Q(IRHO,i,j,k)
                     call vCross(velocity, vorticity, LambVector)
                     data(limits(5)+1:limits(6),i,j,k) = data(limits(5)+1:limits(6),i,j,k) * ratio + LambVector * inv_nsamples_plus_1
                  end if
                  if (self % saveSoundVelocitySquared) then
                     press = Pressure(e % storage % Q(:,i,j,k))
                     data(limits(6)+1:limits(7),i,j,k) = data(limits(6)+1:limits(7),i,j,k) * ratio + thermodynamics % gamma * press / e % storage % Q(IRHO,i,j,k) * inv_nsamples_plus_1
                  end if
                  if (self % saveGradSoundVelocitySquared) then
                     call getVelocityGradients_State(e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), gradvel_x, gradvel_y, gradvel_z)
                     invRho = 1.0_rp / e % storage % Q(IRHO,i,j,k)
                     invRho2 = invRho * invRho
                     uDivRho = invRho * e % storage % Q(IRHOU:IRHOW,i,j,k)
                     press = Pressure(e % storage % Q(:,i,j,k))
                     soundVelocity = sqrt( thermodynamics % gamma * press * invRho )
                     gradPress(1) = thermodynamics % gammaMinus1 * ( e % storage % U_x(IRHOE,i,j,k) - dot_product(uDivRho, gradvel_x) + &
                                 0.5_rp * (e % storage % Q(IRHOU,i,j,k)**2 + e % storage % Q(IRHOV,i,j,k)**2 + e % storage % Q(IRHOW,i,j,k)**2) * invRho2 * e % storage % U_x(IRHO,i,j,k))
                     gradPress(2) = thermodynamics % gammaMinus1 * ( e % storage % U_y(IRHOE,i,j,k) - dot_product(uDivRho, gradvel_y) + &
                                 0.5_rp * (e % storage % Q(IRHOU,i,j,k)**2 + e % storage % Q(IRHOV,i,j,k)**2 + e % storage % Q(IRHOW,i,j,k)**2) * invRho2 * e % storage % U_y(IRHO,i,j,k))
                     gradPress(3) = thermodynamics % gammaMinus1 * ( e % storage % U_z(IRHOE,i,j,k) - dot_product(uDivRho, gradvel_z) + &
                                 0.5_rp * (e % storage % Q(IRHOU,i,j,k)**2 + e % storage % Q(IRHOV,i,j,k)**2 + e % storage % Q(IRHOW,i,j,k)**2) * invRho2 * e % storage % U_z(IRHO,i,j,k))
                     gradSoundVelocity(1) = thermodynamics % gamma * ( gradPress(1) * invRho - press * e % storage % U_x(IRHO,i,j,k) * invRho2 )
                     gradSoundVelocity(2) = thermodynamics % gamma * ( gradPress(2) * invRho - press * e % storage % U_y(IRHO,i,j,k) * invRho2 )
                     gradSoundVelocity(3) = thermodynamics % gamma * ( gradPress(3) * invRho - press * e % storage % U_z(IRHO,i,j,k) * invRho2 )
                     data(limits(7)+1:limits(8),i,j,k) = data(limits(7)+1:limits(8),i,j,k) * ratio + 2.0_rp * soundVelocity * gradSoundVelocity * inv_nsamples_plus_1
                  end if
               end do                  ; end do                   ; end do

               end associate
            end do
#endif 

#ifdef INCNS
         !  if gradients are not saved, limits(2) is equal to limits(5), the latter wont be used
            limits(1) = NO_OF_VARIABLES_Sij + INSRHO
            limits(2) = NO_OF_VARIABLES_Sij + NCONS
            limits(3:6) = limits(2)
            if (self % saveGradients) then
               limits(3) = NO_OF_VARIABLES_Sij + NCONS + NGRAD ! U_x
               limits(4) = NO_OF_VARIABLES_Sij + NCONS + 2*NGRAD ! U_y
               limits(5) = NO_OF_VARIABLES_Sij + NCONS + NDIM * NGRAD ! U_z
               limits(6) = limits(5)
            end if
            if (self % saveLambVector) limits(6) = limits(5) + NDIM

            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            do eID = 1, size(mesh % elements)
               associate(e    => mesh % elements(eID), &
                         data => mesh % elements(eID) % storage % stats % data)

               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / e % storage % Q(INSRHO,i,j,k)
                  rfactor2 = inv_nsamples_plus_1 / POW2( e % storage % Q(INSRHO,i,j,k) )
                  data(U,i,j,k)  = data(U,i,j,k)  * ratio + e % storage % Q(INSRHOU,i,j,k) * rfactor1
                  data(V,i,j,k)  = data(V,i,j,k)  * ratio + e % storage % Q(INSRHOV,i,j,k) * rfactor1
                  data(W,i,j,k)  = data(W,i,j,k)  * ratio + e % storage % Q(INSRHOW,i,j,k) * rfactor1
                  data(UU,i,j,k) = data(UU,i,j,k) * ratio + POW2( e % storage % Q(INSRHOU,i,j,k) ) * rfactor2
                  data(VV,i,j,k) = data(VV,i,j,k) * ratio + POW2( e % storage % Q(INSRHOV,i,j,k) ) * rfactor2
                  data(WW,i,j,k) = data(WW,i,j,k) * ratio + POW2( e % storage % Q(INSRHOW,i,j,k) ) * rfactor2
                  data(UV,i,j,k) = data(UV,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOV,i,j,k) * rfactor2
                  data(UW,i,j,k) = data(UW,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor2
                  data(VW,i,j,k) = data(VW,i,j,k) * ratio + e % storage % Q(INSRHOV,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor2
                  data(limits(1):limits(2),i,j,k) = data(limits(1):limits(2),i,j,k) * ratio + e % storage % Q(:,i,j,k) * inv_nsamples_plus_1
                  if (self % saveGradients) then
                      data(limits(2)+1:limits(3),i,j,k) = data(limits(2)+1:limits(3),i,j,k) * ratio + e % storage % U_x(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(3)+1:limits(4),i,j,k) = data(limits(3)+1:limits(4),i,j,k) * ratio + e % storage % U_y(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(4)+1:limits(5),i,j,k) = data(limits(4)+1:limits(5),i,j,k) * ratio + e % storage % U_z(:,i,j,k) * inv_nsamples_plus_1
                  end if 
                  ! Save Lamb vector: u x w
                  if (self % saveLambVector) then
                     call getVelocityGradients(e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), gradvel_x, gradvel_y, gradvel_z)
                     call ComputeVorticity(gradvel_x, gradvel_y, gradvel_z, vorticity)
                     velocity = e % storage % Q(INSRHOU:INSRHOW,i,j,k) / e % storage % Q(INSRHO,i,j,k)
                     call vCross(velocity, vorticity, LambVector)
                     data(limits(5)+1:limits(6),i,j,k) = data(limits(5)+1:limits(6),i,j,k) * ratio + LambVector * inv_nsamples_plus_1
                  end if
               end do                  ; end do                   ; end do

               end associate
            end do
#endif 

#ifdef MULTIPHASE
         !  if gradients are not saved, limits(2) is equal to limits(5), the latter wont be used
            limits(1) = NO_OF_VARIABLES_Sij + IMC
            limits(2) = NO_OF_VARIABLES_Sij + NCONS + 1 ! Also the density
            limits(3:6) = limits(2)
            if (self % saveGradients) then
               limits(3) = NO_OF_VARIABLES_Sij + NCONS + 1 + NGRAD ! U_x
               limits(4) = NO_OF_VARIABLES_Sij + NCONS + 1 + 2*NGRAD ! U_y
               limits(5) = NO_OF_VARIABLES_Sij + NCONS + 1 + NDIM * NGRAD ! U_z
               limits(6) = limits(5)
            end if
            if (self % saveLambVector) limits(6) = limits(5) + NDIM

            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            do eID = 1, size(mesh % elements)
               associate(e    => mesh % elements(eID), &
                         data => mesh % elements(eID) % storage % stats % data)

               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / sqrt( e % storage % rho(i,j,k) )
                  rfactor2 = inv_nsamples_plus_1 / e % storage % rho(i,j,k)
                  data(U,i,j,k)  = data(U,i,j,k)  * ratio + e % storage % Q(IMSQRHOU,i,j,k) * rfactor1
                  data(V,i,j,k)  = data(V,i,j,k)  * ratio + e % storage % Q(IMSQRHOV,i,j,k) * rfactor1
                  data(W,i,j,k)  = data(W,i,j,k)  * ratio + e % storage % Q(IMSQRHOW,i,j,k) * rfactor1
                  data(UU,i,j,k) = data(UU,i,j,k) * ratio + POW2( e % storage % Q(IMSQRHOU,i,j,k) ) * rfactor2
                  data(VV,i,j,k) = data(VV,i,j,k) * ratio + POW2( e % storage % Q(IMSQRHOV,i,j,k) ) * rfactor2
                  data(WW,i,j,k) = data(WW,i,j,k) * ratio + POW2( e % storage % Q(IMSQRHOW,i,j,k) ) * rfactor2
                  data(UV,i,j,k) = data(UV,i,j,k) * ratio + e % storage % Q(IMSQRHOU,i,j,k) * e % storage % Q(IMSQRHOV,i,j,k) * rfactor2
                  data(UW,i,j,k) = data(UW,i,j,k) * ratio + e % storage % Q(IMSQRHOU,i,j,k) * e % storage % Q(IMSQRHOW,i,j,k) * rfactor2
                  data(VW,i,j,k) = data(VW,i,j,k) * ratio + e % storage % Q(IMSQRHOV,i,j,k) * e % storage % Q(IMSQRHOW,i,j,k) * rfactor2
                  data(limits(1):limits(2)-1,i,j,k) = data(limits(1):limits(2)-1,i,j,k) * ratio + e % storage % Q(:,i,j,k) * inv_nsamples_plus_1
                  data(limits(2)-1:limits(2),i,j,k) = data(limits(2)-1:limits(2),i,j,k) * ratio + e % storage % rho(i,j,k) * inv_nsamples_plus_1
                  if (self % saveGradients) then
                      data(limits(2)+1:limits(3),i,j,k) = data(limits(2)+1:limits(3),i,j,k) * ratio + e % storage % U_x(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(3)+1:limits(4),i,j,k) = data(limits(3)+1:limits(4),i,j,k) * ratio + e % storage % U_y(:,i,j,k) * inv_nsamples_plus_1
                      data(limits(4)+1:limits(5),i,j,k) = data(limits(4)+1:limits(5),i,j,k) * ratio + e % storage % U_z(:,i,j,k) * inv_nsamples_plus_1
                  end if
                  ! Save Lamb vector: u x w
                  if (self % saveLambVector) then
                     call getVelocityGradients(e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), dimensionless, gradvel_x, gradvel_y, gradvel_z)
                     call ComputeVorticity(gradvel_x, gradvel_y, gradvel_z, vorticity)
                     mesh % elements(eID) % storage % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2)) * mesh % elements(eID) % storage % Q(IMC,:,:,:)
                     mesh % elements(eID) % storage % rho = min(max(mesh % elements(eID) % storage % rho, dimensionless % rho_min),dimensionless % rho_max)
                     velocity = e % storage % Q(IMSQRHOU:IMSQRHOW,i,j,k) / sqrt( e % storage % rho(i,j,k) )
                     call vCross(velocity, vorticity, LambVector)
                     data(limits(5)+1:limits(6),i,j,k) = data(limits(5)+1:limits(6),i,j,k) * ratio + LambVector * inv_nsamples_plus_1
                  end if
               end do                  ; end do                   ; end do

               end associate
            end do
#endif

#ifdef ACOUSTIC
         !  if gradients are not saved, limits(2) is equal to limits(5), the latter wont be used
            limits(1) = NO_OF_VARIABLES_Sij + ICAARHO
            limits(2) = NO_OF_VARIABLES_Sij + NCONS
            limits(3:7) = limits(2)
            ! if (self % saveGradients) then
            !    limits(3) = NO_OF_VARIABLES_Sij + NCONS + NGRAD ! U_x
            !    limits(4) = NO_OF_VARIABLES_Sij + NCONS + 2*NGRAD ! U_y
            !    limits(5) = NO_OF_VARIABLES_Sij + NCONS + NDIM * NGRAD ! U_z
            !    limits(6:8) = limits(5)
            ! end if
            if (self % savePressureSquared) limits(6) = limits(5) + 1
            
            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            do eID = 1, size(mesh % elements)
               associate(e    => mesh % elements(eID), &
                         data => mesh % elements(eID) % storage % stats % data)

               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  data(U,i,j,k)  = data(U,i,j,k)  * ratio + e % storage % Q(ICAAU,i,j,k) * inv_nsamples_plus_1
                  data(V,i,j,k)  = data(V,i,j,k)  * ratio + e % storage % Q(ICAAV,i,j,k) * inv_nsamples_plus_1
                  data(W,i,j,k)  = data(W,i,j,k)  * ratio + e % storage % Q(ICAAW,i,j,k) * inv_nsamples_plus_1
                  data(UU,i,j,k) = data(UU,i,j,k) * ratio + POW2( e % storage % Q(ICAAU,i,j,k) ) * inv_nsamples_plus_1
                  data(VV,i,j,k) = data(VV,i,j,k) * ratio + POW2( e % storage % Q(ICAAV,i,j,k) ) * inv_nsamples_plus_1
                  data(WW,i,j,k) = data(WW,i,j,k) * ratio + POW2( e % storage % Q(ICAAW,i,j,k) ) * inv_nsamples_plus_1
                  data(UV,i,j,k) = data(UV,i,j,k) * ratio + e % storage % Q(ICAAU,i,j,k) * e % storage % Q(ICAAV,i,j,k) * inv_nsamples_plus_1
                  data(UW,i,j,k) = data(UW,i,j,k) * ratio + e % storage % Q(ICAAU,i,j,k) * e % storage % Q(ICAAW,i,j,k) * inv_nsamples_plus_1
                  data(VW,i,j,k) = data(VW,i,j,k) * ratio + e % storage % Q(ICAAV,i,j,k) * e % storage % Q(ICAAW,i,j,k) * inv_nsamples_plus_1
                  data(limits(1):limits(2),i,j,k) = data(limits(1):limits(2),i,j,k) * ratio + e % storage % Q(:,i,j,k) * inv_nsamples_plus_1
                  ! if (self % saveGradients) then
                  !     data(limits(2)+1:limits(3),i,j,k) = data(limits(2)+1:limits(3),i,j,k) * ratio + e % storage % U_x(:,i,j,k) * inv_nsamples_plus_1
                  !     data(limits(3)+1:limits(4),i,j,k) = data(limits(3)+1:limits(4),i,j,k) * ratio + e % storage % U_y(:,i,j,k) * inv_nsamples_plus_1
                  !     data(limits(4)+1:limits(5),i,j,k) = data(limits(4)+1:limits(5),i,j,k) * ratio + e % storage % U_z(:,i,j,k) * inv_nsamples_plus_1
                  ! end if
                  if (self % savePressureSquared) then
                     data(limits(5)+1:limits(6),i,j,k) = data(limits(5)+1:limits(6),i,j,k) * ratio + POW2( e % storage % Q(ICAAP,i,j,k) ) * inv_nsamples_plus_1
                  end if
               end do                  ; end do                   ; end do

               end associate
            end do
#endif 

         self % no_of_samples = self % no_of_samples + 1
         
      end subroutine StatisticsMonitor_UpdateValues

      subroutine StatisticsMonitor_GetState(self, reset, dump)
         use ParamfileRegions
         use MPI_Process_Info
         implicit none
         class(StatisticsMonitor_t)    :: self
         logical, intent(out)          :: reset, dump
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: fID, io, ierr
         integer                       :: no_of_lines
         character(len=LINE_LENGTH)    :: paramFile
         character(len=LINE_LENGTH)    :: line
         logical                       :: isInside
         logical                       :: hasStart, hasPause, hasStop, hasDump, hasReset
         integer                       :: rstInt, dumpInt

         reset = .false.
         dump  = .false.
!
!        Only root will read, and it will broadcast the result to the rest
!        -----------------------------------------------------------------
         if ( MPI_Process % isRoot ) then
   
            call get_command_argument(1, paramFile)
!   
!           Open the file and get to the statistics zone
!           --------------------------------------------
            open( newunit = fID, file = trim(paramFile), status = "old", action = "read" )
      
            isInside = .false.
   
            no_of_lines = 0
   
            hasStart = .false.
            hasPause = .false.
            hasStop  = .false.
            hasReset = .false.
            hasDump  = .false.
   
            do
               read(fID,'(A)',iostat=io) line
               if ( io .ne. 0 ) exit
               no_of_lines = no_of_lines + 1 
!   
!              Check whether the line is in or out of the statistics region
!              ------------------------------------------------------------
               if ( index(getSquashedLine(trim(adjustl(line))), "#definestatistics" ) .ne. 0) then
                  isInside = .true.
                  cycle
               end if
   
               if ( isInside .and. (index(getSquashedLine(trim(adjustl(line))),"#end") .ne. 0) ) then
                  isInside = .false.
                  cycle
               end if
!   
!              Cycle if it is outside
!              ----------------------        
               if ( .not. isInside ) cycle
!   
!              It is inside, check whether the marks are present
!              -------------------------------------------------
               select case ( trim(adjustl(line)) )
               
               case(START_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasStart = .true.
   
               case(PAUSE_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasPause = .true.
   
               case(STOP_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasStop = .true.
   
               case(RESET_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasReset = .true.
   
               case(DUMP_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasDump = .true.
   
               end select
   
            end do
!   
!           Update the status of the statistics monitor
!           -------------------------------------------
            reset = hasReset
            dump  = hasDump
            
            if ( (self % state .eq. WAITING_STATE) .or. (self % state .eq. ACTIVE_STATE) ) then
               if ( hasStop ) then
                  dump = .true.
                  reset = .true.
                  self % state = STOP_STATE
               end if
   
               if ( hasPause ) then
                  self % state = STOP_STATE
               end if
   
               if ( hasStart ) then
                  self % state = ACTIVE_STATE
               end if
   
            elseif ( self % state .eq. STOP_STATE ) then
               if ( hasStart )  self % state = ACTIVE_STATE
           
            end if
   
            close(fID)
!   
!           Rewrite the control file to remove the controllers      
!           --------------------------------------------------
            if ( hasStart .or. hasPause .or. hasStop .or. hasReset .or. hasDump ) call rewriteControlFile(paramFile, no_of_lines)
         end if
!
!        Broadcast the result to the rest of the processes
!        -------------------------------------------------
         if ( MPI_Process % DoMPIAction ) then
#ifdef _HAS_MPI_
            if ( reset ) then
               rstInt = 1
            else
               rstInt = 0
            end if

            if ( dump ) then
               dumpInt = 1
            else
               dumpInt = 0
            end if

            call mpi_bcast(rstInt, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif
         end if
         

      end subroutine StatisticsMonitor_GetState

      logical function findIfStatsAreActive(paramFile)
         use ParamfileRegions
         implicit none
         character(len=*), intent(in)     :: paramFile
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fID, io
         character(len=LINE_LENGTH)    :: line

         open(newunit = fID, file = trim(paramFile), action = "read", status = "old" )
         
         findIfStatsAreActive = .false.

         do
            read(fID,'(A)', iostat = io) line
            if ( io .ne. 0 ) exit
      
            if ( index(getSquashedLine(line), '#definestatistics') .ne. 0 ) then
               findIfStatsAreActive = .true.
               close(fID)
               return
            end if
            
         end do

         close(fID)
         
      end function findIfStatsAreActive

      subroutine rewriteControlFile(paramFile, no_of_lines)
         implicit none
         character(len=*), intent(in)  :: paramFile
         integer,          intent(in)  :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fID,i
         character(len=LINE_LENGTH) :: lines(no_of_lines)

         open(newunit = fID, file = trim(paramFIle), action = "read", status = "old" )

         
         do i = 1, no_of_lines

            read(fID,'(A)') lines(i)
            
            select case (trim(adjustl(lines(i))))
            
            case(START_COMMAND)
               lines(i) = "   " // START_COMMAND // "*"
            
            case(PAUSE_COMMAND)
               lines(i) = "   " // PAUSE_COMMAND // "*"

            case(STOP_COMMAND)
               lines(i) = "   " // STOP_COMMAND // "*"

            case(DUMP_COMMAND)
               lines(i) = "   " // DUMP_COMMAND // "*"

            case(RESET_COMMAND)
               lines(i) = "   " // RESET_COMMAND // "*"

            end select
         end do

         close(fID)
   
         open(newunit = fID, file = trim(paramFile), action = "write", status = "old" )

         do i = 1, no_of_lines
            write(fID,'(A)') trim(lines(i))
         end do

         close(fID)

      end subroutine rewriteControlFile

      subroutine StatisticsMonitor_WriteLabel ( self )
         implicit none
         class(StatisticsMonitor_t)             :: self

         if ( self % state .eq. OFF ) return
         write(STD_OUT , '(3X,A10)' , advance = "no") "Stats."

      end subroutine StatisticsMonitor_WriteLabel

      subroutine StatisticsMonitor_WriteValue ( self ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(StatisticsMonitor_t) :: self
         integer                 :: bufferLine
         character(len=10)       :: state

         if ( self % state .eq. OFF ) return

         select case ( self % state)
         case(ACTIVE_STATE)
            state = "Active"
   
         case(WAITING_STATE)
            state = "Waiting"

         case(STOP_STATE)
            state = "Inactive"   

         end select
   
         write(STD_OUT , '(1X,A,1X,A10)' , advance = "no") "|" , trim(state)

      end subroutine StatisticsMonitor_WriteValue 

end module StatisticsMonitor
#endif