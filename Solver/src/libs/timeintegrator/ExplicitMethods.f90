!
!////////////////////////////////////////////////////////////////////////
!
!      RK integrators for DG approximation to conservation
!      laws in 3D
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE ExplicitMethods
   USE SMConstants
   use HexMeshClass
   use TimeIntegratorDefinitions
   use DGSEMClass, only: ComputeTimeDerivative_f
   use ParticlesClass
   use PhysicsStorage, only: CTD_IGNORE_MODE
   use SlidingMeshProcedures
   IMPLICIT NONE

   private
   public   EULER_NAME, RK3_NAME, RK5_NAME, OPTRK_NAME, SSPRK33_NAME, SSPRK43_NAME, EULER_RK3_NAME, LSERK14_4_NAME, MIXED_RK_NAME, ML_RK3_NAME
   public   EULER_KEY, RK3_KEY, RK5_KEY, OPTRK_KEY, SSPRK33_KEY, SSPRK43_KEY, EULER_RK3_KEY, LSERK14_4_KEY, MIXED_RK_KEY, ML_RK3_KEY
   public   TakeExplicitEulerStep, TakeRK3Step, TakeRK5Step, TakeRKOptStep, TakeMLRK3Step
   public   TakeSSPRK33Step, TakeSSPRK43Step, TakeEulerRK3Step, TakeLSERK14_4Step, TakeMixedRKStep
   public   Enable_CTD_AFTER_STEPS, Enable_limiter, CTD_AFTER_STEPS, LIMITED, LIMITER_MIN

   integer,  protected :: eBDF_order = 3
   logical,  protected :: CTD_AFTER_STEPS = .false.
   logical,  protected :: LIMITED = .false.
   real(RP), protected :: LIMITER_MIN = 1e-13_RP
!
!  Implemented integration methods
!  -------------------------------
!  Remember to add them to `Enable_limiter` if they support limiting (or any kind of stage
!  callbacks if that is ever implemented)
!  ---------------------------------------------------------------------------------------
   character(len=*), parameter :: EULER_NAME   = "euler"
   character(len=*), parameter :: RK3_NAME     = "rk3"
   character(len=*), parameter :: RK5_NAME     = "rk5"
   character(len=*), parameter :: OPTRK_NAME   = "optimal rk"
   character(len=*), parameter :: SSPRK33_NAME = "ssprk33"
   character(len=*), parameter :: SSPRK43_NAME = "ssprk43"
   character(len=*), parameter :: EULER_RK3_NAME = "euler rk3"
   character(len=*), parameter :: LSERK14_4_NAME = "lserk14-4"
   character(len=*), parameter :: MIXED_RK_NAME = "mixedrk"
   character(len=*), parameter :: ML_RK3_NAME  = "multi level rk3"


   integer, parameter :: EULER_KEY   = 1
   integer, parameter :: RK3_KEY     = 2
   integer, parameter :: RK5_KEY     = 3
   integer, parameter :: OPTRK_KEY   = 4
   integer, parameter :: SSPRK33_KEY = 5
   integer, parameter :: SSPRK43_KEY = 6
   integer, parameter :: EULER_RK3_KEY = 7
   integer, parameter :: LSERK14_4_KEY = 8
   integer, parameter :: MIXED_RK_KEY  = 9
   integer, parameter :: ML_RK3_KEY    = 10

!========
 CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking an explicit Euler - RK3 step depending on the polynomial order of each element.
!  ------------------------------
 SUBROUTINE TakeEulerRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter, dtAdaptation)
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)   :: t, deltaT, tk
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)


      INTEGER :: i, j, k, id, eID
      integer :: interval

      interval = 10 !Compute time derivative each interval for those elements whose pmax = 1

      if (present(dt_vec)) then   
         
         do k = 1,3
            tk = t + b(k)*deltaT
            if ((k==1) .and. (mod(iter, interval) == 0)) then
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            else
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .true.)
            endif
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1) then
!$omp parallel do schedule(runtime) private(eID)
               do id = 1, SIZE( mesh % LO_Elements )
                  ! Explicit Euler
                  eID = mesh % LO_Elements(id)
#ifdef FLOW 
                  mesh % elements(eID) % storage % Q = mesh % elements(eID) % storage % Q + dt_vec(eID)*mesh % elements(eID) % storage % QDot
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(eID) % storage % c = mesh % elements(eID) % storage % c + dt_vec(eID)*mesh % elements(eID) % storage % cDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime) private(eID)
            do id = 1, SIZE( mesh % HO_Elements )
               ! Runge-Kutta 3
               eID = mesh % HO_Elements(id)
#ifdef FLOW
               mesh % elements(eID) % storage % G_NS = a(k)* mesh % elements(eID) % storage % G_NS  +              mesh % elements(eID) % storage % QDot
               mesh % elements(eID) % storage % Q =       mesh % elements(eID) % storage % Q  + c(k)*dt_vec(eID)* mesh % elements(eID) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
               mesh % elements(eID) % storage % G_CH = a(k)*mesh % elements(eID) % storage % G_CH + mesh % elements(eID) % storage % cDot
               mesh % elements(eID) % storage % c    = mesh % elements(eID) % storage % c         + c(k)*dt_vec(eID)* mesh % elements(eID) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do

         end do ! k

      else !not present dt_vec

         do k = 1,3
            tk = t + b(k)*deltaT
            if ((k==1) .and. (mod(iter, interval) == 0)) then
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            else
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .true.)
            endif
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1) then
!$omp parallel do schedule(runtime) private(eID)
               do id = 1, SIZE( mesh % LO_Elements )
                  ! Explicit Euler
                  eID = mesh % LO_Elements(id)
#ifdef FLOW 
                  mesh % elements(eID) % storage % Q = mesh % elements(eID) % storage % Q + deltaT*mesh % elements(eID) % storage % QDot
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(eID) % storage % c = mesh % elements(eID) % storage % c + deltaT*mesh % elements(eID) % storage % cDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime) private(eID)
            do id = 1, SIZE( mesh % HO_Elements )
               ! Runge-Kutta 3
               eID = mesh % HO_Elements(id)
#ifdef FLOW
               mesh % elements(eID) % storage % G_NS = a(k)* mesh % elements(eID) % storage % G_NS  +              mesh % elements(eID) % storage % QDot
               mesh % elements(eID) % storage % Q =       mesh % elements(eID) % storage % Q  + c(k)*deltaT* mesh % elements(eID) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
               mesh % elements(eID) % storage % G_CH = a(k)*mesh % elements(eID) % storage % G_CH + mesh % elements(eID) % storage % cDot
               mesh % elements(eID) % storage % c    = mesh % elements(eID) % storage % c         + c(k)*deltaT* mesh % elements(eID) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do

         end do ! k

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeEulerRK3Step


SUBROUTINE TakeMixedRKStep( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter, dtAdaptation)
!
!        *****************************************************************************************
!          Uses RK3 in phase1 and LSERK14-4 in phase2, think phase1 is air and phase2 is water. Of course, this will only yield an advantage if: 
!          - you can get away with a time step using the MixedRK that is 5x times larger than that in pure RK3 scheme,
!          - AND if there are at LEAST 20+% air elements (aim for 50% or more)
!          A stand-alone LSERK14-4 implimetnation is given in TakeLSERK14_4Step
!          The MixedRK scheme is valid only for the multiphase solver and throws an error otherwise,  
!          However, the rest of this function still has its implementation consistent
!          with the other RK steppers in case anyone needed to extend it
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      logical, allocatable :: phase1_mask(:)
      logical, allocatable :: phase2_mask(:)
      integer                   :: id, i, j, k
      integer                   :: phase1_count, phase2_count, total_elements
      real(kind=RP)             :: phase1_percentage


#if !defined(MULTIPHASE)
   error stop 'MixedRK time stepping is only implemented for the multiphase solver'
#else

      ! Arrays to store the original solution for all elements
#ifdef FLOW
      real(kind=RP), allocatable :: Q_orig(:,:,:,:,:)
      real(kind=RP), allocatable :: G_NS_orig(:,:,:,:,:)
      ! Arrays to store the updated solution for phase1 elements
      real(kind=RP), allocatable :: Q_phase1(:,:,:,:,:)
      real(kind=RP), allocatable :: G_NS_phase1(:,:,:,:,:)
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
      real(kind=RP), allocatable :: c_orig(:,:,:,:,:)
      real(kind=RP), allocatable :: G_CH_orig(:,:,:,:,:)
      ! Arrays to store the updated solution for phase1 elements
      real(kind=RP), allocatable :: c_phase1(:,:,:,:,:)
      real(kind=RP), allocatable :: G_CH_phase1(:,:,:,:,:)
#endif

      integer, parameter        :: N_STAGES_RK3 = 3
      real(kind=RP), parameter  :: a_RK3(N_STAGES_RK3) = [0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP]
      real(kind=RP), parameter  :: b_RK3(N_STAGES_RK3) = [0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  ]
      real(kind=RP), parameter  :: c_RK3(N_STAGES_RK3) = [1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP ]

      integer, parameter        :: N_STAGES_RK14_4 = 14
      real(kind=RP), parameter  :: a_RK14_4(N_STAGES_RK14_4) = [0.0000000000000000_RP, -0.7188012108672410_RP, -0.7785331173421570_RP, -0.0053282796654044_RP, -0.8552979934029281_RP, -3.9564138245774565_RP, -1.5780575380587385_RP, -2.0837094552574054_RP, -0.7483334182761610_RP, -0.7032861106563359_RP, +0.0013917096117681_RP, -0.0932075369637460_RP, -0.9514200470875948_RP, -7.1151571693922548_RP]
      real(kind=RP), parameter  :: b_RK14_4(N_STAGES_RK14_4) = [0.0000000000000000_RP, 0.0367762454319673_RP, 0.1249685262725025_RP, 0.2446177702277698_RP, 0.2476149531070420_RP, 0.2969311120382472_RP, 0.3978149645802642_RP, 0.5270854589440328_RP, 0.6981269994175695_RP, 0.8190890835352128_RP, 0.8527059887098624_RP, 0.8604711817462826_RP, 0.8627060376969976_RP, 0.8734213127600976_RP]
      real(kind=RP), parameter  :: c_RK14_4(N_STAGES_RK14_4) = [0.0367762454319673_RP, 0.3136296607553959_RP, 0.1531848691869027_RP, 0.0030097086818182_RP, 0.3326293790646110_RP, 0.2440251405350864_RP, 0.3718879239592277_RP, 0.6204126221582444_RP, 0.1524043173028741_RP, 0.0760894927419266_RP, 0.0077604214040978_RP, 0.0024647284755382_RP, 0.0780348340049386_RP, 5.5059777270269628_RP]

      logical :: updateQLowRK

      if (present(dtAdaptation)) then
         updateQLowRK = dtAdaptation
      else
         updateQLowRK = .false.
      end if
      
      allocate(phase1_mask(size(mesh % elements)))
      allocate(phase2_mask(size(mesh % elements)))

      ! Initialize counters
      phase1_count = 0
      phase2_count = 0
      total_elements = size(mesh % elements)

!$omp parallel do schedule(runtime) reduction(+:phase1_count,phase2_count) private(id)
      do id = 1, size(mesh % elements)
         if ( all(mesh % elements(id) % storage % Q(1,:,:,:) > 1.0_RP - 1e-8_RP) ) then
             mesh % elements(id) % MLevel = 1 !required for adaptive time step
             phase1_mask(id) = .true.
             phase2_mask(id) = .false.
             phase1_count = phase1_count + 1
         else
            mesh % elements(id) % MLevel = 4 !required for adaptive time step
            phase1_mask(id) = .false.
            phase2_mask(id) = .true.
            phase2_count = phase2_count + 1
         endif
     end do
!$omp end parallel do

    ! Print phase distribution every 100 iterations for debugging
   !   if (present(iter)) then
   !    if (mod(iter, 100) == 0) then
   !       phase1_percentage = real(phase1_count, RP) / real(total_elements, RP) * 100.0_RP
   !       write(*,*) "Inside explicit method MixedRK"
   !       write(*,'(A,I8,A,F6.2,A,I8,A,I8,A)') 'Iteration ', iter, &
   !             ': Phase1 elements = ', phase1_percentage, '% (', &
   !             phase1_count, ' out of ', total_elements, ' elements)'
   !    end if
   ! end if


      ! Save a copy of the original solution for all elements
#ifdef FLOW
      allocate(Q_orig(size(mesh % elements), size(mesh % elements(1) % storage % Q, 1), &
                      size(mesh % elements(1) % storage % Q, 2), &
                      size(mesh % elements(1) % storage % Q, 3), &
                      size(mesh % elements(1) % storage % Q, 4)))
      allocate(G_NS_orig(size(mesh % elements), size(mesh % elements(1) % storage % G_NS, 1), &
                         size(mesh % elements(1) % storage % G_NS, 2), &
                         size(mesh % elements(1) % storage % G_NS, 3), &
                         size(mesh % elements(1) % storage % G_NS, 4)))
      
      do id = 1, size(mesh % elements)
         Q_orig(id,:,:,:,:) = mesh % elements(id) % storage % Q
         G_NS_orig(id,:,:,:,:) = mesh % elements(id) % storage % G_NS
      end do
      
      ! Allocate arrays for the updated solution for phase1 elements
      allocate(Q_phase1(size(mesh % elements), size(mesh % elements(1) % storage % Q, 1), &
                        size(mesh % elements(1) % storage % Q, 2), &
                        size(mesh % elements(1) % storage % Q, 3), &
                        size(mesh % elements(1) % storage % Q, 4)))
      allocate(G_NS_phase1(size(mesh % elements), size(mesh % elements(1) % storage % G_NS, 1), &
                           size(mesh % elements(1) % storage % G_NS, 2), &
                           size(mesh % elements(1) % storage % G_NS, 3), &
                           size(mesh % elements(1) % storage % G_NS, 4)))
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
      allocate(c_orig(size(mesh % elements), size(mesh % elements(1) % storage % c, 1), &
                      size(mesh % elements(1) % storage % c, 2), &
                      size(mesh % elements(1) % storage % c, 3), &
                      size(mesh % elements(1) % storage % c, 4)))
      allocate(G_CH_orig(size(mesh % elements), size(mesh % elements(1) % storage % G_CH, 1), &
                         size(mesh % elements(1) % storage % G_CH, 2), &
                         size(mesh % elements(1) % storage % G_CH, 3), &
                         size(mesh % elements(1) % storage % G_CH, 4)))
      
      do id = 1, size(mesh % elements)
         c_orig(id,:,:,:,:) = mesh % elements(id) % storage % c
         G_CH_orig(id,:,:,:,:) = mesh % elements(id) % storage % G_CH
      end do
      
      ! Allocate arrays for the updated solution for phase1 elements
      allocate(c_phase1(size(mesh % elements), size(mesh % elements(1) % storage % c, 1), &
                        size(mesh % elements(1) % storage % c, 2), &
                        size(mesh % elements(1) % storage % c, 3), &
                        size(mesh % elements(1) % storage % c, 4)))
      allocate(G_CH_phase1(size(mesh % elements), size(mesh % elements(1) % storage % G_CH, 1), &
                           size(mesh % elements(1) % storage % G_CH, 2), &
                           size(mesh % elements(1) % storage % G_CH, 3), &
                           size(mesh % elements(1) % storage % G_CH, 4)))
#endif

      ! Apply RK3 to phase1 elements and RK14-4 to phase2 elements
      if (present(dt_vec)) then 
         ! Apply RK3 to phase1 elements
         DO k = 1, N_STAGES_RK3

            tk = t + b_RK3(k)*deltaT
            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .false., phase1_mask)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then !For adaptive time step only, update QLowRK
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  if (phase1_mask(id)) then
#ifdef FLOW 
                     mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
#endif
                  end if
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
            DO id = 1, SIZE( mesh % elements )
               if (phase1_mask(id)) then
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a_RK3(k)* mesh % elements(id) % storage % G_NS + mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + c_RK3(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a_RK3(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c = mesh % elements(id) % storage % c + c_RK3(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
               end if
            END DO
!$omp end parallel do
         END DO

         ! Save the updated solution for phase1 elements
#ifdef FLOW
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               Q_phase1(id,:,:,:,:) = mesh % elements(id) % storage % Q
               G_NS_phase1(id,:,:,:,:) = mesh % elements(id) % storage % G_NS
            end if
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               c_phase1(id,:,:,:,:) = mesh % elements(id) % storage % c
               G_CH_phase1(id,:,:,:,:) = mesh % elements(id) % storage % G_CH
            end if
         end do
#endif

         ! Restore original solution for all elements
#ifdef FLOW
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q = Q_orig(id,:,:,:,:)
            mesh % elements(id) % storage % G_NS = G_NS_orig(id,:,:,:,:)
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % c = c_orig(id,:,:,:,:)
            mesh % elements(id) % storage % G_CH = G_CH_orig(id,:,:,:,:)
         end do
#endif

         ! Apply RK14-4 to phase2 elements
         DO k = 1, N_STAGES_RK14_4
            tk = t + b_RK14_4(k)*deltaT
            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .false., phase2_mask)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then !For adaptive time step only, update QLowRK
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  if (phase2_mask(id)) then
#ifdef FLOW 
                     mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
#endif
                  end if
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
            DO id = 1, SIZE( mesh % elements )
               if (phase2_mask(id)) then
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a_RK14_4(k)* mesh % elements(id) % storage % G_NS + mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + c_RK14_4(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a_RK14_4(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c = mesh % elements(id) % storage % c + c_RK14_4(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
               end if
            END DO
!$omp end parallel do
         END DO

         ! Combine the updated solutions
#ifdef FLOW
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               mesh % elements(id) % storage % Q = Q_phase1(id,:,:,:,:)
               mesh % elements(id) % storage % G_NS = G_NS_phase1(id,:,:,:,:)
            end if
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               mesh % elements(id) % storage % c = c_phase1(id,:,:,:,:)
               mesh % elements(id) % storage % G_CH = G_CH_phase1(id,:,:,:,:)
            end if
         end do
#endif

      else
         ! Same procedure but with deltaT instead of dt_vec
         ! Apply RK3 to phase1 elements
         DO k = 1, N_STAGES_RK3
            tk = t + b_RK3(k)*deltaT
            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .false., phase1_mask)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then !For adaptive time step only, update QLowRK
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  if (phase1_mask(id)) then
#ifdef FLOW 
                     mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
#endif
                  end if
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
            DO id = 1, SIZE( mesh % elements )
               if (phase1_mask(id)) then
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a_RK3(k)* mesh % elements(id) % storage % G_NS + mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + c_RK3(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a_RK3(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c = mesh % elements(id) % storage % c + c_RK3(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
               end if
            END DO
!$omp end parallel do
         END DO

         ! Save the updated solution for phase1 elements
#ifdef FLOW
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               Q_phase1(id,:,:,:,:) = mesh % elements(id) % storage % Q
               G_NS_phase1(id,:,:,:,:) = mesh % elements(id) % storage % G_NS
            end if
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               c_phase1(id,:,:,:,:) = mesh % elements(id) % storage % c
               G_CH_phase1(id,:,:,:,:) = mesh % elements(id) % storage % G_CH
            end if
         end do
#endif

         ! Restore original solution for all elements
#ifdef FLOW
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q = Q_orig(id,:,:,:,:)
            mesh % elements(id) % storage % G_NS = G_NS_orig(id,:,:,:,:)
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % c = c_orig(id,:,:,:,:)
            mesh % elements(id) % storage % G_CH = G_CH_orig(id,:,:,:,:)
         end do
#endif

         ! Apply RK14-4 to phase2 elements
         DO k = 1, N_STAGES_RK14_4
            tk = t + b_RK14_4(k)*deltaT
            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .false., phase2_mask)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then !For adaptive time step only, update QLowRK
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
                  if (phase2_mask(id)) then
#ifdef FLOW 
                     mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
#endif
                  end if
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
            DO id = 1, SIZE( mesh % elements )
               if (phase2_mask(id)) then
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a_RK14_4(k)* mesh % elements(id) % storage % G_NS + mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + c_RK14_4(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a_RK14_4(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c = mesh % elements(id) % storage % c + c_RK14_4(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
               end if
            END DO
!$omp end parallel do
         END DO

         ! Combine the updated solutions
#ifdef FLOW
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               mesh % elements(id) % storage % Q = Q_phase1(id,:,:,:,:)
               mesh % elements(id) % storage % G_NS = G_NS_phase1(id,:,:,:,:)
            end if
         end do
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
         do id = 1, size(mesh % elements)
            if (phase1_mask(id)) then
               mesh % elements(id) % storage % c = c_phase1(id,:,:,:,:)
               mesh % elements(id) % storage % G_CH = G_CH_phase1(id,:,:,:,:)
            end if
         end do
#endif
      end if

      ! Deallocate temporary arrays
      deallocate(phase1_mask)
      deallocate(phase2_mask)
#ifdef FLOW
      deallocate(Q_orig)
      deallocate(G_NS_orig)
      deallocate(Q_phase1)
      deallocate(G_NS_phase1)
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
      deallocate(c_orig)
      deallocate(G_CH_orig)
      deallocate(c_phase1)
      deallocate(G_CH_phase1)
#endif

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

! end if flow is multiphase
#endif 

   end subroutine TakeMixedRKStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a RK3 step.
!  ------------------------------
   SUBROUTINE TakeRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter, dtAdaptation)
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
            type(Particles_t)  :: particles
#else
            logical            :: particles
#endif
      REAL(KIND=RP)   :: t, deltaT, tk
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
      REAL(KIND=RP), DIMENSION(3) :: auxiliar = (/1.0_RP/3.0_RP, (3.0_RP/4.0_RP-1.0_RP/3.0_RP) , 1.0_RP/4.0_RP /)
      real(kind=RP) :: center(2)
      real(kind=RP), DIMENSION(3) :: omega = (/0.0_RP, PI/20_RP, 0.0_RP/)
      real(kind=RP) :: rad, angle
      integer :: NB = 7


      INTEGER :: i, j, k, id
      logical :: updateQLowRK

      if (present(dtAdaptation)) then
         updateQLowRK = dtAdaptation
      else
         updateQLowRK = .false.
      end if

      if (present(dt_vec)) then   
         
         do k = 1,3
            tk = t + b(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            !if (mesh % slidingMesh % active .and. mesh % slidingflux) then
            if (mesh % slidingMesh % active) then
               center(1)=0.0_RP
               center(2)=0.0_RP
               rad=1.01_RP
               angle=(omega(2)*deltaT) * auxiliar(k)
               !CALL AdvanceSlidingMesh(mesh=mesh, rotationRadius=rad, rotationCenter=center, numBFacePoints=NB, nodes=mesh%nodeType, useMPI=.FALSE. , angle = angle)
            end if
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
#ifdef FLOW 
                  mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
                  do id = 1, SIZE( mesh % elements )
#ifdef FLOW
                        mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
                        mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif
       
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                        mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                        mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
                  end do ! id
!$omp end parallel do

         end do ! k

      else

         do k = 1,3
            tk = t + b(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
                  !if (mesh % slidingMesh % active .and. mesh % slidingflux) then
            if (mesh % slidingMesh % active) then
               center(1)=0.0_RP
               center(2)=0.0_RP
               rad=1.01_RP
               angle=(omega(2)*deltaT) * auxiliar(k)
               !CALL AdvanceSlidingMesh(mesh=mesh, rotationRadius=rad, rotationCenter=center, numBFacePoints=NB, nodes=mesh%nodeType, useMPI=.FALSE. , angle = angle)
            end if
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1 .and. updateQLowRK) then
!$omp parallel do schedule(runtime)
               do id = 1, SIZE( mesh % elements )
#ifdef FLOW 
                  mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime)
                  do id = 1, SIZE( mesh % elements )
#ifdef FLOW
                        mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
                        mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif
       
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                        mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                        mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
                  end do ! id
!$omp end parallel do
       
               end do ! k
       
            end if
      !
      !     To obtain the updated residuals
            if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)
       
            call checkForNan(mesh, t)
       
         END SUBROUTINE TakeRK3Step

   SUBROUTINE TakeRK5Step( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter, dtAdaptation)
!
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k
      integer, parameter         :: N_STAGES = 5
      real(kind=RP), parameter  :: a(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP ,     -1.697784692471_RP , -1.514183444257_RP ]
      real(kind=RP), parameter  :: b(N_STAGES) = [0.0_RP , 0.1496590219993_RP , 0.3704009573644_RP , 0.6222557631345_RP , 0.9582821306748_RP ]
      real(kind=RP), parameter  :: c(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]


      if (present(dt_vec)) then 

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      else

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeRK5Step

   SUBROUTINE TakeLSERK14_4Step( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter, dtAdaptation)
!
!        *****************************************************************************************
!           This 14 stage LSRK scheme has been extracted from the paper: 
!           "Efﬁcient low-storage Runge–Kutta schemes with optimized stability regions"
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                   :: id, i, j, k
      integer, parameter        :: N_STAGES = 14
      real(kind=RP), parameter  :: a(N_STAGES) = [0.0000000000000000_RP, -0.7188012108672410_RP, -0.7785331173421570_RP, -0.0053282796654044_RP, -0.8552979934029281_RP, -3.9564138245774565_RP, -1.5780575380587385_RP, -2.0837094552574054_RP, -0.7483334182761610_RP, -0.7032861106563359_RP, +0.0013917096117681_RP, -0.0932075369637460_RP, -0.9514200470875948_RP, -7.1151571693922548_RP]
      real(kind=RP), parameter  :: b(N_STAGES) = [0.0000000000000000_RP, 0.0367762454319673_RP, 0.1249685262725025_RP, 0.2446177702277698_RP, 0.2476149531070420_RP, 0.2969311120382472_RP, 0.3978149645802642_RP, 0.5270854589440328_RP, 0.6981269994175695_RP, 0.8190890835352128_RP, 0.8527059887098624_RP, 0.8604711817462826_RP, 0.8627060376969976_RP, 0.8734213127600976_RP]
      real(kind=RP), parameter  :: c(N_STAGES) = [0.0367762454319673_RP, 0.3136296607553959_RP, 0.1531848691869027_RP, 0.0030097086818182_RP, 0.3326293790646110_RP, 0.2440251405350864_RP, 0.3718879239592277_RP, 0.6204126221582444_RP, 0.1524043173028741_RP, 0.0760894927419266_RP, 0.0077604214040978_RP, 0.0024647284755382_RP, 0.0780348340049386_RP, 5.5059777270269628_RP]


      if (present(dt_vec)) then 

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      else

      DO k = 1, N_STAGES

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(k)*deltaT* mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(k)*deltaT* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeLSERK14_4Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a SSP-RK 3-stage 3rd-order step.
!  ------------------------------
   subroutine TakeSSPRK33Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter, dtAdaptation)
      implicit none
!
!     ----------------
!     Input parameters
!     ----------------
!
      type(HexMesh)                               :: mesh
#if defined(FLOW)
      type(Particles_t)                           :: particles
#else
      logical                                     :: particles
#endif
      real(RP)                                    :: t
      real(RP)                                    :: deltaT
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
      real(RP), allocatable, optional, intent(in) :: dt_vec(:)
      logical,               optional, intent(in) :: dts
      real(RP),              optional, intent(in) :: global_dt
      integer,               optional, intent(in) :: iter
      logical,               intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      real(RP), parameter :: a(3) = [1.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP]
      real(RP), parameter :: b(3) = [0.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP]
      real(RP), parameter :: c(3) = [1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP]
      real(RP), parameter :: d(3) = [0.0_RP, 1.0_RP,        0.5_RP]
      real(RP) :: tk
      integer  :: i, j, k, id


      if (present(dt_vec)) then

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do k = 1, 3
            tk = t + d(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(k) * mesh % elements(id) % storage % G_NS &
                                                 + b(k) * mesh % elements(id) % storage % Q    &
                                                 + c(k) * dt_vec(id) * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(k) * mesh % elements(id) % storage % G_CH &
                                                 + b(k) * mesh % elements(id) % storage % c    &
                                                 + c(k) * dt_vec(id) * mesh % elements(id) % storage % cDot
#endif
            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! k

      else

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do k = 1, 3
            tk = t + d(k) * deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(k) * mesh % elements(id) % storage % G_NS &
                                                 + b(k) * mesh % elements(id) % storage % Q    &
                                                 + c(k) * deltaT * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(k) * mesh % elements(id) % storage % G_CH &
                                                 + b(k) * mesh % elements(id) % storage % c    &
                                                 + c(k) * deltaT * mesh % elements(id) % storage % cDot
#endif
            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! k

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeSSPRK33Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a SSP-RK 4-stage 3rd-order step.
!  ------------------------------
   subroutine TakeSSPRK43Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter, dtAdaptation)
      implicit none
!
!     ----------------
!     Input parameters
!     ----------------
!
      type(HexMesh)                               :: mesh
#if defined(FLOW)
      type(Particles_t)                           :: particles
#else
      logical                                     :: particles
#endif
      real(RP)                                    :: t
      real(RP)                                    :: deltaT
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
      real(RP), allocatable, optional, intent(in) :: dt_vec(:)
      logical,               optional, intent(in) :: dts
      real(RP),              optional, intent(in) :: global_dt
      integer,               optional, intent(in) :: iter
      logical,               intent(in), optional :: dtAdaptation
      
!
!     ---------------
!     Local variables
!     ---------------
!
      real(RP), parameter :: a(4) = [1.0_RP, 0.0_RP, 2.0_RP/3.0_RP, 0.0_RP]
      real(RP), parameter :: b(4) = [0.0_RP, 1.0_RP, 1.0_RP/3.0_RP, 1.0_RP]
      real(RP), parameter :: c(4) = [0.5_RP, 0.5_RP, 1.0_RP/6.0_RP, 0.5_RP]
      real(RP), parameter :: d(4) = [0.0_RP, 0.5_RP, 1.0_RP,        0.5_RP]
      real(RP) :: tk
      integer  :: i, j, k, id


      if (present(dt_vec)) then

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do k = 1, 4
            tk = t + d(k)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(k) * mesh % elements(id) % storage % G_NS &
                                                 + b(k) * mesh % elements(id) % storage % Q    &
                                                 + c(k) * dt_vec(id) * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(k) * mesh % elements(id) % storage % G_CH &
                                                 + b(k) * mesh % elements(id) % storage % c    &
                                                 + c(k) * dt_vec(id) * mesh % elements(id) % storage % cDot
#endif
            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! k

      else

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do k = 1, 4
            tk = t + d(k) * deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(k) * mesh % elements(id) % storage % G_NS &
                                                 + b(k) * mesh % elements(id) % storage % Q    &
                                                 + c(k) * deltaT * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(k) * mesh % elements(id) % storage % G_CH &
                                                 + b(k) * mesh % elements(id) % storage % c    &
                                                 + c(k) * deltaT * mesh % elements(id) % storage % cDot
#endif

            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! k

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeSSPRK43Step

   SUBROUTINE TakeExplicitEulerStep( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter, dtAdaptation)
!
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, k

      CALL ComputeTimeDerivative( mesh, particles, t, CTD_IGNORE_MODE)
      if ( present(dts) ) then
         if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
      end if

      if (present(dt_vec)) then
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do
      else
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
         END DO
!$omp end parallel do
      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeExplicitEulerStep

   subroutine TakeExplicitBDFStep(mesh, particles, t, deltaT, ComputeTimeDerivative)
      implicit none
      type(HexMesh)                      :: mesh
#ifdef FLOW
      type(Particles_t)                  :: particles
#else
      logical                            :: particles
#endif
      REAL(KIND=RP)                      :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                  :: id
      real(kind=RP), parameter :: invGamma1 = 1.0_RP, invGamma2 = 2.0_RP/3.0_RP, invGamma3 = 6.0_RP / 11.0_RP
      logical, save            :: isFirst = .true., isSecond = .false., isThird = .false.

      if (isThird) then
!
!        Perform the third order stages
!        ------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
!           Set y^{*,n+1} in Q and downgrade y^n and y^{n-1}
!           ------------------------------------------------
            mesh % elements(id) % storage % QDot = mesh % elements(id) % storage % prevQ(2) % Q
            mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q
            mesh % elements(id) % storage % Q =   3.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                  - 3.0_RP * mesh % elements(id) % storage % prevQ(2) % Q &
                                                  + mesh % elements(id) % storage % QDot
         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q =   2.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                - 0.5_RP * mesh % elements(id) % storage % prevQ(2) % Q &
                                                + (1.0_RP/3.0_RP) * mesh % elements(id) % storage % Q &
                                                + deltaT * mesh % elements(id) % storage % QDot
            mesh % elements(id) % storage % Q = invGamma3 * mesh % elements(id) % storage % Q
         end do
!$omp end parallel do


      elseif (isSecond) then
!
!        Perform the second order stages
!        -------------------------------
         if (eBDF_ORDER > 2) then
!$omp parallel do schedule(runtime)
            do id = 1, size(mesh % elements)
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            end do
!$omp end parallel do
         end if

!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
!
!           Set y^{*,n+1} in Q
!           ------------------
            mesh % elements(id) % storage % Q = 2.0_RP*mesh % elements(id) % storage % Q -mesh % elements(id) % storage % prevQ(1) % Q
!
!           Set y^{n} in prevQ
!           --------------------
            mesh % elements(id) % storage % prevQ(1) % Q = 0.5_RP*(mesh % elements(id) % storage % Q + mesh % elements(id) % storage % prevQ(1) % Q)

         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q =   mesh % elements(id) % storage % prevQ(1) % Q + 0.5_RP*mesh % elements(id) % storage % Q  &
                                                + deltaT * mesh % elements(id) % storage % QDot
            mesh % elements(id) % storage % Q = invGamma2 * mesh % elements(id) % storage % Q

         end do
!$omp end parallel do

         if (eBDF_ORDER > 1 ) then
!
!           Move to third order
!           -------------------
            isFirst = .false.
            isSecond = .false.
            isThird = .true.
         end if

      elseif (isFirst) then
!
!        Perform the first order stages
!        ------------------------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            if (eBDF_ORDER > 1 ) then
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q

            end if

            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot

         END DO
!$omp end parallel do

            if (eBDF_ORDER > 1 ) then
!
!              Move to second order
!              --------------------
               isFirst = .false.
               isSecond = .true.
            end if

      end if




   end subroutine TakeExplicitBDFStep

   subroutine Enable_CTD_AFTER_STEPS()
      implicit none
      CTD_AFTER_STEPS = .true.
   end subroutine Enable_CTD_AFTER_STEPS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE TakeRKOptStep( mesh, particles, t, deltaT, ComputeTimeDerivative , N_STAGES, dt_vec, dts, global_dt, iter, dtAdaptation)
!
!        *****************************************************************************************
!       Optimal RK coefficients from Bassi2009
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      integer, intent(in)         :: N_STAGES
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k
      real(kind=RP), dimension(6,7) :: Am, Bm
      real(kind=RP) :: a(N_STAGES), b(N_STAGES)

      Bm(1,:) = (/ 0.9998596842476678, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Bm(2,:) = (/ 0.528003175664866, 0.5193233361621609, 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Bm(3,:) = (/ 0.4062766523561424, 0.3590274668186006, 0.2782786562184366, 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Bm(4,:) = (/ 0.32451232607547, 0.2850381110111294, 0.2299189950459751, 0.3245118697485674, 0.1925289659886354, 0.0, 0.0 /)
      Bm(5,:) = (/ 0.2712469842000987, 0.2506886071464794, 0.1571621623659122, 0.2281484031198761, 0.2523511867383585, 0.1918765315676819, 0.0 /)
      Bm(6,:) = (/ 0.2328811281838825, 0.2007437175473038, 0.1577268986576998, 0.2052094755795549, 0.2138853585222901, 0.192146045128098, 0.1428487872093191 /)

      Am(1,:) = (/ 0.0 / 1.0, -0.4998596842476678 / 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Am(2,:) = (/ 0.0 / 1.0, -0.1645541151603977 / 0.315637725010225, -0.20368561115193584 / 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Am(3,:) = (/ 0.0 / 1.0, -0.11076800268308828 / 0.1937832814991212, -0.1652441853194794 / 0.207607995049003, &
            -0.0706706611694336 / 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Am(4,:) = (/ 0.0 / 1.0, -0.11517541944711737 / 0.1896693024156952, -0.0953688085954342 / 0.2159361297115306, &
            -0.01398286533444451 / 0.1925286952557861, -0.1319831744927813 / 0.1925289659886354, 0.0, 0.0 /)
      Am(5,:) = (/ 0.0 / 1.0, -0.056511008554826186 / 0.1229547684000589, -0.1277338387464205 / 0.1094339401033201, &
            -0.047728222262592115 / 0.1824888900957738, -0.045659513024102316 / 0.1785098941878928,-0.0738412925504657 / 0.1918765315676819, 0.0 /)
      Am(6,:) = (/ 0.0 / 1.0, -0.044038551801784204 / 0.1267393084745009, -0.0740044090728029 / 0.1061061571001407, &
            -0.051620741557559094 / 0.1525740388611983, -0.052635436718356604 / 0.1650271578073516,-0.04885820071493849 / 0.1178619741653911, &
             -0.0742840709627069 / 0.1428487872093191 /)

      a = Am(N_STAGES-1,1:N_STAGES)
      b = Bm(N_STAGES-1,1:N_STAGES)

      tk = t + deltaT

      if (present(dt_vec)) then 

      DO k = 1, N_STAGES

         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  +  dt_vec(id) * mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q   +  b(k) * mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH  +  dt_vec(id) * mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c          +  b(k) * mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      else

      DO k = 1, N_STAGES

         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(k)* mesh % elements(id) % storage % G_NS  + deltaT * mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + b(k) * mesh % elements(id) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(k)*mesh % elements(id) % storage % G_CH + deltaT * mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + b(k) * mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeRKOptStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for Multi Level RK3
!  ------------------------------
   SUBROUTINE TakeMLRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter, dtAdaptation)
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)      :: t, deltaT
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
      logical, intent(in), optional :: dtAdaptation
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(4) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  ,     1.0_RP     /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
	  REAL(KIND=RP), DIMENSION(3) :: d = (/1.0_RP/3.0_RP,  5.0_RP/12.0_RP,    1.0_RP/4.0_RP /)        ! Temporal step of stages
	  


      INTEGER       :: i, j, id, locLevel, k2, k3, k1, nLevel
	  INTEGER, ALLOCATABLE :: k(:)
	  INTEGER, parameter :: kmax = 3
	  REAL(KIND=RP) :: corrector
	  REAL(KIND=RP), ALLOCATABLE :: cL(:) , tk(:), deltaTL(:) ! 

      logical :: updateQLowRK

      if (present(dtAdaptation)) then
         updateQLowRK = dtAdaptation
      else
         updateQLowRK = .false.
      end if
	  
	  nLevel = mesh % MLRK % maxLevel
	  allocate(k(nLevel), cL(nLevel), tk(nLevel), deltaTL(nLevel))
	  
	  deltaTL(:) = deltaT		! Local timestep, for level 1 is the deltaT
	  tk(:)      = t
	  
	  associate ( MLIter_eID => mesh % MLRK % MLIter_eID, MLIter => mesh % MLRK % MLIter  )
	  
      k(:) = kmax                ! A counter array to track at which stage and level , 3 is the total number of stages
      do k1 = 1,kmax
!           LEVEL 1
            tk(:) = t + b(k1)*deltaT    											! Actual time at each level 
            call ComputeTimeDerivative( mesh, particles, tk(1), CTD_IGNORE_MODE)  	! Compute Qdot all elements
			locLevel = 1															! Level locator

         if (k1==1 .and. updateQLowRK) then !For adaptive time step only, update QLowRK
!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
#ifdef FLOW 
               mesh % elements(id) % storage % QLowRK = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot
#endif
            end do ! id
!$omp end parallel do
         end if

!           -------------------------------------------------------------------------------------------------------------------------------
!           LEVEL 2-LEVEL N-1
!           -------------------------------------------------------------------------------------------------------------------------------
			do k2 = 1, max(kmax**(nLevel-2),1)				! loop for the intermediate levels nstages**(nLevel-2)
				k(nLevel-1) = k(nLevel-1)+1					! Update the counter level-stages by adding by 1 for nLevel-1
				do i=nLevel-1,1,-1							! Check if stages already exceed total number of stages
					if (k(i).gt.kmax) then
						k(i)=1
						if (i.ne.1) then
							k(i-1)=k(i-1) +1
						end if 
					else
						exit
					end if						
				end do 
				
				do i=2, nLevel 								! Update the local timestep for level 2: nLevel
					deltaTL(i) = deltaTL(i-1) * d(k(i-1))
				end do 
!               -------------------------------------------------------------------------------------------------------------------------------
!               LEVEL N 
!               -------------------------------------------------------------------------------------------------------------------------------
				do k3 = 1,kmax
					k(nLevel) = k3							! Update the counter level-stages for nLevel 
!           		Update G_NS/G_CH from QDot - Depend on level etc
!$omp parallel do schedule(runtime)	private(id, corrector)	
					do i = 1, MLIter(locLevel,1)            ! Loop all elements at the active level to update the G / Qdot contribution
						id = MLIter_eID(i)					! ID of true element at partition
	  
						if (locLevel.eq.nLevel) then		! Variable a, weight of the previous Qdot contribution
							corrector = a(k(nLevel))		
					    elseif(i.gt.MLIter(min(locLevel+1,nLevel),1)) then 
							corrector = a(k(locLevel))		! If locLevel is j then all element with level > j will be on first stage corrector=0
						else
							corrector = 0.0_RP
						end if 
						
						associate(storage => mesh % elements(id) % storage)
#ifdef FLOW        
						storage % G_NS = corrector * storage % G_NS  +  storage % QDot
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
						storage % G_CH = corrector * storage % G_CH + storage % cDot
#endif
						end associate
					end do 
!$omp end parallel do	
				
!           		Marching in time, all elements
					corrector = 1.0_RP
					do i = 1, nLevel
						corrector = corrector * d(k(i))
					end do 
                    do i = 1, nLevel
						cL(i) = c(k(i)) * corrector/d(k(i)) ! Timestep scaling for local timesteping at level N for all level 
					end do 

!$omp parallel do schedule(runtime)	private(id, corrector)			
					do i = 1, MLIter(1,1)
						do j = nLevel, 1, -1
						    corrector = cL(j)
							if (i.le.MLIter(j,1)) exit
						end do
						
						id = MLIter_eID(i)
						
						associate(storage => mesh % elements(id) % storage)
#ifdef FLOW             
						storage % Q    = storage % Q  + corrector * deltaT * storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
						storage % c    = storage % c + corrector * deltaT * storage % G_CH
#endif
						end associate
					end do
!$omp end parallel do	
                    if (all(k(2:nLevel) == kmax)) then
						exit
					else
						do i=nLevel,2,-1
						    locLevel =i                   	! The level for next loop
							if (k(i).ne.kmax) exit
						end do 
					end if 
					
					tk(locLevel)=tk(locLevel-1)+b(k(locLevel)+1) * deltaTL(locLevel)
					call ComputeTimeDerivative( mesh, particles, tk(locLevel), CTD_IGNORE_MODE, Level=locLevel) ! Update Qdot
				end do 
			end do 

       end do ! k1
	   
	   end associate
	   
	   call ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE) ! Necessary for residual computation
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)
	  
	  deallocate(k, cL, tk, deltaTL)

   END SUBROUTINE TakeMLRK3Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Enable_limiter(integrator, minimum)
!
!     ---------
!     Interface
!     ---------
      integer,            intent(in) :: integrator
      real(RP), optional, intent(in) :: minimum


      LIMITED = (integrator == SSPRK33_KEY) .or. (integrator == SSPRK43_KEY)

      if (present(minimum)) then
         LIMITER_MIN = minimum
      end if

   end subroutine Enable_limiter

#if defined(NAVIERSTOKES) || defined(SPALARTALMARAS)
   subroutine stage_limiter(mesh)
!
!     -------
!     Modules
!     -------
!
      use ElementClass,      only: Element
      use NodalStorageClass, only: NodalStorage
      use FluidData,         only: thermodynamics
!
!     ---------
!     Interface
!     ---------
!
      type(HexMesh), target, intent(inout) :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      real(RP)               :: m
      real(RP)               :: Q(5), Qavg(5)
      real(RP)               :: rho, minrho
      real(RP)               :: p, pavg, minp
      real(RP)               :: theta
      real(RP)               :: gm1
      integer                :: eID
      integer                :: i, j, k
      type(Element), pointer :: e
      real(RP),      pointer :: wx(:), wy(:), wz(:)


      gm1 = thermodynamics % gammaMinus1

!$omp parallel do default(private) shared(mesh, NodalStorage) firstprivate(gm1, LIMITER_MIN)
      do eID = 1, mesh % no_of_elements
         e  => mesh % elements(eID)
         wx => NodalStorage(e % Nxyz(1)) % w
         wy => NodalStorage(e % Nxyz(2)) % w
         wz => NodalStorage(e % Nxyz(3)) % w

         ! Compute averages
         Qavg = 0.0_RP
         do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)
            Qavg = Qavg + e % storage % Q(:,i,j,k) * wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k)
         end do               ; end do               ; end do
         Qavg = Qavg / e % geom % volume

         ! Density first
         minrho = huge(1.0_RP)
         do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)
            rho = e % storage % Q(1,i,j,k)
            if (rho < minrho) minrho = rho
         end do               ; end do               ; end do

         if (Qavg(1) /= minrho) then
            m = min(LIMITER_MIN, Qavg(1))
            theta = abs((Qavg(1) - m) / (Qavg(1) - minrho))
            if (theta <= 1.0_RP) then
               e % storage % Q(1,:,:,:) = theta * (e % storage % Q(1,:,:,:) - Qavg(1)) + Qavg(1)
            end if
         end if

         ! Pressure now (Jensen's inequality is NOT conservative for the pressure)
         minp = huge(1.0_RP)
         pavg = 0.0_RP
         do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)
            Q = e % storage % Q(:,i,j,k)
            p = gm1 * (Q(5) - 0.5_RP * (Q(2)**2 + Q(3)**2 + Q(4)**2) / Q(1))
            pavg = pavg + p * wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k)
            if (p < minp) minp = p
         end do               ; end do               ; end do
         pavg = pavg / e % geom % volume

         if (pavg /= minp) then
            m = min(LIMITER_MIN, pavg)
            theta = abs((pavg - m) / (pavg - minp))
            if (theta <= 1.0_RP) then
               do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)
                  e % storage % Q(:,i,j,k) = theta * (e % storage % Q(:,i,j,k) - Qavg) + Qavg
               end do               ; end do               ; end do
            end if
         end if

      end do
!$omp end parallel do

      nullify(e)
      nullify(wx)
      nullify(wy)
      nullify(wz)

   end subroutine stage_limiter
#else
   subroutine stage_limiter(mesh)
      type(HexMesh), intent(in) :: mesh
   end subroutine stage_limiter
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   Subroutine checkForNan(mesh, t)
!
!        **************************************************************************************************************
!        Look if there is a nan in the solution, if at least one if found, stops and creates a hsol for user inspection
!        **************************************************************************************************************
!
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      type(HexMesh), intent(in)     :: mesh
      real(kind=RP), intent(in)     :: t

      !local variables
      integer                    :: eID
      CHARACTER(len=LINE_LENGTH) :: FinalName      !  Final name for particular restart file
      logical                    :: NanNotFound, allNan
      integer                    :: ierr

      ! use not found instead of found as OMP reduction initialized the private value as true
      ! this is redundant for the OMP but left for non parallel compilations
      NanNotFound = .TRUE.

!$omp parallel do reduction(.AND.:NanNotFound) schedule(runtime)
      do eID=1, mesh % no_of_elements
         if ( any(isnan(mesh % elements(eID) % storage % Q))) then
            NanNotFound = .FALSE.
         endif
      end do
!$omp end parallel do

#ifdef _HAS_MPI_
      call mpi_allreduce(NanNotFound, allNan, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      NanNotFound = allNan
#endif

      if (.not. NanNotFound) then
          if ( MPI_Process % isRoot ) print*, "Numerical divergence obtained in solver."
          WRITE(FinalName,'(A,ES11.5,A)')  'RESULTS/horses_divergence_',t,'.hsol'
          if ( MPI_Process % isRoot ) print *, "Writing failed solution: ", FinalName
          call mesh % SaveSolution(0, t, FinalName, .FALSE.)
      end if

#ifdef _HAS_MPI_
     call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

      if (.not. NanNotFound) call exit(99)


   End Subroutine checkForNan
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods