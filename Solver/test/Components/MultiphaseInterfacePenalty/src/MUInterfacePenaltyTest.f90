#include "Includes.h"
!
!  Checks the BR1 interface penalty of the multiphase viscous fluxes (computeElementInterfaceFlux_MU).
!  The mobility flux of the concentration differentiates mu, so its penalty must act on [[mu]], not on [[c]].
!
!  Only ViscousDiscretization % sigma is changed between evaluations; the Cahn-Hilliard Laplacian uses
!  CHDiscretization, so mu and every other term of QDot stay the same, and the difference
!  D = QDot(sigma=1) - QDot(sigma=0) is exactly the penalty contribution. With the penalty on mu,
!     sum_nodes w J mu D_c = - sum_faces sum_ij w_i w_j Jf_ij sigma0 M0 (mu_L - mu_R)**2
!
program MUInterfacePenaltyTest
   use SMConstants
   use FTValueDictionaryClass
   use PhysicsStorage
   use SharedBCModule
   use DGSEMClass
   use MPI_Process_Info
   use SpatialDiscretization
   use pAdaptationClass,        only: GetMeshPolynomialOrders
   use NodalStorageClass
   use FluidData
   use FileReaders,             only: ReadControlFile
   use InterpolationMatrices,   only: Initialize_InterpolationMatrices
   use EllipticDiscretizations, only: ViscousDiscretization
   use MeshTypes,               only: HMESH_INTERIOR
   implicit none

   type QDot_t
      real(kind=RP), allocatable :: v(:,:,:,:)
   end type QDot_t

   type(FTValueDictionary)   :: controlVariables
   type(DGSem)               :: sem
   type(QDot_t), allocatable :: Q0(:), Q1(:), Q2(:)
   integer, allocatable      :: Nx(:), Ny(:), Nz(:)
   integer                   :: Nmax, nFailed, eID, fID, i, j, k, N
   logical                   :: success
   real(kind=RP)             :: x(NDIM), r, w, sigma0, dmu, dc
   real(kind=RP)             :: maxD_c, maxD_other, maxQDot, maxLin, sumD, sumAbsD
   real(kind=RP)             :: energy, expected_mu, expected_c
   real(kind=RP), allocatable :: D(:,:,:,:)

   call SetSolver(MULTIPHASE_SOLVER)
   call MPI_Process % Init
   call controlVariables % initWithSize(16)
   call ConstructSharedBCModule
   call ReadControlFile(controlVariables)
   call ConstructPhysicsStorage(controlVariables, success)
   if (.not. success) error stop "Physics parameters input error"
   call GetMeshPolynomialOrders(controlVariables, Nx, Ny, Nz, Nmax)
   call InitializeNodalStorage(controlVariables, Nmax)
   call Initialize_InterpolationMatrices(Nmax)
   call sem % construct(controlVariables = controlVariables, Nx_ = Nx, Ny_ = Ny, Nz_ = Nz, success = success)
   if (.not. success) error stop "Mesh construction error"
   call Initialize_SpaceAndTimeMethods(controlVariables, sem % mesh)

   do fID = 1, size(sem % mesh % faces)
      if (sem % mesh % faces(fID) % faceType /= HMESH_INTERIOR) error stop "The test needs a fully periodic mesh"
   end do
!
!  Smooth state plus a node-wise perturbation of c, so that [[c]] and [[mu]] are nonzero and not proportional
!  ---------------------------------------------------------------------------------------------------------
!
   do eID = 1, size(sem % mesh % elements)
      associate(e => sem % mesh % elements(eID))
      do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
         x = e % geom % x(:,i,j,k)
         r = modulo(sin(12.9898_RP*(eID + 7*i + 13*j + 29*k)) * 43758.5453_RP, 1.0_RP) - 0.5_RP
         e % storage % QNS(IMC,i,j,k)      = 0.5_RP + 0.3_RP*sin(PI*x(1))*cos(PI*x(2))*sin(PI*x(3)) + 0.04_RP*r
         e % storage % QNS(IMSQRHOU,i,j,k) = 0.2_RP*sin(PI*x(2))
         e % storage % QNS(IMSQRHOV,i,j,k) = 0.1_RP*cos(PI*x(3))
         e % storage % QNS(IMSQRHOW,i,j,k) = 0.1_RP*sin(PI*x(1))
         e % storage % QNS(IMP,i,j,k)      = 0.5_RP*cos(PI*x(1))*cos(PI*x(3))
      end do                ; end do                ; end do
      end associate
   end do

   call EvaluateWithPenalty(0.0_RP, Q0)
   call EvaluateWithPenalty(2.0_RP, Q2)
   call EvaluateWithPenalty(1.0_RP, Q1)   ! last, so that storage % mu and the face states belong to this state
!
!  Element sums
!  ------------
!
   maxD_c = 0.0_RP ; maxD_other = 0.0_RP ; maxQDot = 0.0_RP ; maxLin = 0.0_RP
   sumD = 0.0_RP ; sumAbsD = 0.0_RP ; energy = 0.0_RP
   do eID = 1, size(sem % mesh % elements)
!     Explicit source/mold: assigning the difference to an unallocated D would give it 1-based bounds
      if (allocated(D)) deallocate(D)
      allocate(D, mold = Q1(eID) % v)
      D = Q1(eID) % v - Q0(eID) % v
      associate(e => sem % mesh % elements(eID))
      maxD_c     = max(maxD_c,     maxval(abs(D(IMC,:,:,:))))
      maxD_other = max(maxD_other, maxval(abs(D(IMSQRHOU:IMP,:,:,:))))
      maxQDot    = max(maxQDot,    maxval(abs(Q0(eID) % v)))
      maxLin     = max(maxLin,     maxval(abs(Q2(eID) % v - Q0(eID) % v - 2.0_RP*D)))
      do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
         w = NodalStorage(e % Nxyz(1)) % w(i) * NodalStorage(e % Nxyz(2)) % w(j) &
           * NodalStorage(e % Nxyz(3)) % w(k) * e % geom % jacobian(i,j,k)
         sumD    = sumD    + w * D(IMC,i,j,k)
         sumAbsD = sumAbsD + w * abs(D(IMC,i,j,k))
         energy  = energy  + w * e % storage % mu(1,i,j,k) * D(IMC,i,j,k)
      end do                ; end do                ; end do
      end associate
   end do
!
!  Face sums: expected value with the penalty on mu, and the value it would have with the penalty on c
!  ---------------------------------------------------------------------------------------------------
!
   expected_mu = 0.0_RP ; expected_c = 0.0_RP
   do fID = 1, size(sem % mesh % faces)
      associate(f => sem % mesh % faces(fID))
      N = maxval(f % Nf)
      sigma0 = 0.5_RP * N * (N + 1) / f % geom % h
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
         w   = NodalStorage(f % Nf(1)) % w(i) * NodalStorage(f % Nf(2)) % w(j) * f % geom % jacobian(i,j)
         dmu = f % storage(1) % mu(1,i,j) - f % storage(2) % mu(1,i,j)
         dc  = f % storage(1) % Q(IMC,i,j) - f % storage(2) % Q(IMC,i,j)
         expected_mu = expected_mu - w * sigma0 * multiphase % M0 * dmu * dmu
         expected_c  = expected_c  - w * sigma0 * multiphase % M0 * dc  * dmu
      end do              ; end do
      end associate
   end do

   write(STD_OUT,'(/,A)')        "  Multiphase interface penalty test"
   write(STD_OUT,'(A,ES24.16)')  "     max |QDot|                      = ", maxQDot
   write(STD_OUT,'(A,ES24.16)')  "     max |D_c|                       = ", maxD_c
   write(STD_OUT,'(A,ES24.16)')  "     sum w J mu D_c                  = ", energy
   write(STD_OUT,'(A,ES24.16)')  "     expected, penalty on [[mu]]     = ", expected_mu
   write(STD_OUT,'(A,ES24.16,/)')"     value if penalty were on [[c]]  = ", expected_c

   nFailed = 0
   call Check("1. D is zero in momentum and pressure", maxD_other <= 1.0e-13_RP * maxQDot)
   call Check("2. QDot is affine in the penalty parameter", maxLin <= 1.0e-10_RP * maxD_c)
   call Check("3. The penalty conserves c", abs(sumD) <= 1.0e-12_RP * sumAbsD)
   call Check("4. Energy identity: the penalty acts on [[mu]]", &
              abs(energy - expected_mu) <= 1.0e-10_RP * abs(expected_mu))
   call Check("5. The penalty is active (D_c /= 0)", maxD_c > 1.0e-10_RP * maxQDot)

   if (nFailed > 0) then
      write(STD_OUT,'(/,I0,A)') nFailed, " check(s) failed"
      error stop 99
   end if
   write(STD_OUT,'(/,A)') "  All checks passed"

contains

   subroutine EvaluateWithPenalty(penalty, QDot)
      real(kind=RP),             intent(in)  :: penalty
      type(QDot_t), allocatable, intent(out) :: QDot(:)
      integer :: eID

      ViscousDiscretization % sigma = penalty
      call ComputeTimeDerivative(sem % mesh, sem % particles, 0.0_RP, CTD_IGNORE_MODE)
      allocate(QDot(size(sem % mesh % elements)))
      do eID = 1, size(sem % mesh % elements)
         QDot(eID) % v = sem % mesh % elements(eID) % storage % QDot
      end do
   end subroutine EvaluateWithPenalty

   subroutine Check(name, passed)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: passed

      if (passed) then
         write(STD_OUT,'(A,A,A)') "     ", name, ": passed"
      else
         write(STD_OUT,'(A,A,A)') "     ", name, ": FAILED"
         nFailed = nFailed + 1
      end if
   end subroutine Check

end program MUInterfacePenaltyTest
