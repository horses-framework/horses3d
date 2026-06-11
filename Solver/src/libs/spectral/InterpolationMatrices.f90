#include "Includes.h"
module InterpolationMatrices
   use PolynomialInterpAndDerivsModule
   use NodalStorageClass
   use SMConstants
   implicit none

   private
   public InterpolationMatrix_t
   public Tset, TsetM
   public Interp3DArrays, Interp3DArraysOneDir
   public Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
!
!  ---------------------------------
!  Interpolation matrix derived type
!  ---------------------------------
!
   type InterpolationMatrix_t
      logical                    :: Constructed = .false.
      real(kind=RP), allocatable :: T(:,:)
      contains
         procedure :: construct => InterpolationMatrix_Construct
         procedure :: destruct  => InterpolationMatrix_Destruct
   end type InterpolationMatrix_t

!  **************************************************
!  The set of interpolation matrices is stored here.
!  >> The following criteria is adopted:
!        * Tset(N,M) with N=M is the identity matrix.
!        * Tset(N,M) with N<M is a forward matrix.
!        * Tset(N,M) with N>M is a backwards matrix.
!  **************************************************
   type(InterpolationMatrix_t),   allocatable, target :: Tset(:,:)
   type(InterpolationMatrix_t),   allocatable, target :: TsetM(:,:,:,:)

   interface Interp3DArrays
      module procedure Interp3DArrays_Tset, Interp3DArrays_interp
   end interface Interp3DArrays
!========
 contains
!========

   subroutine InterpolationMatrix_Construct(this, Norigin, Ndest, offset_, scale_)
!
!     ****************************************************************
!        This subroutine computes the interpolation matrix
!        from Norigin to Ndest
!     ****************************************************************
!
      implicit none
      class(InterpolationMatrix_t), intent(inout) :: this
      integer                     , intent(in)    :: Norigin !<  Origin polynomial order
      integer                     , intent(in)    :: Ndest   !<  Destination polynomial order
      real(kind=RP)      ,optional, intent(in)    :: offset_ ! Offset for h-mortars or sliding meshes
      real(kind=RP)      ,optional, intent(in)    :: scale_  ! Scale for h-mortars or sliding meshes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                       :: i, j
      real(kind=RP), allocatable    :: Ttemp(:,:)
      type(NodalStorage_t), pointer :: spAo
      type(NodalStorage_t), pointer :: spAd
      real(kind=RP)                 :: offset, scale

!
!     Evaluate if it's necessary to construct the matrix
!     **************************************************

      if ( this % Constructed ) return

      !
      ! Check if offset and scale are present (used for h-mortars and sliding mortars)
      !
      if (present(offset_)) then
         offset = offset_
      else
         offset = 0.0_rp
      end if
      if (present(scale_)) then
         scale = scale_
      else
         scale = 1.0_rp
      end if

!
!     Matrix construction
!     *******************

      spAo => NodalStorage(Norigin)
      spAd => NodalStorage(Ndest)

!
!     Allocate memory
!     ---------------
      allocate( this % T(0:Ndest, 0:Norigin) )
!
!     Construct the forward interpolation matrix
!     ------------------------------------------
      if (Norigin < Ndest) then
!
!        Prolongation matrix
!        -------------------
            call PolynomialInterpolationMatrix(Norigin, Ndest, spAo % x, spAo % wb, offset + scale * spAd % x, this % T)
         
      else
!
!        Restriction (backwards) matrix
!        ------------------------------
         allocate( Ttemp(0:Norigin, 0:Ndest) )

         call PolynomialInterpolationMatrix(Ndest, Norigin, spAd % x, spAd % wb, offset + scale * spAo % x, Ttemp)

         this % T = transpose(Ttemp)

         do j = 0, Norigin ; do i = 0, Ndest
            this % T(i,j) = this % T(i,j) * spAo % w(j) / spAd % w(i)
         end do            ; end do
      end if
!
!     Set constructed flag
!     ********************
      this % Constructed = .TRUE.

      nullify (spAo)
      nullify (spAd)

   end subroutine InterpolationMatrix_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------
!  Class destructor
!  ----------------

   elemental subroutine InterpolationMatrix_Destruct(this)
      implicit none
      class(InterpolationMatrix_t), intent(inout) :: this

      safedeallocate (this % T)
      this % Constructed = .FALSE.
   end subroutine InterpolationMatrix_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for interpolating an array of 3D arrays
!  --------------------------------------------------
   pure subroutine Interp3DArrays_Tset(Nvars, Nin, inArray, Nout, outArray )
      implicit none
      !-------------------------------------------------------
      integer                                                        , intent(in)  :: Nvars
      integer      , dimension(3)                                    , intent(in)  :: Nin
      integer      , dimension(3)                                    , intent(in)  :: Nout
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)  :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out) :: outArray
      !-------------------------------------------------------
      integer :: i,j,k,l,m,n
      !-------------------------------------------------------

      outArray = 0.0_RP

      do n = 0, Nin(3)  ; do k = 0, Nout(3)
         do m = 0, Nin(2)  ; do j = 0, Nout(2)
            do l = 0, Nin(1)  ; do i = 0, Nout(1)
               outArray(:,i,j,k) = outArray(:,i,j,k) +   Tset(Nin(1),Nout(1)) % T (i,l) &
                                                       * Tset(Nin(2),Nout(2)) % T (j,m) &
                                                       * Tset(Nin(3),Nout(3)) % T (k,n) &
                                                       * inArray(:,l,m,n)
            end do             ; end do
         end do             ; end do
      end do             ; end do

   end subroutine Interp3DArrays_Tset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for interpolating an array of 3D arrays
!  --------------------------------------------------
   pure subroutine Interp3DArrays_interp(Nvars, Nin, inArray, Nout, outArray, interpXi, interpEta, interpZeta )
      implicit none
      !-------------------------------------------------------
      integer                                                        , intent(in)  :: Nvars
      integer      , dimension(3)                                    , intent(in)  :: Nin
      integer      , dimension(3)                                    , intent(in)  :: Nout
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)  :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out) :: outArray
      real(kind=RP), dimension(0:Nout(1), 0:Nin(1)), target          , intent(in)  :: interpXi
      real(kind=RP), dimension(0:Nout(2), 0:Nin(2)), target          , intent(in)  :: interpEta
      real(kind=RP), dimension(0:Nout(3), 0:Nin(3)), target          , intent(in)  :: interpZeta
      !-------------------------------------------------------
      integer :: i,j,k,l,m,n
      !-------------------------------------------------------

      outArray = 0.0_RP

      do n = 0, Nin(3)  ; do k = 0, Nout(3)
         do m = 0, Nin(2)  ; do j = 0, Nout(2)
            do l = 0, Nin(1)  ; do i = 0, Nout(1)
               outArray(:,i,j,k) = outArray(:,i,j,k) +   interpXi   (i,l) &
                                                       * interpEta  (j,m) &
                                                       * interpZeta (k,n) &
                                                       * inArray(:,l,m,n)
            end do             ; end do
         end do             ; end do
      end do             ; end do

   end subroutine Interp3DArrays_interp
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for interpolating an array of 3D arrays
!  --------------------------------------------------
   subroutine Interp3DArraysOneDir(Nvars, Nin, inArray, Nout, outArray, Dir )
      implicit none
      !-------------------------------------------------------
      integer                                                        , intent(in)   :: Nvars
      integer      , dimension(3)                                    , intent(in)   :: Nin
      integer      , dimension(3)                                    , intent(in)   :: Nout
      integer                                                        , intent(in)   :: Dir
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)   :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out)  :: outArray
      !-------------------------------------------------------
      integer :: i,j,k,l
      !-------------------------------------------------------

      outArray = 0.0_RP

      select case(Dir)
         case (1)

            do k = 0, Nout(3)
               do j = 0, Nout(2)
                  do l = 0, Nin(1)  ; do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + Tset(Nin(Dir),Nout(Dir)) % T(i,l) * inArray(:,l,j,k)
                  end do             ; end do
               end do
            end do

         case (2)

            do k = 0, Nout(3)
               do l = 0, Nin(2)  ; do j = 0, Nout(2)
                  do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + Tset(Nin(Dir),Nout(Dir)) % T(j,l) * inArray(:,i,l,k)
                  end do
               end do             ; end do
            end do

         case (3)

            do l = 0, Nin(3)  ; do k = 0, Nout(3)
               do j = 0, Nout(2)
                  do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + Tset(Nin(Dir),Nout(Dir)) % T(k,l) * inArray(:,i,j,l)
                  end do
               end do
            end do             ; end do

      end select

   end subroutine Interp3DArraysOneDir
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Initialize_InterpolationMatrices(Nmax)
      implicit none
      integer, intent(in) :: Nmax

      allocate ( Tset(0:Nmax,0:Nmax) )
      allocate ( TsetM(0:Nmax,0:Nmax, 4, 2) )
   end subroutine Initialize_InterpolationMatrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Finalize_InterpolationMatrices
      implicit none

      call Tset % destruct
      call TsetM % destruct
      deallocate ( Tset )
      deallocate ( TsetM )
   end subroutine Finalize_InterpolationMatrices
end module InterpolationMatrices