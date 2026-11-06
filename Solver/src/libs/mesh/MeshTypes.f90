      Module MeshTypes
      use SMConstants
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      private
      public   HMESH_NONE, HMESH_UNDEFINED, HMESH_INTERIOR
      public   HMESH_BOUNDARY, HMESH_MPI, emptyBCName
      public   EFRONT, EBACK, EBOTTOM, ERIGHT, ETOP, ELEFT
      public   SPECMESH, HOPRMESH, GMSHMESH, MORTAR_NONE
      public   MORTAR_BIG, MORTAR_SMALL4, MORTAR_SLIDING
      public   small2big, big2small
   
      public indexesOnOtherFace, leftIndexes2Right, coordRotation

      integer, parameter :: EFRONT = 1, EBACK = 2, EBOTTOM = 3
      integer, parameter :: ERIGHT = 4, ETOP = 5, ELEFT = 6

   
      integer, parameter :: HMESH_NONE      = 0       ! Not constructed
      integer, parameter :: HMESH_UNDEFINED = -1      ! Constructed but undefined
      integer, parameter :: HMESH_INTERIOR  = 1       ! Interior face
      integer, parameter :: HMESH_BOUNDARY  = 2       ! Physical boundary face
      integer, parameter :: HMESH_MPI       = 3       ! MPI face

      integer, parameter :: MORTAR_NONE     = 0       ! No mortar, h-conforming face
      integer, parameter :: MORTAR_BIG      = 1       ! Big mortar (master) face 
      integer, parameter :: MORTAR_SMALL4   = 2       ! Small mortar (slave) face 
      integer, parameter :: MORTAR_SLIDING  = 3       ! Sliding mortar

      integer, parameter :: small2big       = 0
      integer, parameter :: big2small       = 1

      CHARACTER(LEN=3), parameter   :: emptyBCName = "---"
      
      integer, parameter :: SPECMESH = 1
      integer, parameter :: HOPRMESH = 2
      integer, parameter :: GMSHMESH = 3
      
      contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------
!  indexesOnOtherFace: 
!     This routine takes indexes on the a Face of a mortar and
!     outputs the corresponding indexes on the other Face 
!  -------------------------------------------------------------------------
   subroutine indexesOnOtherFace(i,j,Nx,Ny,rotation,side,ii,jj)
      implicit none
      !-arguments---------------------------------------------
      integer, intent(in)  :: i,j       !<  Input indexes
      integer, intent(in)  :: Nx, Ny    !<  Polynomial orders
      integer, intent(in)  :: rotation  !<  Face rotation
      integer, intent(in)  :: side  !<  Face rotation
      integer, intent(out) :: ii,jj     !>  Output indexes
      !-------------------------------------------------------
      
      select case (side)
         case (LEFT)
            call leftIndexes2Right(i,j,Nx,Ny,rotation,ii,jj)
         case (RIGHT)
            call rightIndexes2Left(i,j,Nx,Ny,rotation,ii,jj)
         case default 
            print *, "indexesOnOtherFace ERROR: Unknown side of face"
      end select
      
   end subroutine indexesOnOtherFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------
!  leftIndexes2Right: 
!     This routine takes indexes on the left Face of a mortar and
!     outputs the corresponding indexes on the right Face 
!  -------------------------------------------------------------------------
   subroutine leftIndexes2Right(i,j,Nx,Ny,rotation,ii,jj)
      implicit none
      !-arguments---------------------------------------------
      integer, intent(in)  :: i,j       !<  Input indexes
      integer, intent(in)  :: Nx, Ny    !<  Polynomial orders
      integer, intent(in)  :: rotation  !<  Face rotation
      integer, intent(out) :: ii,jj     !>  Output indexes
      !-------------------------------------------------------
      
      select case (rotation)
         case (0)
            ii = i
            jj = j
         case (1)
            ii = Ny - j
            jj = i
         case (2)
            ii = Nx - i
            jj = Ny - j
         case (3)
            ii = j
            jj = Nx - i
         case (4)
            ii = j
            jj = i
         case (5)
            ii = Nx - i
            jj = j
         case (6)
            ii = Ny - j
            jj = Nx - i
         case (7)
            ii = i
            jj = Ny - j
         case default 
            print *, "ERROR: Unknown rotation in element faces"
      end select
      
   end subroutine leftIndexes2Right
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------
!  rightIndexes2Left: 
!     This routine takes indexes on the tight Face of a mortar and
!     outputs the corresponding indexes on the left Face 
!  -------------------------------------------------------------------------
   subroutine rightIndexes2Left(i,j,Nx,Ny,rotation,ii,jj)
      implicit none
      !-arguments---------------------------------------------
      integer, intent(in)  :: i,j       !<  Input indexes
      integer, intent(in)  :: Nx, Ny    !<  Polynomial orders
      integer, intent(in)  :: rotation  !<  Face rotation
      integer, intent(out) :: ii,jj     !>  Output indexes
      !-------------------------------------------------------
      
      select case (rotation)
         case (0)
            ii = i
            jj = j
         case (1)
            jj = Nx - i
            ii = j
         case (2)
            ii = Nx - i
            jj = Ny - j
         case (3)
            jj = i
            ii = Ny - j
         case (4)
            jj = i
            ii = j
         case (5)
            ii = Nx - i
            jj = j
         case (6)
            jj = Nx - i
            ii = Ny - j
         case (7)
            ii = i
            jj = Ny - j
         case default 
            print *, "ERROR: Unknown rotation in element faces"
      end select
      
   end subroutine rightIndexes2Left
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  --------------------------------------------------------------------------
!  coordRotation: Takes the coordinates on the master face (left) and returns
!                 the corresponding coordinates on the slave face (right)
!  --------------------------------------------------------------------------
   SUBROUTINE coordRotation(xi,eta,rotation,xiRot, etaRot)
      IMPLICIT NONE
      real(kind=RP), intent(in)   :: xi, eta       ! Master coords 
      integer,       intent(in)   :: rotation      ! Face rotation
      real(kind=RP), intent(out)  :: xiRot, etaRot ! Slave coords
      
      SELECT CASE (rotation)
      CASE (0)
         xiRot  = xi
         etaRot = eta
      CASE (1)
         xiRot  = eta
         etaRot = -xi
      CASE (2)
         xiRot  = - xi
         etaRot = - eta
      CASE (3)
         xiRot  = -eta
         etaRot =  xi
      CASE (4)
         xiRot  = eta
         etaRot = xi
      CASE (5)
         xiRot  = - xi
         etaRot =  eta
      CASE (6)
         xiRot  = - eta
         etaRot = - xi
      CASE (7)
         xiRot  = - eta
         etaRot =  xi
      CASE DEFAULT 
         PRINT *, "ERROR: Unknown rotation in element faces"
      end SELECT
      
   end SUBROUTINE coordRotation



      END Module MeshTypes