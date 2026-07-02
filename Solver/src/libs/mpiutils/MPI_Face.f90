#include "Includes.h"
module MPI_Face_Class
   use SMConstants
   use MPI_Process_Info
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   private
   public  MPI_FacesSet_t

   public  ConstructMPIFaces, DestructMPIFaces
   public  ConstructMPIFacesStorage

   type MPI_Face_t
      integer                    :: nDOFs
      integer                    :: no_of_faces
      integer, allocatable       :: faceIDs(:)
      integer, allocatable       :: elementSide(:)
      integer                    :: Nrecv_req
      integer                    :: Qrecv_req
      integer                    :: gradQrecv_req
      integer                    :: AviscFluxRecv_req
      integer                    :: QBaseRecv_req
      integer                    :: transportVRecv_req
      integer                    :: transportDRecv_req
      integer                    :: sizeQ
      integer                    :: sizeU_xyz
      integer                    :: sizeAviscFlux
      integer      , allocatable :: Nsend(:)          ! Information to send: [fNxi, fNeta, eNxi, eNeta, eNzeta, eGlobID]
      integer      , allocatable :: Nrecv(:)
      real(kind=RP), allocatable :: Qsend(:)
      real(kind=RP), allocatable :: U_xyzsend(:)
      real(kind=RP), allocatable :: Qrecv(:)
      real(kind=RP), allocatable :: U_xyzrecv(:)
      real(kind=RP), allocatable :: AviscFluxSend(:)
      real(kind=RP), allocatable :: AviscFluxRecv(:)
      real(kind=RP), allocatable :: QbaseSend(:)
      real(kind=RP), allocatable :: QbaseRecv(:)
      real(kind=RP), allocatable :: transportVSend(:)
      real(kind=RP), allocatable :: transportVRecv(:)
      real(kind=RP), allocatable :: transportDSend(:)
      real(kind=RP), allocatable :: transportDRecv(:)
      contains
         procedure   :: Construct        => MPI_Face_Construct
         procedure   :: Destruct         => MPI_Face_Destruct
         procedure   :: SendN            => MPI_Face_SendN
         procedure   :: RecvN            => MPI_Face_RecvN
         procedure   :: SendQ            => MPI_Face_SendQ
         procedure   :: RecvQ            => MPI_Face_RecvQ
         procedure   :: SendU_xyz        => MPI_Face_SendU_xyz
         procedure   :: RecvU_xyz        => MPI_Face_RecvU_xyz
         procedure   :: SendAviscFlux    => MPI_Face_SendAviscFlux
         procedure   :: RecvAviscFlux    => MPI_Face_RecvAviscFlux
         procedure   :: WaitForN         => MPI_Face_WaitForN
         procedure   :: WaitForSolution  => MPI_Face_WaitForSolution
         procedure   :: WaitForGradients => MPI_Face_WaitForGradients
         procedure   :: WaitForAviscFlux => MPI_Face_WaitForAviscFlux
         procedure   :: SendQbase        => MPI_Face_SendQbase
         procedure   :: RecvQbase        => MPI_Face_RecvQbase
         procedure   :: WaitForSolutionB => MPI_Face_WaitForSolutionB
         procedure   :: SendTransportV   => MPI_Face_SendTransportV
         procedure   :: RecvTransportV   => MPI_Face_RecvTransportV
         procedure   :: WaitForTransportV => MPI_Face_WaitForTransportV
         procedure   :: SendTransportD   => MPI_Face_SendTransportD
         procedure   :: RecvTransportD   => MPI_Face_RecvTransportD
         procedure   :: WaitForTransportD => MPI_Face_WaitForTransportD
   end type MPI_Face_t

   type MPI_FacesSet_t
      logical                          :: constructed = .false.
      integer                          :: nDomainShared
      type(MPI_Face_t), allocatable    :: faces(:)
      integer         , allocatable    :: listDomain(:)
   end type MPI_FacesSet_t

   interface MPI_Face_t
      module procedure Construct_MPI_Face
   end interface MPI_Face_t

   contains

      subroutine ConstructMPIFaces(facesSet)
         implicit none
         type(MPI_FacesSet_t)    :: facesSet
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain

         if ( MPI_Process % doMPIAction ) then
            allocate(facesSet % faces(MPI_Process % nProcs))
            allocate(facesSet % listDomain(MPI_Process % nProcs))
            
            facesSet % nDomainShared = MPI_Process % nProcs
            
            do domain = 1, MPI_Process % nProcs
               facesSet % faces(domain) = MPI_Face_t()
               facesSet % faces(domain) % no_of_faces = 0
               facesSet % listDomain (domain) = domain
            end do
         end if

         facesSet % constructed = .TRUE.

      end subroutine ConstructMPIFaces

      subroutine ConstructMPIFacesStorage(facesSet, NCONS, NGRAD, NDOFS, NCONSB_in, transportData)
!
!        ***************************************************
!           Allocates buffers to send and receive.
!           This subroutine is called from:
!               mesh % SetElementConnectivitiesAndLinkFaces
!        ***************************************************
!
         implicit none
         type(MPI_FacesSet_t)    :: facesSet
         integer, intent(in)     :: NCONS, NGRAD
         integer, intent(in)     :: NDOFS(MPI_Process % nProcs)
         integer, intent(in), optional     :: NCONSB_in
         logical, intent(in), optional     :: transportData
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain, NCONSB

         if (present(NCONSB_in)) then
             NCONSB = NCONSB_in
         else
             NCONSB = 0
         end if 

         do domain = 1, MPI_Process % nProcs
            facesSet % faces(domain) % nDOFs         = NDOFS(domain)
            facesSet % faces(domain) % sizeQ         = NCONS * NDOFS(domain)
            facesSet % faces(domain) % sizeU_xyz     = NDIM * NGRAD * NDOFS(domain)
            facesSet % faces(domain) % sizeAviscFlux = NCONS * NDOFS(domain)

            if ( NDOFS(domain) .gt. 0 ) then

               safedeallocate(facesSet % faces(domain) % Qsend)
               safedeallocate(facesSet % faces(domain) % U_xyzsend)
               safedeallocate(facesSet % faces(domain) % Qrecv)
               safedeallocate(facesSet % faces(domain) % U_xyzrecv)
               safedeallocate(facesSet % faces(domain) % AviscFluxSend)
               safedeallocate(facesSet % faces(domain) % AviscFluxRecv)
               safedeallocate(facesSet % faces(domain) % QBaseSend)
               safedeallocate(facesSet % faces(domain) % QBaseRecv)

               allocate( facesSet % faces(domain) % Qsend(NCONS * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % U_xyzsend(NDIM * NGRAD * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % Qrecv(NCONS * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % U_xyzrecv(NDIM * NGRAD * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % AviscFluxSend(NCONS * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % AviscFluxRecv(NCONS * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % QBaseSend(NCONSB * NDOFS(domain)) )
               allocate( facesSet % faces(domain) % QBaseRecv(NCONSB * NDOFS(domain)) )
               if (transportData) then
                  safedeallocate(facesSet % faces(domain) % transportVSend)
                  safedeallocate(facesSet % faces(domain) % transportVRecv)
                  safedeallocate(facesSet % faces(domain) % transportDSend)
                  safedeallocate(facesSet % faces(domain) % transportDRecv)
                  allocate( facesSet % faces(domain) % transportVSend(NDIM * NDOFS(domain)) )
                  allocate( facesSet % faces(domain) % transportVRecv(NDIM * NDOFS(domain)) )
                  allocate( facesSet % faces(domain) % transportDSend(NDIM * NDIM * NDOFS(domain)) )
                  allocate( facesSet % faces(domain) % transportDRecv(NDIM * NDIM * NDOFS(domain)) )
               end if
            end if
         end do
      end subroutine ConstructMPIFacesStorage
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------
!     Subroutine to send the polynomial order
!     ---------------------------------------
      subroutine MPI_Face_SendN(self, domain, req)
         implicit none
         !-arguments-----------------------------------------------
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % Nsend, 6 * self % no_of_faces, MPI_INT, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendN
!     ------------------------------------------
!     Subroutine to receive the polynomial order
!     ------------------------------------------
      subroutine MPI_Face_RecvN(self, domain, req)
         implicit none
         !-arguments-----------------------------------------------
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % Nrecv, 6 * self % no_of_faces, MPI_INT, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, req, ierr)
            self % Nrecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvN
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------
!     Subroutine to send the variable Q
!     ---------------------------------
      subroutine MPI_Face_SendQ(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % Qsend, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendQ
!     ------------------------------------
!     Subroutine to receive the variable Q
!     ------------------------------------
      subroutine MPI_Face_RecvQ(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         self % Qrecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % Qrecv, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, req, ierr)
            self % Qrecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------
!     Subroutine to send the variable gradQ
!     -------------------------------------
      subroutine MPI_Face_SendU_xyz(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
		 integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
		 req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % U_xyzsend, nEqn * NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendU_xyz
!     ----------------------------------------
!     Subroutine to receive the variable gradQ
!     ----------------------------------------
      subroutine MPI_Face_RecvU_xyz(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
		 integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
		 req = MPI_REQUEST_NULL
		 self % gradQrecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % U_xyzrecv, nEqn * NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, req, ierr)
			self % gradQrecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvU_xyz
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------
!     Subroutine to send the AviscFlux
!     --------------------------------
      subroutine MPI_Face_SendAviscFlux(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
		 integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
		 req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % AviscFluxSend, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, req, ierr)
         end if
#endif

      end subroutine MPI_Face_SendAviscFlux
!     -----------------------------------
!     Subroutine to receive the AviscFlux
!     -----------------------------------
      subroutine MPI_Face_RecvAviscFlux(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
		 integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
		 req = MPI_REQUEST_NULL
		 self % AviscFluxRecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % AviscFluxRecv, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, &
                           MPI_ANY_TAG, MPI_COMM_WORLD, req, ierr)
			self % AviscFluxRecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvAviscFlux
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------
!     Subroutine to send the variable Qbase
!     ---------------------------------
      subroutine MPI_Face_SendQBase(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % QBaseSend, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendQBase
!     ------------------------------------
!     Subroutine to receive the variable Qbase
!     ------------------------------------
      subroutine MPI_Face_RecvQBase(self, domain, nEqn, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(in)    :: nEqn
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
        req = MPI_REQUEST_NULL
        self % QBaseRecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % QBaseRecv, nEqn * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, req, ierr)
            self % QBaseRecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvQBase
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------
!     Subroutine to send the variable transportV
!     ---------------------------------
      subroutine MPI_Face_SendTransportV(self, domain, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % transportVSend, NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendTransportV
!     ------------------------------------
!     Subroutine to receive the variable transportV
!     ------------------------------------
      subroutine MPI_Face_RecvTransportV(self, domain, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
        req = MPI_REQUEST_NULL
        self % transportVRecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % transportVRecv, NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, req, ierr)
            self % transportVRecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvTransportV
!     ---------------------------------
!     Subroutine to send the variable transportD
!     ---------------------------------
      subroutine MPI_Face_SendTransportD(self, domain, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
         req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % transportDSend, NDIM * NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, req, ierr)
         end if
#endif
      end subroutine MPI_Face_SendTransportD
!     ------------------------------------
!     Subroutine to receive the variable transportD
!     ------------------------------------
      subroutine MPI_Face_RecvTransportD(self, domain, req)
         implicit none
         class(MPI_Face_t)      :: self
         integer, intent(in)    :: domain
         integer, intent(out)   :: req
         !-local-variables-----------------------------------------
         integer  :: ierr
         !---------------------------------------------------------
#ifdef _HAS_MPI_
        req = MPI_REQUEST_NULL
        self % transportDRecv_req = MPI_REQUEST_NULL
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % transportDRecv, NDIM * NDIM * self % nDOFs, MPI_DOUBLE_PRECISION, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, req, ierr)
            self % transportDRecv_req = req
         end if
#endif
      end subroutine MPI_Face_RecvTransportD
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------
!     Wait until the polynomial order is received
!     -------------------------------------------
      subroutine MPI_Face_WaitForN(self)
         implicit none
         !-arguments-----------------------------------------------
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
         !-local-variables-----------------------------------------
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)
         !---------------------------------------------------------
         call mpi_wait(self % Nrecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForN
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForSolution(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)

!
!        -----------------------------------
!        Wait until the solution is received
!        -----------------------------------
!
         call mpi_wait(self % Qrecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForGradients(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr, status(MPI_STATUS_SIZE)

!
!        ---------------------------------
!        Wait until gradients are received
!        --------------------------------
!
         call mpi_wait(self % gradQrecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForGradients
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForAviscFlux(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)

!
!        -----------------------------------
!        Wait until the solution is received
!        -----------------------------------
!
         call mpi_wait(self % AviscFluxRecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForAviscFlux
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForSolutionB(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)

!
!        -----------------------------------
!        Wait until the solution is received
!        -----------------------------------
!
         call mpi_wait(self % QBaseRecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForSolutionB
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForTransportV(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)

!
!        -----------------------------------
!        Wait until the solution is received
!        -----------------------------------
!
         call mpi_wait(self % transportVRecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForTransportV
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine MPI_Face_WaitForTransportD(self)
         implicit none
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr
         integer              :: status(MPI_STATUS_SIZE)

!
!        -----------------------------------
!        Wait until the solution is received
!        -----------------------------------
!
         call mpi_wait(self % transportDRecv_req, status, ierr)
#endif

      end subroutine MPI_Face_WaitForTransportD
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DestructMPIFaces(facesSet)
         implicit none
         type(MPI_FacesSet_t)    :: facesSet
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain

         if ( MPI_Process % doMPIAction .and. (allocated(facesSet % faces))) then
            do domain = 1, MPI_Process % nProcs
               call facesSet % faces(domain) % Destruct()
            end do
            safedeallocate(facesSet % faces)
            safedeallocate(facesSet % listDomain)
         end if
         facesSet % constructed = .FALSE.

      end subroutine DestructMPIFaces

      type(MPI_Face_t) function Construct_MPI_Face()
         implicit none

         call MPI_Face_Destruct(Construct_MPI_Face)

      end function Construct_MPI_Face

      subroutine MPI_Face_Construct(self, no_of_faces)
         implicit none
         class(MPI_Face_t) :: self
         integer, intent(in)  :: no_of_faces

         self % no_of_faces = no_of_faces
         allocate(self % faceIDs(no_of_faces))
         allocate(self % elementSide(no_of_faces))

         allocate( self % Nsend(6 * no_of_faces) )
         allocate( self % Nrecv(6 * no_of_faces) )

         self % faceIDs       = -1
         self % elementSide   = -1
#ifdef _HAS_MPI_
         self % Nrecv_req         = MPI_REQUEST_NULL
         self % Qrecv_req         = MPI_REQUEST_NULL
         self % gradQrecv_req     = MPI_REQUEST_NULL
         self % AviscFluxRecv_req = MPI_REQUEST_NULL
         self % QBaseRecv_req     = MPI_REQUEST_NULL
         self % transportVRecv_req = MPI_REQUEST_NULL
         self % transportDRecv_req = MPI_REQUEST_NULL
#endif

      end subroutine MPI_Face_Construct

      subroutine MPI_Face_Destruct(self)
         implicit none
         class(MPI_Face_t) :: self

         self % no_of_faces = 0
         safedeallocate(self % faceIDs)
         safedeallocate(self % elementSide)
         safedeallocate(self % Nsend)
         safedeallocate(self % Nrecv)
         safedeallocate(self % Qsend)
         safedeallocate(self % U_xyzsend)
         safedeallocate(self % Qrecv)
         safedeallocate(self % U_xyzrecv)
         safedeallocate(self % AviscFluxSend)
         safedeallocate(self % AviscFluxRecv)
         safedeallocate(self % QBaseSend)
         safedeallocate(self % QBaseRecv)
         safedeallocate(self % transportVSend)
         safedeallocate(self % transportVRecv)
         safedeallocate(self % transportDSend)
         safedeallocate(self % transportDRecv)
#ifdef _HAS_MPI_
         self % Nrecv_req         = MPI_REQUEST_NULL
         self % Qrecv_req         = MPI_REQUEST_NULL
         self % gradQrecv_req     = MPI_REQUEST_NULL
         self % AviscFluxRecv_req = MPI_REQUEST_NULL
         self % QBaseRecv_req     = MPI_REQUEST_NULL
         self % transportVRecv_req = MPI_REQUEST_NULL
         self % transportDRecv_req = MPI_REQUEST_NULL
#endif

      end subroutine MPI_Face_Destruct

end module MPI_Face_Class
