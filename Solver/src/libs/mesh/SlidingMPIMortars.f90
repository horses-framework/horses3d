!
!////////////////////////////////////////////////////////////////////////////////
!
!     SlidingMPIMortars
!     -----------------
!
!     Exchanges data across the sliding interface when the stator side and the
!     rotor side of a mortar are on different ranks (strategy 2).
!
!     Who does what
!        - The stator rank holds the full mortar: it has the geometry and it
!          computes the numerical flux.
!        - The rotor rank holds a light copy of the same mortar, without
!          geometry. It sends its solution and gradients, and receives the flux.
!        - This module only sends and receives values. The mortars themselves
!          are built by SlidingConstructSplitMortars, in SlidingMeshProcedures.
!
!     How the two ranks agree without sending any index
!        Both ranks build the same list of shared mortars, in the same order,
!        from the interface directory (identical on every rank) and the current
!        sector. A message therefore carries values only: the k-th block received
!        belongs to the k-th mortar of the list. If the two lists were ever in a
!        different order, the values would land on the wrong mortars and nothing
!        would crash. SlidingCheckCommAlignment checks this.
!
!     Who calls whom
!     --------------
!        setup, once
!         |
!         +-- SlidingConstructMPIMortars
!               +-- SlidingSizeMortarBuffers        sweeps the N sectors
!               +-- ConstructMPIFaces               allocates the per-rank array
!               +-- ConstructMPIFacesStorage        allocates the buffers
!
!        every sector change
!         |
!         +-- SlidingBuildCommLists
!         |     +-- RowPos
!         |     +-- SortCommLinks
!         +-- SlidingCheckCommAlignment
!
!        every Runge-Kutta stage, in this order
!         |
!         +-- SlidingUpdateMPIMortarSolution / Gather...      rotor -> stator
!         +-- SlidingUpdateMPIMortarGradFlux / Gather...      stator -> rotor
!         +-- SlidingUpdateMPIMortarGradients / Gather...     rotor -> stator
!         +-- SlidingUpdateMPIMortarFlux     / Gather...      stator -> rotor
!
!        end of run
!         |
!         +-- SlidingDestructMPIMortars
!
!////////////////////////////////////////////////////////////////////////////////
!

#include "Includes.h"

module SlidingMPIMortars
   use SMConstants
   use SlidingMeshClass
   use HexMeshClass
   use MPI_Face_Class
   use MPI_Process_Info
   use PhysicsStorage
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   private

   public :: SlidingConstructMPIMortars, SlidingDestructMPIMortars
   public :: SlidingBuildCommLists, SlidingCheckCommAlignment
   public :: SlidingUpdateMPIMortarSolution,  SlidingGatherMPIMortarSolution
   public :: SlidingUpdateMPIMortarFlux,      SlidingGatherMPIMortarFlux
   public :: SlidingUpdateMPIMortarGradients, SlidingGatherMPIMortarGradients
   public :: SlidingUpdateMPIMortarGradFlux,  SlidingGatherMPIMortarGradFlux

!  ---------------------------------------------------------------------------
!   One mortar shared with another rank.
!
!   ringRow and Mortarpos identify the mortar on both ranks: each rank finds
!   them in the directory, for the current sector. mortarIndex is the position
!   of the mortar in the local mortar_faces, so it differs between the two
!   ranks and is never sent.
!
!   iAmOwner is true on the stator rank, which holds the full mortar and
!   computes the flux.
!  ---------------------------------------------------------------------------
   type SlidingMortarLink
      integer :: otherRank   = -1      ! rank holding the other side
      integer :: ringRow     = 0       ! STATIC ring row, the common reference
      integer :: Mortarpos   = 0       ! 0 or 1
      integer :: mortarIndex = 0       ! local index into mesh % mortar_faces
      integer :: nDOF        = 0       ! (Nf(1)+1) * (Nf(2)+1) on this mortar
      logical :: iAmOwner    = .false. ! static side: builds geometry, solves Riemann
   end type SlidingMortarLink

   type(SlidingMortarLink), allocatable :: mortarComm(:)
   integer                              :: nComm    = 0
   integer                              :: nCommMax = 0
   type(MPI_FacesSet_t)                 :: slidingFacesSet

contains
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingSizeMortarBuffers
!
!  For each partner rank, the largest number of values exchanged over all
!  sectors. The buffers are allocated once at this size and never resized.
!
!  The mortars do not exist yet when this runs, so their size is taken from the
!  polynomial orders stored in the directory. This assumes a uniform order on
!  the interface, which SlidingBuildGlobalRings checks.
! ---------------------------------------------------------------------------
subroutine SlidingSizeMortarBuffers(mesh, nDOFsMax)
   implicit none
   class(HexMesh), intent(in)  :: mesh
   integer,        intent(out) :: nDOFsMax(MPI_Process % nProcs)

   integer :: s, j, r, Mortarpos, rows(0:1), d, n
   integer :: nDOFs(MPI_Process % nProcs)

   associate(SM => mesh % SlidingMesh)

      nDOFsMax = 0
      nCommMax = 0

      do s = 0, SM % numElemsPerLayer - 1
         nDOFs = 0
         n = 0

         do j = 1, SM % nRingGlobal
            call SlidingSlaveRows(SM, j, s, rows)

            do Mortarpos = 0, 1
               r = rows(Mortarpos)

               if (SM % statorLocalElem(j) >  0 .and. SM % rotorLocalElem(r) == 0) then
                  d = SM % rotorRing(r)  % rank + 1
               else if (SM % statorLocalElem(j) == 0 .and. SM % rotorLocalElem(r) >  0) then
                  d = SM % statorRing(j) % rank + 1
               else
                  cycle
               end if

               nDOFs(d) = nDOFs(d) + product(SM % statorRing(j) % Nel + 1)
               n = n + 1
            end do
         end do

         nDOFsMax = max(nDOFsMax, nDOFs)
         nCommMax = max(nCommMax, n)
      end do

   end associate
end subroutine SlidingSizeMortarBuffers
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingConstructMPIMortars
!
!  Sets up the partner ranks and the buffers of the sliding exchanges. A
!  partner is any rank that shares at least one mortar with this rank, in any
!  sector. Called once, at setup.
!
!  nDOFsMax is passed twice, as NDOFS and NMDOFS: ConstructMPIFacesStorage only
!  allocates the flux buffers when NMDOFS is present.
! ---------------------------------------------------------------------------
subroutine SlidingConstructMPIMortars(mesh)
   implicit none
   class(HexMesh), intent(inout) :: mesh

   integer :: d, nDOFsMax(MPI_Process % nProcs), sharedDomain(MPI_Process % nProcs)

   if (.not. MPI_Process % doMPIAction) return

   call SlidingSizeMortarBuffers(mesh, nDOFsMax)

   call ConstructMPIFaces(slidingFacesSet)

   slidingFacesSet % nDomainShared = 0
   sharedDomain = 0

   do d = 1, MPI_Process % nProcs
      if (nDOFsMax(d) <= 0) cycle
      call slidingFacesSet % faces(d) % Construct(nCommMax)
      slidingFacesSet % nDomainShared = slidingFacesSet % nDomainShared + 1
      sharedDomain(slidingFacesSet % nDomainShared) = d
   end do

   safedeallocate(slidingFacesSet % listDomain)
   allocate(slidingFacesSet % listDomain(slidingFacesSet % nDomainShared))
   slidingFacesSet % listDomain = sharedDomain(1:slidingFacesSet % nDomainShared)

#if defined(NAVIERSTOKES)
   call ConstructMPIFacesStorage(slidingFacesSet, NCONS, NGRAD, nDOFsMax, nDOFsMax)
#elif defined(INCNS)
   call ConstructMPIFacesStorage(slidingFacesSet, NCONS, NCONS, nDOFsMax, nDOFsMax)
#elif defined(CAHNHILLIARD) && !defined(MULTIPHASE)
   call ConstructMPIFacesStorage(slidingFacesSet, NCOMP, NCOMP, nDOFsMax, nDOFsMax)
#elif defined(MULTIPHASE)
   call ConstructMPIFacesStorage(slidingFacesSet, NCONS, NCONS, nDOFsMax, nDOFsMax)
#endif

   safedeallocate(mortarComm)
   allocate(mortarComm(max(nCommMax,1)))
   nComm = 0

   slidingFacesSet % constructed = .true.

end subroutine SlidingConstructMPIMortars
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingBuildCommLists
!
!  Lists the mortars this rank shares with other ranks in the current sector,
!  and sorts them. Both ranks of a pair sort with the same keys, so they get
!  the same order without sending anything:
!
!    otherRank   groups the mortars by partner, one message per partner
!    Ax_Id       axial position of the stator face
!    Azim_Id     azimuthal position of the stator face
!    Mortarpos   which of the two mortars of that face
!
!  Called at setup and at every sector change. No communication.
! ---------------------------------------------------------------------------
subroutine SlidingBuildCommLists(mesh, sectorID)
   implicit none
   class(HexMesh), intent(in) :: mesh
   integer,        intent(in) :: sectorID

   integer :: j, r, m, Mortarpos, rows(0:1), k, d, domain
   integer :: no_of_links(MPI_Process % nProcs)

   if (.not. MPI_Process % doMPIAction) return

   associate(SM => mesh % SlidingMesh)

      nComm = 0

      do j = 1, SM % nRingGlobal
         call SlidingSlaveRows(SM, j, sectorID, rows)

         do Mortarpos = 0, 1
            r = rows(Mortarpos)

            if (SM % statorLocalElem(j) >  0 .and. SM % rotorLocalElem(r) == 0) then
               ! I own the static side: I hold the full mortar
               m = RowPos(SM % myMasterRows, SM % nMasterRows, j)
               nComm = nComm + 1
               mortarComm(nComm) % otherRank = SM % rotorRing(r) % rank
               mortarComm(nComm) % mortarIndex = MortarSlot(SM, Mortarpos, m)
               mortarComm(nComm) % iAmOwner = .true.

            else if (SM % statorLocalElem(j) == 0 .and. SM % rotorLocalElem(r) >  0) then
               ! I own the rotating side: light copy, the full mortar is elsewhere
               m = RowPos(SM % mySlaveRows, SM % nSlaveRows, j)
               nComm = nComm + 1
               mortarComm(nComm) % otherRank = SM % statorRing(j) % rank
               mortarComm(nComm) % mortarIndex = MortarSlot(SM, Mortarpos, m)
               mortarComm(nComm) % iAmOwner = .false.

            else
               cycle
            end if

            mortarComm(nComm) % ringRow = j
            mortarComm(nComm) % Mortarpos = Mortarpos
            mortarComm(nComm) % nDOF = &
            product(mesh % mortar_faces(mortarComm(nComm) % mortarIndex) % Nf + 1)
         end do
      end do

      call SortCommLinks(SM, mortarComm, nComm)

      ! Number of shared mortars per partner, for this sector
      no_of_links = 0
      do k = 1, nComm
         no_of_links(mortarComm(k) % otherRank + 1) = &
              no_of_links(mortarComm(k) % otherRank + 1) + 1
      end do

      do d = 1, slidingFacesSet % nDomainShared
         domain = slidingFacesSet % listDomain(d)
         slidingFacesSet % faces(domain) % no_of_faces = no_of_links(domain)
      end do

   end associate
end subroutine SlidingBuildCommLists
!
!////////////////////////////////////////////////////////////////////////////////
!  
! Position of ring row j in a list of owned rows. Zero when absent.
pure integer function RowPos(rows, n, j) result(m)
   implicit none
   integer, intent(in) :: rows(:), n, j
   integer :: k

   m = 0
   do k = 1, n
      if (rows(k) == j) then
         m = k
         return
      end if
   end do

end function RowPos
!
!////////////////////////////////////////////////////////////////////////////////
!  
! Insertion sort of the links, with the keys of LinkBefore.
subroutine SortCommLinks(SM, link, n)
   implicit none
   type(SlidingMesh),       intent(in)    :: SM
   type(SlidingMortarLink), intent(inout) :: link(:)
   integer,                 intent(in)    :: n

   integer                 :: i, k
   type(SlidingMortarLink) :: key

   do i = 2, n
      key = link(i)
      k = i-1
      do while (k >= 1)
         if (.not. LinkBefore(SM, key, link(k))) exit
         link(k+1) = link(k)
         k = k-1
      end do
      link(k+1) = key
   end do
end subroutine SortCommLinks
!
!////////////////////////////////////////////////////////////////////////////////
!  
! True if link a comes before link b. Ax_Id and Azim_Id are read from the
! directory, not stored in the link.
pure logical function LinkBefore(SM, a, b) result(before)
   implicit none
   type(SlidingMesh),       intent(in) :: SM
   type(SlidingMortarLink), intent(in) :: a, b
   integer :: aAx, bAx, aAz, bAz

   before = .false.

   if (a % otherRank /= b % otherRank) then
      before = (a % otherRank < b % otherRank) 
      return
   end if

   aAx = SM % RingAxId(a % ringRow) 
   bAx = SM % RingAxId(b % ringRow)
   if (aAx /= bAx) then
      before = (aAx < bAx) 
      return
   end if

   aAz = SM % RingAzimId(a % ringRow) 
   bAz = SM % RingAzimId(b % ringRow)
   if (aAz /= bAz) then
      before = (aAz < bAz) 
      return
   end if

   before = (a % Mortarpos < b % Mortarpos)
end function LinkBefore
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingCheckCommAlignment
!
!  Sends the sorting keys and the size of each shared mortar to the partner,
!  and stops the run if the two lists differ. Without this check, a mismatch
!  would put values on the wrong mortars without any error.
!
!  All receives are posted before any send, so the order in which ranks walk
!  their partners cannot cause a deadlock.
! ---------------------------------------------------------------------------
subroutine SlidingCheckCommAlignment(mesh)
   implicit none
   class(HexMesh), intent(in) :: mesh
#ifdef _HAS_MPI_
   integer, parameter   :: NKEY = 4      ! Ax_Id, Azim_Id, Mortarpos, nDOF
   integer :: d, domain, i, k, n, ierr, nShared, nreqs, idx_send
   integer, allocatable :: sendbuf(:,:), recvbuf(:,:), all_reqs(:)

   if (.not. MPI_Process % doMPIAction)     return
   if (.not. slidingFacesSet % constructed) return

   nShared = slidingFacesSet % nDomainShared
   if (nShared <= 0) return

   allocate(sendbuf(NKEY*nCommMax, nShared))
   allocate(recvbuf(NKEY*nCommMax, nShared))
   sendbuf = 0
   recvbuf = 0

   nreqs = 2 * nShared
   allocate(all_reqs(nreqs))
   all_reqs = MPI_REQUEST_NULL

   associate(SM => mesh % SlidingMesh)

   ! Post the receives
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      n = slidingFacesSet % faces(domain) % no_of_faces
      if (n <= 0) cycle

      call mpi_irecv(recvbuf(1,d), NKEY*n, MPI_INTEGER, domain-1, 0, &
      MPI_COMM_WORLD, all_reqs(d), ierr)
   end do

   ! Pack the keys, in list order, and send
   idx_send = nShared + 1
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      n = slidingFacesSet % faces(domain) % no_of_faces

      if (n <= 0) then
         idx_send = idx_send + 1
         cycle
      end if

      k = 0
      do i = 1, nComm
         if (mortarComm(i) % otherRank + 1 /= domain) cycle
         sendbuf(NKEY*k+1, d) = SM % RingAxId  (mortarComm(i) % ringRow)
         sendbuf(NKEY*k+2, d) = SM % RingAzimId(mortarComm(i) % ringRow)
         sendbuf(NKEY*k+3, d) = mortarComm(i) % Mortarpos
         sendbuf(NKEY*k+4, d) = mortarComm(i) % nDOF
         k = k + 1
      end do

      if (k /= n) then
         write(STD_OUT,*) 'sliding: rank', MPI_Process % rank, 'counted', k, &
         'links towards rank', domain-1, 'but no_of_faces says', n
         error stop
      end if

      call mpi_isend(sendbuf(1,d), NKEY*n, MPI_INTEGER, domain-1, 0, &
         MPI_COMM_WORLD, all_reqs(idx_send), ierr)
      idx_send = idx_send + 1
   end do

   call MPI_Waitall(nreqs, all_reqs, MPI_STATUSES_IGNORE, ierr)

   ! Compare term by term
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      n = slidingFacesSet % faces(domain) % no_of_faces
      if (n <= 0) cycle

      do k = 1, n
         if (any(sendbuf(NKEY*(k-1)+1:NKEY*k, d) /= recvbuf(NKEY*(k-1)+1:NKEY*k, d))) then
            write(STD_OUT,*) 'sliding: comm lists out of step between ranks', &
                              MPI_Process % rank, 'and', domain-1, 'at link', k
            write(STD_OUT,*) '   mine   (Ax_Id, Azim_Id, Mortarpos, nDOF) =', &
                              sendbuf(NKEY*(k-1)+1:NKEY*k, d)
            write(STD_OUT,*) '   theirs (Ax_Id, Azim_Id, Mortarpos, nDOF) =', &
                              recvbuf(NKEY*(k-1)+1:NKEY*k, d)
            error stop
         end if
      end do
   end do

   end associate

   deallocate(sendbuf, recvbuf, all_reqs)
#endif
end subroutine SlidingCheckCommAlignment
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingUpdateMPIMortarSolution
!
!  Rotor -> stator. The rotor rank sends the solution on its mortars, the
!  stator rank receives it. The values are already on the mortar points, so
!  nothing needs to be interpolated on arrival.
! ---------------------------------------------------------------------------
subroutine SlidingUpdateMPIMortarSolution(mesh, nEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn

#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter, ierr, nShared, nreqs, idx_send
   integer, allocatable :: all_reqs(:)

   if (.not. MPI_Process % doMPIAction) return

   nShared = slidingFacesSet % nDomainShared
   if (nShared <= 0) return

   nreqs = 2 * nShared
   allocate(all_reqs(nreqs))
   all_reqs = MPI_REQUEST_NULL

   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces > 0) then
         call slidingFacesSet % faces(domain) % RecvQ(domain, nEqn, all_reqs(d))
      end if
   end do

   idx_send = nShared + 1
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)

      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) then
         idx_send = idx_send + 1
         cycle
      end if

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (mortarComm(k) % iAmOwner) cycle   ! the owner receives

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            slidingFacesSet % faces(domain) % Qsend(counter:counter+nEqn-1) = &
                  f % storage(2) % Q(1:nEqn,i,j)
            counter = counter + nEqn
         end do; end do
         end associate
      end do

      call slidingFacesSet % faces(domain) % SendQ(domain, nEqn, all_reqs(idx_send))
      idx_send = idx_send + 1
   end do

   call MPI_Waitall(nreqs, all_reqs, MPI_STATUSES_IGNORE, ierr)
   deallocate(all_reqs)
#endif

end subroutine SlidingUpdateMPIMortarSolution
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingGatherMPIMortarSolution
!
!  Copies the received rotor solution into the rotor side of each full mortar.
!  The k-th block goes to the k-th mortar of the sorted list.
! ---------------------------------------------------------------------------
subroutine SlidingGatherMPIMortarSolution(mesh, nEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn

#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter

   if (.not. MPI_Process % doMPIAction) return

   do d = 1, slidingFacesSet % nDomainShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) cycle

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (.not. mortarComm(k) % iAmOwner) cycle   ! only the owner receives

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            f % storage(2) % Q(1:nEqn,i,j) = slidingFacesSet % faces(domain) % Qrecv(counter:counter+nEqn-1)
            counter = counter + nEqn
         end do; end do
         end associate
      end do
   end do
#endif

end subroutine SlidingGatherMPIMortarSolution
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingUpdateMPIMortarFlux
!
!  Stator -> rotor. The stator rank sends the rotor-side flux it has just
!  computed on each shared mortar.
! ---------------------------------------------------------------------------
subroutine SlidingUpdateMPIMortarFlux(mesh, nEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter, ierr, nShared, nreqs, idx_send
   integer, allocatable :: all_reqs(:)

   if (.not. MPI_Process % doMPIAction) return

   nShared = slidingFacesSet % nDomainShared
   if (nShared <= 0) return

   nreqs = 2 * nShared
   allocate(all_reqs(nreqs))
   all_reqs = MPI_REQUEST_NULL

   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces > 0) then
         call slidingFacesSet % faces(domain) % RecvMortarFlux(domain, nEqn, all_reqs(d))
      end if
   end do

   idx_send = nShared + 1
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)

      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) then
         idx_send = idx_send + 1
         cycle
      end if

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (.not. mortarComm(k) % iAmOwner) cycle   ! only the owner sends a flux

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            slidingFacesSet % faces(domain) % Flux_M_Send(counter:counter+nEqn-1) = &
                  f % storage(2) % fStar(1:nEqn,i,j)
            counter = counter + nEqn
         end do; end do
         end associate
      end do

      call slidingFacesSet % faces(domain) % SendMortarFlux(domain, nEqn, all_reqs(idx_send))
      idx_send = idx_send + 1
   end do

   call MPI_Waitall(nreqs, all_reqs, MPI_STATUSES_IGNORE, ierr)
   deallocate(all_reqs)
#endif
end subroutine SlidingUpdateMPIMortarFlux
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingGatherMPIMortarFlux
!
!  Projects each received flux onto the rotor face. A rotor face is covered by
!  two mortars, so each projection is added to the face, never copied over it.
! ---------------------------------------------------------------------------
subroutine SlidingGatherMPIMortarFlux(mesh, nEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn

#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter

   if (.not. MPI_Process % doMPIAction) return

   do d = 1, slidingFacesSet % nDomainShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) cycle

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (mortarComm(k) % iAmOwner) cycle   ! the owner does not receive

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            f % storage(2) % fStar(1:nEqn,i,j) = &
                  slidingFacesSet % faces(domain) % Flux_M_Recv(counter:counter+nEqn-1)
            counter = counter + nEqn
         end do; end do

         call mesh % faces(f % Mortar(2)) % ProjectMortarFluxToElements(nEqn, [2, 0], f, &
                                       f % storage(2) % fStar(1:nEqn,:,:), sliding = .true.)
         end associate
      end do
   end do
#endif

end subroutine SlidingGatherMPIMortarFlux
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingUpdateMPIMortarGradients
!
!  Rotor -> stator, like the solution: the stator rank needs the rotor
!  gradients on its mortars to compute the viscous flux.
! ---------------------------------------------------------------------------
subroutine SlidingUpdateMPIMortarGradients(mesh, nGradEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nGradEqn
   
#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter, ierr, nShared, nreqs, idx_send
   integer, allocatable :: all_reqs(:)

   if (.not. MPI_Process % doMPIAction) return

   nShared = slidingFacesSet % nDomainShared
   if (nShared <= 0) return

   nreqs = 2 * nShared
   allocate(all_reqs(nreqs))
   all_reqs = MPI_REQUEST_NULL

   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces > 0) then
         call slidingFacesSet % faces(domain) % RecvU_xyz(domain, nGradEqn, all_reqs(d))
      end if
   end do

   idx_send = nShared + 1
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)

      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) then
         idx_send = idx_send + 1
         cycle
      end if

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (mortarComm(k) % iAmOwner) cycle   ! the owner receives

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex), &
                     buf => slidingFacesSet % faces(domain) % U_xyzsend)
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            buf(counter              : counter+  nGradEqn-1) = f % storage(2) % U_x(1:nGradEqn,i,j)
            buf(counter+  nGradEqn   : counter+2*nGradEqn-1) = f % storage(2) % U_y(1:nGradEqn,i,j)
            buf(counter+2*nGradEqn   : counter+3*nGradEqn-1) = f % storage(2) % U_z(1:nGradEqn,i,j)
            counter = counter + 3*nGradEqn
         end do; end do
         end associate
      end do

      call slidingFacesSet % faces(domain) % SendU_xyz(domain, nGradEqn, all_reqs(idx_send))
      idx_send = idx_send + 1
   end do

   call MPI_Waitall(nreqs, all_reqs, MPI_STATUSES_IGNORE, ierr)
   deallocate(all_reqs)
#endif

end subroutine SlidingUpdateMPIMortarGradients
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingGatherMPIMortarGradients
!
!  Copies the received rotor gradients into the rotor side of each full
!  mortar, in the order of the sorted list.
! ---------------------------------------------------------------------------
subroutine SlidingGatherMPIMortarGradients(mesh, nGradEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nGradEqn
#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter

   if (.not. MPI_Process % doMPIAction) return

   do d = 1, slidingFacesSet % nDomainShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) cycle

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (.not. mortarComm(k) % iAmOwner) cycle   ! only the owner receives

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex), &
                     buf => slidingFacesSet % faces(domain) % U_xyzrecv)
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            f % storage(2) % U_x(1:nGradEqn,i,j) = buf(counter              : counter+  nGradEqn-1)
            f % storage(2) % U_y(1:nGradEqn,i,j) = buf(counter+  nGradEqn   : counter+2*nGradEqn-1)
            f % storage(2) % U_z(1:nGradEqn,i,j) = buf(counter+2*nGradEqn   : counter+3*nGradEqn-1)
            counter = counter + 3*nGradEqn
         end do; end do
         end associate
      end do
   end do
#endif
end subroutine SlidingGatherMPIMortarGradients
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingUpdateMPIMortarGradFlux
!
!  Stator -> rotor. The stator rank sends the BR1 gradient flux it has kept
!  for the rotor side of each shared mortar.
! ---------------------------------------------------------------------------
subroutine SlidingUpdateMPIMortarGradFlux(mesh, nEqn)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter, ierr, nShared, nreqs, idx_send, nBlock
   integer, allocatable :: all_reqs(:)

   if (.not. MPI_Process % doMPIAction) return

   nShared = slidingFacesSet % nDomainShared
   if (nShared <= 0) return

   nBlock = nEqn * NDIM
   nreqs = 2 * nShared
   allocate(all_reqs(nreqs))
   all_reqs = MPI_REQUEST_NULL

   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces > 0) then
         call slidingFacesSet % faces(domain) % RecvMortarGradFlux(domain, nEqn, all_reqs(d))
      end if
   end do

   idx_send = nShared + 1
   do d = 1, nShared
      domain = slidingFacesSet % listDomain(d)

      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) then
         idx_send = idx_send + 1
         cycle
      end if

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (.not. mortarComm(k) % iAmOwner) cycle   ! only the owner sends

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            slidingFacesSet % faces(domain) % GradFlux_M_Send(counter:counter+nBlock-1) = &
                  reshape(f % storage(2) % unStar(1:nEqn,1:NDIM,i,j), [nBlock])
            counter = counter + nBlock
         end do; end do
         end associate
      end do

      call slidingFacesSet % faces(domain) % SendMortarGradFlux(domain, nEqn, all_reqs(idx_send))
      idx_send = idx_send + 1
   end do

   call MPI_Waitall(nreqs, all_reqs, MPI_STATUSES_IGNORE, ierr)
   deallocate(all_reqs)
#endif
end subroutine SlidingUpdateMPIMortarGradFlux
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingGatherMPIMortarGradFlux
!
!  Projects each received gradient flux onto the rotor face, adding to the
!  contribution of the other mortar. factor is the same as in the local call
!  for the rotor side.
! ---------------------------------------------------------------------------
subroutine SlidingGatherMPIMortarGradFlux(mesh, nEqn, factor)
   implicit none
   class(HexMesh)      :: mesh
   integer, intent(in) :: nEqn
   integer, intent(in) :: factor
#ifdef _HAS_MPI_
   integer :: k, d, domain, i, j, counter, nBlock

   if (.not. MPI_Process % doMPIAction) return

   nBlock = nEqn * NDIM

   do d = 1, slidingFacesSet % nDomainShared
      domain = slidingFacesSet % listDomain(d)
      if (slidingFacesSet % faces(domain) % no_of_faces <= 0) cycle

      counter = 1
      do k = 1, nComm
         if (mortarComm(k) % otherRank + 1 /= domain) cycle
         if (mortarComm(k) % iAmOwner) cycle   ! the owner does not receive

         associate(f => mesh % mortar_faces(mortarComm(k) % mortarIndex))
         do j = 0, f % Nf(2); do i = 0, f % Nf(1)
            f % storage(2) % unStar(1:nEqn,1:NDIM,i,j) = &
                  reshape(slidingFacesSet % faces(domain) % GradFlux_M_Recv(counter:counter+nBlock-1), [nEqn, NDIM])
            counter = counter + nBlock
         end do; end do

         call mesh % faces(f % Mortar(2)) % ProjectMortarGradientFluxToElements(nEqn = nEqn, slaveFace = f, &
                  Hflux = f % storage(2) % unStar(1:nEqn,1:NDIM,:,:), whichElements = [2, 0], &
                  factor = factor, sliding = .true.)
         end associate
      end do
   end do
#endif
end subroutine SlidingGatherMPIMortarGradFlux
!
!////////////////////////////////////////////////////////////////////////////////
!  
! ---------------------------------------------------------------------------
!  SlidingDestructMPIMortars
!
!  Frees the buffers of the sliding exchanges. Called at the end of the run,
!  from DestructDGSem.
! ---------------------------------------------------------------------------
subroutine SlidingDestructMPIMortars()
   implicit none

   safedeallocate(mortarComm)
   nComm = 0
   nCommMax = 0

   if (slidingFacesSet % constructed) then
      call DestructMPIFaces(slidingFacesSet)
   end if
   slidingFacesSet % constructed = .false.

end subroutine SlidingDestructMPIMortars

end module SlidingMPIMortars