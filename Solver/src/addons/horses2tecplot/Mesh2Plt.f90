module Mesh2PltModule
   use SMConstants
   use Storage
   use FileReadingUtilities      , only: getFileName
   implicit none

   private
   public   Mesh2Plt

#define PRECISION_FORMAT "(E18.10)"

   contains
      subroutine Mesh2Plt(meshFile, mode)
         use getTask, only: MODE_MULTIZONE, MODE_FINITEELM
         implicit none
         character(len=*), intent(in)     :: meshFile
         integer,          intent(in)     :: mode
!  
!        ---------------
!        Local variables   
!        ---------------
!
         integer                    :: eID, i, j, k, fid
         type(Mesh_t)               :: mesh
         character(len=LINE_LENGTH) :: meshPltName
!
!        Read the mesh from the *.hmesh file
!        -----------------------------------
         call mesh % ReadMesh(meshFile)
         mesh % isSurface = .false.
!
!        Create the tecplot mesh file name
!        ---------------------------------
         meshPltName = trim(getFileName(meshFile)) // ".tec"
         print *, "meshPltName: ", meshPltName
!
!        Open the file
!        -------------
         open(newunit=fid, file=trim(meshPltName), action="write", status="unknown")
!
!        Write the title and variables
!        -----------------------------
         write(fid,'(A,A,A)') 'TITLE = "',trim(meshFile),'"'
         write(fid,'(A)') 'VARIABLES = "x","y","z"'

!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM) then
            call WriteSingleFluidZoneToTecplot(fid,mesh)
         else
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
!
!              Write the tecplot file
!              ----------------------
               call WriteElementToTecplot(fid, e)
               end associate
            end do
         end if
!
!        Close the file
!        --------------
         close(fid)
   
      end subroutine Mesh2Plt


      subroutine WriteElementToTecplot(fid, e)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use OutputVariables
         use SolutionFile
         implicit none
         integer,            intent(in)    :: fid
         type(Element_t),    intent(inout) :: e 
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i,j,k,var
         character(len=LINE_LENGTH) :: formatout
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % Nmesh(1)+1,", J=",e % Nmesh(2)+1, &
                                            ", K=",e % Nmesh(3)+1,", F=POINT"

         formatout = getFormat()

         do k = 0, e % Nmesh(3)   ; do j = 0, e % Nmesh(2)    ; do i = 0, e % Nmesh(1)
            write(fid,trim(formatout)) e % x(:,i,j,k)
         end do               ; end do                ; end do

      end subroutine WriteElementToTecplot

!
!     Writes a single fluid zone using the FE Tecplot format
!     -> This format is more efficiently read by paraview and tecplot.
!     ------------------------------------------------------
      subroutine WriteSingleFluidZoneToTecplot(fid,mesh)
         use Storage
         use OutputVariables
         implicit none
         integer     , intent(in)        :: fid
         type(Mesh_t), intent(inout)     :: mesh
         !---------
         integer :: numOfPoints  ! Number of plot points
         integer :: numOfFElems  ! Number of FINITE elements
         integer :: firstPoint(size(mesh % elements))
         integer :: eID          ! Spectral (not finite) element counter
         integer :: i, j, k
         integer :: N(3)
         integer :: corners(8), cornersFace(4)
         character(len=LINE_LENGTH) :: formatout
         !---------
         
!        Definitions
!        -----------
         formatout = getFormat() ! format for point data
         
         ! Count points and elements
         numOfPoints   = product(mesh % elements(1) % Nmesh + 1)
         if (mesh % isSurface) then
             numOfFElems   = product(mesh % elements(1) % Nmesh(1:2))
         else
             numOfFElems   = product(mesh % elements(1) % Nmesh)
         end if
         firstPoint(1) = 1
         do eID = 2, size(mesh % elements)
            associate ( e => mesh % elements(eID) )
            firstPoint(eID) = numOfPoints + 1
            numOfPoints = numOfPoints + product(e % Nmesh + 1)
            if (mesh % isSurface) then
                numOfFElems   = numOfFElems + product(e % Nmesh(1:2))
            else
                numOfFElems   = numOfFElems + product(e % Nmesh)
            end if
            ! numOfFElems = numOfFElems + product(e % Nmesh    )
            end associate
         end do
         
         if (mesh % isSurface) then
             write(fid,'(A,I0,A,I0,A)') 'ZONE T="FLUID" N=',numOfPoints,' E=',numOfFElems,' ET=QUADRILATERAL, F=FEPOINT'
         else
             write(fid,'(A,I0,A,I0,A)') 'ZONE T="FLUID" N=',numOfPoints,' E=',numOfFElems,' ET=BRICK, F=FEPOINT'
         end if
         
!        Write the points
!        ----------------
         do eID = 1, size(mesh % elements)
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nmesh(3) ; do j = 0, e % Nmesh(2) ; do i = 0, e % Nmesh(1)
               write(fid,trim(formatout)) e % x(:,i,j,k)
            end do                ; end do                ; end do
            end associate
         end do
         
!        Write the elems connectivity
!        ----------------------------
         ! surface mesh case
         if (mesh % isSurface) then
             do eID = 1, size(mesh % elements)
                associate ( e => mesh % elements(eID) )

                do j = 0, e % Nmesh(2) - 1 ; do i = 0, e % Nmesh(1) - 1
                   cornersFace =  [ ij2localDOF(i,j,e%Nmesh(1:2)), ij2localDOF(i+1,j,e%Nmesh(1:2)), ij2localDOF(i+1,j+1,e%Nmesh(1:2)), ij2localDOF(i,j+1,e%Nmesh(1:2)) ] + firstPoint(eID)
                   write(fid,*) cornersFace
                end do                  ; end do

                end associate
             end do
         ! normal elements case
         else
             do eID = 1, size(mesh % elements)
                associate ( e => mesh % elements(eID) )
                
                do k = 0, e % Nmesh(3) - 1 ; do j = 0, e % Nmesh(2) - 1 ; do i = 0, e % Nmesh(1) - 1
                   corners =  [ ijk2localDOF(i,j,k  ,e%Nmesh), ijk2localDOF(i+1,j,k  ,e%Nmesh), ijk2localDOF(i+1,j+1,k  ,e%Nmesh), ijk2localDOF(i,j+1,k  ,e%Nmesh), &
                                ijk2localDOF(i,j,k+1,e%Nmesh), ijk2localDOF(i+1,j,k+1,e%Nmesh), ijk2localDOF(i+1,j+1,k+1,e%Nmesh), ijk2localDOF(i,j+1,k+1,e%Nmesh)  ] + firstPoint(eID)
                   write(fid,*) corners
                end do                    ; end do                    ; end do
                
                end associate
             end do
         end if
      end subroutine WriteSingleFluidZoneToTecplot

      character(len=LINE_LENGTH) function getFormat()
         use OutputVariables
         implicit none

         getFormat = ""

         write(getFormat,'(A,I0,A,A)') "(",3,PRECISION_FORMAT,")"

      end function getFormat

!     --------------------------------------------------------------
!     ij2localDOF:
!     Returns the local DOF index for a face in zero-based numbering
!     --------------------------------------------------------------
      function ij2localDOF(i,j,Nout) result(idx)
         implicit none
         
         integer, intent(in)   :: i, j, Nout(2)
         integer               :: idx
         
         IF (i < 0 .OR. i > Nout(1))     error stop 'error in ijk2local, i has wrong value'
         IF (j < 0 .OR. j > Nout(2))     error stop 'error in ijk2local, j has wrong value'
         
         idx = j*(Nout(1)+1) + i
      end function ij2localDOF

!     ------------------------------------------------------------------
!     ijk2localDOF:
!     Returns the local DOF index for an element in zero-based numbering
!     ------------------------------------------------------------------
      function ijk2localDOF(i,j,k,Nout) result(idx)
         implicit none
         
         integer, intent(in)   :: i, j, k, Nout(3)
         integer               :: idx
         
         IF (i < 0 .OR. i > Nout(1))     error stop 'error in ijk2local, i has wrong value'
         IF (j < 0 .OR. j > Nout(2))     error stop 'error in ijk2local, j has wrong value'
         IF (k < 0 .OR. k > Nout(3))     error stop 'error in ijk2local, k has wrong value'
         
         idx = k*(Nout(1)+1)*(Nout(2)+1) + j*(Nout(1)+1) + i
      end function ijk2localDOF

end module Mesh2PltModule