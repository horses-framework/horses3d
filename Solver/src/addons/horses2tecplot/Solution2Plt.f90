module Solution2PltModule
   use SMConstants
   use SolutionFile
   use InterpolationMatrices
   use FileReadingUtilities      , only: getFileName
   use getTask
   implicit none

   private
   public   Solution2Plt, WriteBoundaryToTecplot

#define PRECISION_FORMAT "(E18.10)"

   contains
      subroutine Solution2Plt(meshName, solutionName, fixedOrder, basis, Nout, mode)
         use getTask
         use Headers
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: basis
         logical,          intent(in)     :: fixedOrder
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode

         write(STD_OUT,'(/)')
         call SubSection_Header("Job description")
         
         select case (mode)
         case(MODE_FINITEELM)
            write(STD_OUT,'(30X,A3,A)') "->", " Output mode: Tecplot FE"
         case(MODE_MULTIZONE)
            write(STD_OUT,'(30X,A3,A)') "->", " Output mode: Tecplot Multi-Zone"
         end select
         
         select case ( basis )

         case(EXPORT_GAUSS)

            if ( fixedOrder ) then
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points with fixed order"
               write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                                Nout(1),",",Nout(2),",",Nout(3),"]."
               call Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout, mode)
   
            else
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points"
               call Solution2Plt_GaussPoints(meshName, solutionName, mode)

            end if

         case(EXPORT_HOMOGENEOUS)
            
            write(STD_OUT,'(30X,A3,A)') "->", " Export to homogeneous points"
            write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                        Nout(1),",",Nout(2),",",Nout(3),"]."
            call Solution2Plt_Homogeneous(meshName, solutionName, Nout, mode)

         end select

      end subroutine Solution2Plt
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_GaussPoints(meshName, solutionName, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                    :: mesh
         character(len=LINE_LENGTH)      :: meshPltName
         character(len=LINE_LENGTH)      :: solutionFile
         character(len=1024)             :: title
         integer                         :: no_of_elements, eID
         integer                         :: fid, bID
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
         no_of_elements = mesh % no_of_elements
!
!        Transform zones to the output variables
!        ---------------------------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = e % Nsol
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Project mesh and solution
!           -------------------------
            call ProjectStorageGaussPoints(e, spA, e % Nmesh, e % Nsol, mesh % hasGradients, mesh % isStatistics)

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         call getOutputVariables()

         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM) then
            call WriteSingleFluidZoneToTecplot(fid,mesh)
         else
            do eID = 1, no_of_elements
               associate ( e => mesh % elements(eID) )
!
!              Write the tecplot file
!              ----------------------
               call WriteElementToTecplot(fid, e, mesh % refs, &
                                          mesh % hasGradients, mesh % isStatistics, mesh % hasSensor)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if
!
!        Close the file
!        --------------
         close(fid)
      
      end subroutine Solution2Plt_GaussPoints

      subroutine ProjectStorageGaussPoints(e, spA, NM, NS, hasGradients, hasStats)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t),  intent(in)  :: spA(0:)
         integer           ,  intent(in)  :: NM(3)
         integer           ,  intent(in)  :: NS(3)
         logical,             intent(in)  :: hasGradients
         logical,             intent(in)  :: hasStats
         
         e % Nout = e % Nsol
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, NM, NS)

         end if

         e % Qout(1:,0:,0:,0:) => e % Q
         if ( hasGradients ) then
            e % U_xout(1:,0:,0:,0:) => e % U_x
            e % U_yout(1:,0:,0:,0:) => e % U_y
            e % U_zout(1:,0:,0:,0:) => e % U_z
         end if
         e % QDot_out(1:,0:,0:,0:) => e % QDot
         if (hasStats) e % statsout(1:,0:,0:,0:) => e % stats
         if (hasUt_NS) e % ut_NSout(1:,0:,0:,0:) => e % ut_NS
         if (hasMu_NS) e % mu_NSout(1:,0:,0:,0:) => e % mu_NS
         if (hasWallY) e % wallYout(1:,0:,0:,0:) => e % wallY
         if (hasMu_sgs) e % mu_sgsout(1:,0:,0:,0:) => e % mu_sgs
         if (hasSource) e % sourceout(1:,0:,0:,0:) => e % source

      end subroutine ProjectStorageGaussPoints
!
!//////////////////////////////////////////////////////////////////////////////////
!
!     Gauss points with fixed order procedures
!     ----------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Allocate the output spectral basis
!        ----------------------------------
         call spA(Nout(1)) % Construct(GAUSS, Nout(1))
         call spA(Nout(2)) % Construct(GAUSS, Nout(2))
         call spA(Nout(3)) % Construct(GAUSS, Nout(3))
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol, mesh % nodeType)
!
!           Construct interpolation matrices
!           --------------------------------
            associate( spAoutXi   => spA(Nout(1)), &
                       spAoutEta  => spA(Nout(2)), &
                       spAoutZeta => spA(Nout(3)) )
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), spAoutXi   % x)   ! TODO: check why it was Nsol(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), spAoutEta  % x)   ! TODO: check why it was Nsol(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), spAoutZeta % x)   ! TODO: check why it was Nsol(1)
            end associate
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageGaussPoints_FixedOrder(e, spA, e % Nmesh, e % Nsol, e % Nout, &
                                                                    Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                                    Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                                    Tset(e % Nout(3), e % Nsol(3)) % T, &
                                                                    mesh % hasGradients, mesh % isStatistics, mesh % hasTimeDeriv)

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         call getOutputVariables()
         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM) then
            call WriteSingleFluidZoneToTecplot(fid,mesh)
         else
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )

               call WriteElementToTecplot(fid, e, mesh % refs, &
                                          mesh % hasGradients, mesh % isStatistics, mesh % hasSensor)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if
!
!        Close the file
!        --------------
         close(fid)

      end subroutine Solution2Plt_GaussPoints_FixedOrder

      subroutine ProjectStorageGaussPoints_FixedOrder(e, spA, NM, NS, Nout, Tx, Ty, Tz, hasGradients, hasStats, hasTimeDeriv)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t),  intent(in)  :: spA(0:)
         integer           ,  intent(in)  :: NM(3)
         integer           ,  intent(in)  :: NS(3)
         integer           ,  intent(in)  :: Nout(3)
         real(kind=RP),       intent(in)  :: Tx(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: Ty(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: Tz(0:e % Nout(3), 0:e % Nsol(3))
         logical,             intent(in)  :: hasGradients
         logical,             intent(in)  :: hasStats
         logical,             intent(in)  :: hasTimeDeriv
!
!        Project mesh
!        ------------         
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, NM, Nout)

         end if
!
!        Project the solution
!        --------------------
         if ( all( e % Nsol .eq. e % Nout ) ) then
            e % Qout(1:,0:,0:,0:) => e % Q

            if ( hasGradients ) then
               e % U_xout(1:,0:,0:,0:) => e % U_x
               e % U_yout(1:,0:,0:,0:) => e % U_y
               e % U_zout(1:,0:,0:,0:) => e % U_z
            end if


            e % QDot_out(1:,0:,0:,0:) => e % QDot
            if (hasStats) e % statsout(1:,0:,0:,0:) => e % stats
            if (hasUt_NS) e % ut_NSout(1:,0:,0:,0:) => e % ut_NS
            if (hasMu_NS) e % mu_NSout(1:,0:,0:,0:) => e % mu_NS
            if (hasWallY) e % wallYout(1:,0:,0:,0:) => e % wallY
            if (hasMu_sgs) e % mu_sgsout(1:,0:,0:,0:) => e % mu_sgs
            if (hasSource) e % sourceout(1:,0:,0:,0:) => e % source

         else
            allocate( e % Qout(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongSolutionToGaussPoints(NVARS, e % Nsol, e % Q, e % Nout, e % Qout, Tx, Ty, Tz)
   
            if ( hasGradients ) then
               allocate( e % U_xout(1:NGRADVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
               allocate( e % U_yout(1:NGRADVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
               allocate( e % U_zout(1:NGRADVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
               call prolongSolutionToGaussPoints(NGRADVARS, e % Nsol, e % U_x, e % Nout, e % U_xout, Tx, Ty, Tz)
               call prolongSolutionToGaussPoints(NGRADVARS, e % Nsol, e % U_y, e % Nout, e % U_yout, Tx, Ty, Tz)
               call prolongSolutionToGaussPoints(NGRADVARS, e % Nsol, e % U_z, e % Nout, e % U_zout, Tx, Ty, Tz)
            end if

            if (hasStats) then
                allocate( e % statsout(1:NSTAT,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(NSTAT, e % Nsol, e % stats, e % Nout, e % statsout, Tx, Ty, Tz)
            end if

            if (hasUt_NS) then
                allocate( e % ut_NSout(1:NSTAT,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(1, e % Nsol, e % ut_NS, e % Nout, e % ut_NSout, Tx, Ty, Tz)
            end if

            if (hasMu_NS) then
                allocate( e % mu_NSout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(1, e % Nsol, e % mu_NS, e % Nout, e % mu_NSout, Tx, Ty, Tz)            
            end if

            if (hasWallY) then
                allocate( e % wallYout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(1, e % Nsol, e % WallY, e % Nout, e % wallYout, Tx, Ty, Tz)            
            end if

            if (hasMu_sgs) then
                allocate( e % mu_sgsout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(1, e % Nsol, e % mu_sgs, e % Nout, e % mu_sgsout, Tx, Ty, Tz)            
            end if

            if (hasSource) then
                allocate( e % sourceout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
                call prolongSolutionToGaussPoints(1, e % Nsol, e % source, e % Nout, e % sourceout, Tx, Ty, Tz)            
            end if

            allocate( e % QDot_out(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            if ( hasTimeDeriv ) then
                call prolongSolutionToGaussPoints(NVARS, e % Nsol, e % QDot, e % Nout, e % QDot_out, Tx, Ty, Tz)
            end if

         end if

      end subroutine ProjectStorageGaussPoints_FixedOrder
!
!////////////////////////////////////////////////////////////////////////////
!
!     Homogeneous procedures
!     ----------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_Homogeneous(meshName, solutionName, Nout, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
         real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
         integer                                    :: i
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Set homogeneous nodes
!        ---------------------
         xi   = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(1),i=0,Nout(1)) /), (/ Nout(1)+1 /) )
         eta  = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(2),i=0,Nout(2)) /), (/ Nout(2)+1 /) )
         zeta = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(3),i=0,Nout(3)) /), (/ Nout(3)+1 /) )
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)      ! TODO: check why it was Nmesh(1) 
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)     ! TODO: check why it was Nmesh(1)

!
!           Construct interpolation matrices for the solution
!           -------------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), eta)        ! TODO: check why it was Nout(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), zeta)       ! TODO: check why it was Nout(1)
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                     Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                     Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                     Tset(e % Nout(3), e % Nsol(3)) % T, &
                                                     mesh % hasGradients, mesh % isStatistics, mesh % hasTimeDeriv)


            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         call getOutputVariables()
         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM) then
            call WriteSingleFluidZoneToTecplot(fid,mesh)
         else
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )

               call WriteElementToTecplot(fid, e, mesh % refs, &
                                          mesh % hasGradients, mesh % isStatistics, mesh % hasSensor)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if

!
!        Close the file
!        --------------
         close(fid)

      end subroutine Solution2Plt_Homogeneous

      subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol, hasGradients, hasStats, hasTimeDeriv)
         use Storage
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
         real(kind=RP),       intent(in)  :: TxSol(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: TySol(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: TzSol(0:e % Nout(3), 0:e % Nsol(3))
         logical,             intent(in)  :: hasGradients
         logical,             intent(in)  :: hasStats
         logical,             intent(in)  :: hasTimeDeriv
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l, m, n
!
!        Project mesh
!        ------------         
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

!
!        Project the solution
!        --------------------
         allocate( e % Qout(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % Qout = 0.0_RP

         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % Qout(:,i,j,k) = e % Qout(:,i,j,k) + e % Q(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do


         allocate( e % QDot_out(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % QDot_out = 0.0_RP
         if ( hasTimeDeriv ) then

             do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
                do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                   e % QDot_out(:,i,j,k) = e % QDot_out(:,i,j,k) + e % QDot(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
                end do            ; end do            ; end do
             end do            ; end do            ; end do
         end if

         if ( hasGradients ) then
            allocate( e % U_xout(1:NGRADVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)))
            e % U_xout = 0.0_RP
   
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_xout(:,i,j,k) = e % U_xout(:,i,j,k) + e % U_x(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

          allocate( e % U_yout(1:NGRADVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % U_yout = 0.0_RP
   
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_yout(:,i,j,k) = e % U_yout(:,i,j,k) + e % U_y(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

          allocate( e % U_zout(1:NGRADVARS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % U_zout = 0.0_RP
   
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % U_zout(:,i,j,k) = e % U_zout(:,i,j,k) + e % U_z(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do

         end if

         if (hasStats) then
            allocate( e % statsout(1:NSTAT,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % statsout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % statsout(:,i,j,k) = e % statsout(:,i,j,k) + e % stats(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do
         end if
         
         if (hasUt_NS) then
            allocate( e % ut_NSout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % ut_NSout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % ut_NSout(:,i,j,k) = e % ut_NSout(:,i,j,k) + e % ut_NS(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do
         end if
         
         if (hasMu_NS) then
            allocate( e % mu_NSout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % mu_NSout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % mu_NSout(:,i,j,k) = e % mu_NSout(:,i,j,k) + e % mu_NS(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do  
         end if
         
         if (hasWallY) then
            allocate( e % wallYout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % wallYout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % wallYout(:,i,j,k) = e % wallYout(:,i,j,k) + e % WallY(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do  
         end if
         
         if (hasMu_sgs) then
            allocate( e % mu_sgsout(1,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % mu_sgsout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % mu_sgsout(:,i,j,k) = e % mu_sgsout(:,i,j,k) + e % mu_sgs(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do  
         end if

         if (hasSource) then
            allocate( e % sourceout(1:NVARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % sourceout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % sourceout(:,i,j,k) = e % sourceout(:,i,j,k) + e % source(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do  
         end if

      end subroutine ProjectStorageHomogeneousPoints
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Write solution
!     --------------
!
!/////////////////////////////////////////////////////////////////////////////
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
         numOfPoints   = product(mesh % elements(1) % Nout + 1)
         if (mesh % isSurface) then
             numOfFElems   = product(mesh % elements(1) % Nout(1:2))
         else
             numOfFElems   = product(mesh % elements(1) % Nout)
         end if
         firstPoint(1) = 1
         do eID = 2, size(mesh % elements)
            associate ( e => mesh % elements(eID) )
            firstPoint(eID) = numOfPoints + 1
            numOfPoints = numOfPoints + product(e % Nout + 1)
            if (mesh % isSurface) then
                numOfFElems   = numOfFElems + product(e % Nout(1:2))
            else
                numOfFElems   = numOfFElems + product(e % Nout)
            end if
            ! numOfFElems = numOfFElems + product(e % Nout    )
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
            N = e % Nout
            allocate (e % outputVars(1:no_of_outputVariables, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call ComputeOutputVariables(no_of_outputVariables, outputVariableNames, e % Nout, e, e % outputVars, mesh % refs, &
                                        mesh % hasGradients, mesh % isStatistics, mesh % hasSensor)
            
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               write(fid,trim(formatout)) e % xOut(:,i,j,k), e % outputVars(:,i,j,k)
            end do                ; end do                ; end do
            end associate
         end do
         
!        Write the elems connectivity
!        ----------------------------
         ! surface mesh case
         if (mesh % isSurface) then
             do eID = 1, size(mesh % elements)
                associate ( e => mesh % elements(eID) )

                do j = 0, e % Nout(2) - 1 ; do i = 0, e % Nout(1) - 1
                   cornersFace =  [ ij2localDOF(i,j,e%Nout(1:2)), ij2localDOF(i+1,j,e%Nout(1:2)), ij2localDOF(i+1,j+1,e%Nout(1:2)), ij2localDOF(i,j+1,e%Nout(1:2)) ] + firstPoint(eID)
                   write(fid,*) cornersFace
                end do                  ; end do

                end associate
             end do
         ! normal elements case
         else
             do eID = 1, size(mesh % elements)
                associate ( e => mesh % elements(eID) )
                
                do k = 0, e % Nout(3) - 1 ; do j = 0, e % Nout(2) - 1 ; do i = 0, e % Nout(1) - 1
                   corners =  [ ijk2localDOF(i,j,k  ,e%Nout), ijk2localDOF(i+1,j,k  ,e%Nout), ijk2localDOF(i+1,j+1,k  ,e%Nout), ijk2localDOF(i,j+1,k  ,e%Nout), &
                                ijk2localDOF(i,j,k+1,e%Nout), ijk2localDOF(i+1,j,k+1,e%Nout), ijk2localDOF(i+1,j+1,k+1,e%Nout), ijk2localDOF(i,j+1,k+1,e%Nout)  ] + firstPoint(eID)
                   write(fid,*) corners
                end do                    ; end do                    ; end do
                
                end associate
             end do
         end if
      end subroutine WriteSingleFluidZoneToTecplot
!
!////////////////////////////////////////////////////////////////////////////
!
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
!
!//////////////////////////////////////////////////////////////////////////////
!
!     Writes a single boundary zone using the FE Tecplot format
!     -> This format is more efficiently read by paraview and tecplot.
!     ------------------------------------------------------
      subroutine WriteSingleBoundaryZoneToTecplot(fd,boundary, elements)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use OutputVariables
         use SolutionFile
         implicit none
         !-arguments-------------------------------------------
         integer         , intent(in) :: fd
         type(Boundary_t), intent(in) :: boundary
         type(Element_t) , intent(in) :: elements(:)
         !-local-variables-------------------------------------
         integer :: numOfPoints  ! Number of plot points
         integer :: numOfFElems  ! Number of FINITE elements
         integer :: fID, side
         integer :: corners(4)
         integer :: i,j,k
         integer :: N(3)
         integer :: firstPoint(boundary % no_of_faces)
         integer :: Nf      (2,boundary % no_of_faces)
         character(len=LINE_LENGTH) :: formatout
         !-----------------------------------------------------
         
         formatout = getFormat()
         
         ! Count points and elements
         numOfPoints   = 0
         numOfFElems   = 0
         
         do fID = 1, boundary % no_of_faces
            associate (e => elements( boundary % elements(fID) ))
            side = boundary % elementSides(fID)
            
            select case (side)
               case(1,2) ; Nf(:,fID) = [e % Nout(1), e % Nout(3)]
               case(3,5) ; Nf(:,fID) = [e % Nout(1), e % Nout(2)]
               case(4,6) ; Nf(:,fID) = [e % Nout(2), e % Nout(3)]
            end select
            
            firstPoint(fID) = numOfPoints + 1
            numOfPoints     = numOfPoints + product(Nf(:,fID)+1)
            numOfFElems     = numOfFElems + product(Nf(:,fID)  )
            end associate
         end do

         ! don't write if boundary doesn't have elements associated, happens for periodic conditions
         if (numOfFElems .eq. 0) return
         
         write(fd,'(A,I0,A,I0,A,A,A)') "ZONE N=", numOfPoints,", E=", numOfFElems, &
                                                  ',ET=QUADRILATERAL, F=FEPOINT, T="boundary_', trim(boundary % Name), '"'
                  
!        Write the points
!        ----------------
         do fID=1, boundary % no_of_faces
            
            associate (e => elements( boundary % elements(fID) ))
            side = boundary % elementSides(fID)
            N = e % Nout
            select case (side)
            
               case(1)
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,0,k), e % outputVars(:,i,0,k)
                  end do                ; end do
                  
               case(2)
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,e % Nout(2),k), e % outputVars(:,i,e % Nout(2),k)
                  end do                ; end do
               
               case(3)
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,0), e % outputVars(:,i,j,0)
                  end do                ; end do
                  
               case(4)
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,e % Nout(1),j,k), e % outputVars(:,e % Nout(1),j,k)
                  end do                ; end do
                  
               case(5)
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,e % Nout(3)), e % outputVars(:,i,j,e % Nout(3))
                  end do                ; end do
                  
               case(6)
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,0,j,k), e % outputVars(:,0,j,k)
                  end do                ; end do
                  
            end select
            
            end associate
         end do
         
!        Write the elems connectivity
!        ----------------------------
         do fID = 1, boundary % no_of_faces
            
            do j = 0, Nf(2,fID) - 1 ; do i = 0, Nf(1,fID) - 1
               corners =  [ ij2localDOF(i,j,Nf(:,fID)), ij2localDOF(i+1,j,Nf(:,fID)), ij2localDOF(i+1,j+1,Nf(:,fID)), ij2localDOF(i,j+1,Nf(:,fID)) ] + firstPoint(fID)
               write(fd,*) corners
            end do                  ; end do
            
         end do
         
      end subroutine WriteSingleBoundaryZoneToTecplot
!
!////////////////////////////////////////////////////////////////////////////
!
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
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine WriteElementToTecplot(fid,e,refs, hasGradients, hasStats, hasSensor)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use OutputVariables
         use SolutionFile
         implicit none
         integer,            intent(in)    :: fid
         type(Element_t),    intent(inout) :: e 
         real(kind=RP),      intent(in)    :: refs(NO_OF_SAVED_REFS)
         logical,            intent(in)    :: hasGradients
         logical,            intent(in)    :: hasStats
         logical,            intent(in)    :: hasSensor
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i,j,k,var
         character(len=LINE_LENGTH) :: formatout
!
!        Get output variables
!        --------------------
         
         allocate (e % outputVars(1:no_of_outputVariables, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         call ComputeOutputVariables(no_of_outputVariables, outputVariableNames, e % Nout, e, e % outputVars, refs, hasGradients, hasStats, hasSensor)
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                            ", K=",e % Nout(3)+1,", F=POINT"

         formatout = getFormat()

         do k = 0, e % Nout(3)   ; do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
            write(fid,trim(formatout)) e % xOut(:,i,j,k), e % outputVars(:,i,j,k)
         end do               ; end do                ; end do

      end subroutine WriteElementToTecplot
      
      subroutine WriteBoundaryToTecplot(fd,boundary, elements)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use OutputVariables
         use SolutionFile
         implicit none
         !-arguments-------------------------------------------
         integer         , intent(in) :: fd
         type(Boundary_t), intent(in) :: boundary
         type(Element_t) , intent(in) :: elements(:)
         !-local-variables-------------------------------------
         integer :: fID, side
         integer :: i,j,k
         character(len=LINE_LENGTH) :: formatout
         !-----------------------------------------------------
         
         formatout = getFormat()
         
         do fID=1, boundary % no_of_faces
            
            associate (e => elements( boundary % elements(fID) ))
            side = boundary % elementSides(fID)
            
            select case (side)
            
               case(1)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(3)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,0,k), e % outputVars(:,i,0,k)
                  end do                ; end do
                  
               case(2)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(3)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,e % Nout(2),k), e % outputVars(:,i,e % Nout(2),k)
                  end do                ; end do
               
               case(3)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,0), e % outputVars(:,i,j,0)
                  end do                ; end do
                  
               case(4)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(2)+1,", J=",e % Nout(3)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,e % Nout(1),j,k), e % outputVars(:,e % Nout(1),j,k)
                  end do                ; end do
                  
               case(5)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,e % Nout(3)), e % outputVars(:,i,j,e % Nout(3))
                  end do                ; end do
                  
               case(6)
                  
                  write(fd,'(A,I0,A,I0,A,I0,A,A,I0,A)') "ZONE I=",e % Nout(2)+1,", J=",e % Nout(3)+1, &
                                                  ", K=",1,', F=POINT, T="boundary_', trim(boundary % Name), fID, '"'
                  
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,0,j,k), e % outputVars(:,0,j,k)
                  end do                ; end do
                  
            end select
            
            end associate
         end do
         
      end subroutine WriteBoundaryToTecplot
      

      character(len=LINE_LENGTH) function getFormat()
         use OutputVariables
         implicit none

         getFormat = ""

         write(getFormat,'(A,I0,A,A)') "(",3+no_of_outputVariables,PRECISION_FORMAT,")"

      end function getFormat

end module Solution2PltModule
