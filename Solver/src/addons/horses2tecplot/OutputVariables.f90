#include "Includes.h"
module OutputVariables
!
!//////////////////////////////////////////////////////////////////////////////////
!
!        This module selects user-defined output variables for the
!     tecplot file. The user may add extra variables, this can be done by
!     following the steps:
!
!           * Adding a new ID for the variable NEWVAR_V (following the given order).
!           * Adding a new key for the variable, NEWVARKey.
!           * Adding the new key to the keys array.
!           * Adding the "select case" procedure that computes the output
!                 variable.
!
!//////////////////////////////////////////////////////////////////////////////////
!
   use SMConstants
   use PhysicsStorage
   use Headers
   use Storage, only: NVARS, hasMPIranks
   use PhysicsStorage_NSTPT, only: IC
   use PhysicsStorage_NSSA, only: IRHOTHETA
   use PhysicsStorage_iNS, only: INSRHO, INSRHOU, INSRHOV, INSRHOW, INSP
   use PhysicsStorage_MU, only: IMC, IMSQRHOU, IMSQRHOV, IMSQRHOW, IMP
   use PhysicsStorage_CAA, only: ICAARHO, ICAAU, ICAAV, ICAAW, ICAAP
   use PhysicsStorage_CAA, only: NCONSB, IBRHO, IBU, IBV, IBW, IBP, IBA2

   private
   public   no_of_outputVariables, preliminarNoOfVariables, askedVariables, getNoOfCommas
   public   outputVariableNames, preliminarVariables, variableNames
   public   getOutputVariables, ComputeOutputVariables, getOutputVariablesLabel
   public   getOutputVariablesList, outputVariablesForVariable, OutputVariablesForPreliminarVariable
   public   outScale, hasVariablesFlag, Lreference
				 
   integer, parameter   :: STR_VAR_LEN = 16
!
!  ***************************
!  Variables without gradients
!  ***************************
!
   enum, bind(C)
      enumerator :: Q_V = 1, QDot_V, RHO_V, U_V, V_V, W_V
      enumerator :: P_V, P0_V , RHODOT_V, RHOUDOT_V, RHOVDOT_V, RHOWDOT_V
      enumerator :: RHOEDOT_V, CDOT_V, T_V, Mach_V, S_V, Vabs_V
      enumerator :: Vvec_V, Ht_V, RHOU_V, RHOV_V
      enumerator :: RHOW_V, RHOE_V, C_V, Cp_V
      enumerator :: SQRHOU_V, SQRHOV_V, SQRHOW_V
      enumerator :: Nxi_V, Neta_V, Nzeta_V, Nav_V, N_V
      enumerator :: Xi_V, Eta_V, Zeta_V, ThreeAxes_V, Axes_V, eID_V
      enumerator :: LAMBvec_V, LAMBx_V, LAMBy_V, LAMBz_V
      enumerator :: QBase_V, RHOBase_V, UBase_V, VBase_V, WBase_V, PBase_V, A2Base_V
      enumerator :: MPIRANK_V
      enumerator :: GRADV_V, UX_V, VX_V, WX_V
      enumerator :: UY_V, VY_V, WY_V, UZ_V, VZ_V, WZ_V
      enumerator :: CX_V, CY_V, CZ_V
      enumerator :: OMEGA_V, OMEGAX_V, OMEGAY_V, OMEGAZ_V
      enumerator :: OMEGAABS_V, QCRIT_V, DIV_V
      enumerator :: Vvec_Vmean, U_Vmean, V_Vmean, W_Vmean, ReST 
      enumerator :: ReSTxx, ReSTxy, ReSTxz, ReSTyy, ReSTyz, ReSTzz
      enumerator :: Vfvec_Vrms, Uf_Vrms, Vf_Vrms, Wf_Vrms
      enumerator :: U_TAU_V, WallY_V, Tauw_V, MU, YPLUS, Cf_V, MUTMINF
      enumerator :: MU_sgs_V, Source_V, RHO_S_V, RHOU_S_V, RHOV_S_V
      enumerator :: RHOW_S_V, RHOE_S_V, C_S_V, SENSOR_V
      enumerator :: LASTVARIABLE
   end enum

   integer, parameter   :: NO_OF_VARIABLES = LASTVARIABLE-1
   integer, parameter   :: NO_OF_INVISCID_VARIABLES = MPIRANK_V

   character(len=STR_VAR_LEN), parameter  :: QKey          = "Q"
   character(len=STR_VAR_LEN), parameter  :: QDOTKey       = "QDot"
   character(len=STR_VAR_LEN), parameter  :: RHOKey        = "rho"
   character(len=STR_VAR_LEN), parameter  :: UKey          = "u"
   character(len=STR_VAR_LEN), parameter  :: VKey          = "v"
   character(len=STR_VAR_LEN), parameter  :: WKey          = "w"
   character(len=STR_VAR_LEN), parameter  :: PKey          = "p"
   character(len=STR_VAR_LEN), parameter  :: P0Key         = "pt"
   character(len=STR_VAR_LEN), parameter  :: RHODOTKey     = "rhoDot"
   character(len=STR_VAR_LEN), parameter  :: RHOUDOTKey    = "rhouDot"
   character(len=STR_VAR_LEN), parameter  :: RHOVDOTKey    = "rhovDot"
   character(len=STR_VAR_LEN), parameter  :: RHOWDOTKey    = "rhowDot"
   character(len=STR_VAR_LEN), parameter  :: RHOEDOTKey    = "rhoeDot"
   character(len=STR_VAR_LEN), parameter  :: CDOTKey       = "cDot"
   character(len=STR_VAR_LEN), parameter  :: TKey          = "T"
   character(len=STR_VAR_LEN), parameter  :: MachKey       = "Mach"
   character(len=STR_VAR_LEN), parameter  :: SKey          = "s"
   character(len=STR_VAR_LEN), parameter  :: VabsKey       = "Vabs"
   character(len=STR_VAR_LEN), parameter  :: VvecKey       = "V"
   character(len=STR_VAR_LEN), parameter  :: HtKey         = "Ht"
   character(len=STR_VAR_LEN), parameter  :: RHOUKey       = "rhou"
   character(len=STR_VAR_LEN), parameter  :: RHOVKey       = "rhov"
   character(len=STR_VAR_LEN), parameter  :: RHOWKey       = "rhow"
   character(len=STR_VAR_LEN), parameter  :: RHOEKey       = "rhoe"
   character(len=STR_VAR_LEN), parameter  :: cKey          = "c"
   character(len=STR_VAR_LEN), parameter  :: CpKey         = "Cp"
   character(len=STR_VAR_LEN), parameter  :: SQRHOUKey     = "sqrhou"
   character(len=STR_VAR_LEN), parameter  :: SQRHOVKey     = "sqrhov"
   character(len=STR_VAR_LEN), parameter  :: SQRHOWKey     = "sqrhow"
   character(len=STR_VAR_LEN), parameter  :: NxiKey        = "Nxi"
   character(len=STR_VAR_LEN), parameter  :: NetaKey       = "Neta"
   character(len=STR_VAR_LEN), parameter  :: NzetaKey      = "Nzeta"
   character(len=STR_VAR_LEN), parameter  :: NavKey        = "Nav"
   character(len=STR_VAR_LEN), parameter  :: NKey          = "N"
   character(len=STR_VAR_LEN), parameter  :: XiKey         = "Ax_Xi"
   character(len=STR_VAR_LEN), parameter  :: EtaKey        = "Ax_Eta"
   character(len=STR_VAR_LEN), parameter  :: ZetaKey       = "Ax_Zeta"
   character(len=STR_VAR_LEN), parameter  :: ThreeAxesKey  = "ThreeAxes"
   character(len=STR_VAR_LEN), parameter  :: AxesKey       = "Axes"
   character(len=STR_VAR_LEN), parameter  :: eIDKey        = "eID"
   character(len=STR_VAR_LEN), parameter  :: LambKey       = "lamb"
   character(len=STR_VAR_LEN), parameter  :: LambxKey      = "lamb_x"
   character(len=STR_VAR_LEN), parameter  :: LambyKey      = "lamb_y"
   character(len=STR_VAR_LEN), parameter  :: LambzKey      = "lamb_z"
   character(len=STR_VAR_LEN), parameter  :: QBaseKey      = "Qbase"
   character(len=STR_VAR_LEN), parameter  :: RHOBaseKey    = "rhobase"
   character(len=STR_VAR_LEN), parameter  :: UBaseKey      = "ubase"
   character(len=STR_VAR_LEN), parameter  :: VBaseKey      = "vbase"
   character(len=STR_VAR_LEN), parameter  :: WBaseKey      = "wbase"
   character(len=STR_VAR_LEN), parameter  :: PBaseKey      = "pbase"
   character(len=STR_VAR_LEN), parameter  :: A2BaseKey     = "A2base"
   character(len=STR_VAR_LEN), parameter  :: mpiRankKey    = "mpi_rank"
   character(len=STR_VAR_LEN), parameter  :: gradVKey      = "gradV"
   character(len=STR_VAR_LEN), parameter  :: uxKey         = "u_x"
   character(len=STR_VAR_LEN), parameter  :: vxKey         = "v_x"
   character(len=STR_VAR_LEN), parameter  :: wxKey         = "w_x"
   character(len=STR_VAR_LEN), parameter  :: uyKey         = "u_y"
   character(len=STR_VAR_LEN), parameter  :: vyKey         = "v_y"
   character(len=STR_VAR_LEN), parameter  :: wyKey         = "w_y"
   character(len=STR_VAR_LEN), parameter  :: uzKey         = "u_z"
   character(len=STR_VAR_LEN), parameter  :: vzKey         = "v_z"
   character(len=STR_VAR_LEN), parameter  :: wzKey         = "w_z"
   character(len=STR_VAR_LEN), parameter  :: cxKey         = "c_x"
   character(len=STR_VAR_LEN), parameter  :: cyKey         = "c_y"
   character(len=STR_VAR_LEN), parameter  :: czKey         = "c_z"
   character(len=STR_VAR_LEN), parameter  :: omegaKey      = "omega"
   character(len=STR_VAR_LEN), parameter  :: omegaxKey     = "omega_x"
   character(len=STR_VAR_LEN), parameter  :: omegayKey     = "omega_y"
   character(len=STR_VAR_LEN), parameter  :: omegazKey     = "omega_z"
   character(len=STR_VAR_LEN), parameter  :: omegaAbsKey   = "omega_abs"
   character(len=STR_VAR_LEN), parameter  :: QCriterionKey = "Qcrit"
   character(len=STR_VAR_LEN), parameter  :: divVKey       = "divV"
   character(len=STR_VAR_LEN), parameter  :: VvecMeanKey   = "Vmean"
   character(len=STR_VAR_LEN), parameter  :: UMeanKey      = "umean"
   character(len=STR_VAR_LEN), parameter  :: VMeanKey      = "vmean"
   character(len=STR_VAR_LEN), parameter  :: WMeanKey      = "wmean"
   character(len=STR_VAR_LEN), parameter  :: ReSTKey       = "Sij"
   character(len=STR_VAR_LEN), parameter  :: ReSTxxKey     = "Sxx"
   character(len=STR_VAR_LEN), parameter  :: ReSTxyKey     = "Sxy"
   character(len=STR_VAR_LEN), parameter  :: ReSTxzKey     = "Sxz"
   character(len=STR_VAR_LEN), parameter  :: ReSTyyKey     = "Syy"
   character(len=STR_VAR_LEN), parameter  :: ReSTyzKey     = "Syz"
   character(len=STR_VAR_LEN), parameter  :: ReSTzzKey     = "Szz"
   character(len=STR_VAR_LEN), parameter  :: VfvecRmsKey    = "Vfrms"
   character(len=STR_VAR_LEN), parameter  :: UfRmsKey       = "ufrms"
   character(len=STR_VAR_LEN), parameter  :: VfRmsKey       = "vfrms"
   character(len=STR_VAR_LEN), parameter  :: WfRmsKey       = "wfrms"
   character(len=STR_VAR_LEN), parameter  :: UTAUKey       = "u_tau"
   character(len=STR_VAR_LEN), parameter  :: WallYKey      = "wall_distance"
   character(len=STR_VAR_LEN), parameter  :: TauwKey       = "wall_shear"
   character(len=STR_VAR_LEN), parameter  :: muKey         = "mu_ns"
   character(len=STR_VAR_LEN), parameter  :: yplusKey      = "yplus"
   character(len=STR_VAR_LEN), parameter  :: cfKey         = "Cf"
   character(len=STR_VAR_LEN), parameter  :: mutminfKey    = "mutminf"
   character(len=STR_VAR_LEN), parameter  :: muSGSKey      = "mu_sgs"
   character(len=STR_VAR_LEN), parameter  :: sourceKey     = "source"
   character(len=STR_VAR_LEN), parameter  :: RhoSKey       = "rho_source"
   character(len=STR_VAR_LEN), parameter  :: RhoUSKey      = "rhou_source"
   character(len=STR_VAR_LEN), parameter  :: RhoVSKey      = "rhov_source"
   character(len=STR_VAR_LEN), parameter  :: RhoWSKey      = "rhow_source"
   character(len=STR_VAR_LEN), parameter  :: RhoESKey      = "rhoe_source"
   character(len=STR_VAR_LEN), parameter  :: CSKey         = "c_source"
   character(len=STR_VAR_LEN), parameter  :: sensorKey     = "sensor"

   character(len=STR_VAR_LEN), dimension(NO_OF_VARIABLES), parameter  :: variableNames = (/ QKey,QDOTKey, RHOKey, UKey, VKey, WKey, &
                                                                            PKey, P0Key, RHODOTKey, RHOUDOTKey, RHOVDOTKey, RHOWDOTKey, RHOEDOTKey, &
                                                                            CDOTKey, TKey, MachKey, SKey, VabsKey, &
                                                                            VvecKey, HtKey, RHOUKey, RHOVKey, RHOWKey, &
                                                                            RHOEKey, cKey, CpKey, &
                                                                            SQRHOUKey, SQRHOVKey, SQRHOWKey, &
                                                                            NxiKey, NetaKey, NzetaKey, NavKey, NKey, &
                                                                            XiKey, EtaKey, ZetaKey, ThreeAxesKey, AxesKey, eIDKey, &
                                                                            LambKey, LambxKey, LambyKey, LambzKey, &
                                                                            QBaseKey, RHOBaseKey, UBaseKey, VBaseKey, WBaseKey, PBaseKey, A2BaseKey, &
                                                                            mpiRankKey, &
                                                                            gradVKey, uxKey, vxKey, wxKey, &
                                                                            uyKey, vyKey, wyKey, uzKey, vzKey, wzKey, &
                                                                            cxKey, cyKey, czKey, &
                                                                            omegaKey, omegaxKey, omegayKey, omegazKey, &
                                                                            omegaAbsKey, QCriterionKey, divVKey, &
                                                                            VvecMeanKey, UMeanKey, VMeanKey, WMeanKey, ReSTKey, &
                                                                            ReSTxxKey, ReSTxyKey, ReSTxzKey, ReSTyyKey, ReSTyzKey, ReSTzzKey, &
                                                                            VfvecRmsKey, UfRmsKey, VfRmsKey, WfRmsKey, &
                                                                            UTAUKey, WallYKey, TauwKey, muKey, yplusKey, &
                                                                            cfKey, mutminfKey, muSGSKey, sourceKey, RhoSKey, &
                                                                            RhoUSKey, RhoVSKey, RhoWSKey, RhoESKey, CSKey, sensorKey /)
                                                                        
                                                                        
                                                               
   integer                          :: no_of_outputVariables, preliminarNoOfVariables
   integer, allocatable             :: outputVariableNames(:), preliminarVariables(:)
   logical                          :: outScale
   logical                          :: hasVariablesFlag
   character(len=LINE_LENGTH)       :: askedVariables   
   real(kind=RP)                    :: Lreference

   contains
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine getOutputVariables()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: pos, pos2, i
         character(len=STR_VAR_LEN)   :: inputVar
!
!        ***********************************************************
!              Read the variables. They are first loaded into
!           a "preliminar" variables, prior to introduce them in
!           the real outputVariables. This is because some of
!           the output variables lead to multiple variables (e.g.
!           Q, V or N).
!        ***********************************************************
!
         if ( .not. hasVariablesFlag ) then
!
!           Default: export Q
!           -------
            preliminarNoOfVariables = 1
            if (.not. allocated(preliminarVariables)) allocate( preliminarVariables(preliminarNoOfVariables) )
            preliminarVariables(1) = Q_V

         else
            pos = 0
!
!           Prepare to read the variable names
!           ----------------------------------
            preliminarNoOfVariables = getNoOfCommas(trim(askedVariables)) + 1
            if (.not. allocated(preliminarVariables)) allocate( preliminarVariables(preliminarNoOfVariables) )

            if ( preliminarNoOfVariables .eq. 1 ) then
               read(askedVariables(pos+1:len_trim(askedVariables)),*) inputVar
               preliminarVariables(1) = outputVariableForName(adjustl(trim(inputVar)))
            else
               do i = 1, preliminarNoOfVariables-1
                  pos2 = index(trim(askedVariables(pos+1:)),",") + pos
                  read(askedVariables(pos+1:pos2),*) inputVar
                  preliminarVariables(i) = outputVariableForName(adjustl(trim(inputVar)))
                  pos = pos2
               end do
            
               pos = index(trim(askedVariables),",",BACK=.true.)
               preliminarVariables(preliminarNoOfVariables) = outputVariableForName(TRIM(ADJUSTL(askedVariables(pos+1:))))
               
            end if
         end if
!
!        *******************************************************
!        Convert the preliminary variables into output variables
!        *******************************************************
!
         no_of_outputVariables = 0

         do i = 1, preliminarNoOfVariables
            no_of_outputVariables = no_of_outputVariables + outputVariablesForVariable(preliminarVariables(i))
         end do

         if (hasMPIranks)  no_of_outputVariables = no_of_outputVariables + 1

         safedeallocate( outputVariableNames) ; allocate( outputVariableNames(no_of_outputVariables) )

         pos = 1
         do i = 1, preliminarNoOfVariables
            pos2 = pos + outputVariablesForVariable(preliminarVariables(i)) - 1
            call outputVariablesForPreliminarVariable(preliminarVariables(i), outputVariableNames(pos:pos2) )

            pos = pos + outputVariablesForVariable(preliminarVariables(i))

         end do

         if (hasMPIranks) outputVariableNames(pos) = MPIRANK_V
!
!        *****************************
!        Describe the output variables
!        *****************************
!
         write(STD_OUT,'(/)')
         call Section_Header("Output variables")
         write(STD_OUT,'(/)')
         call Subsection_Header("Selected output variables")

         do i = 1, no_of_outputVariables
            write(STD_OUT,'(30X,A,A)') "* ",trim(variableNames(outputVariableNames(i)))

         end do

         if ( outScale ) then
            write(STD_OUT,'(30X,A,A)') "-> Variables are exported with dimensions."

         else
            write(STD_OUT,'(30X,A,A)') "-> Dimensionless mode."

         end if

      end subroutine getOutputVariables

      subroutine ComputeOutputVariables(noOutput, outputVarNames, N, e, output, refs, hasGradients, hasStats, hasSensor)
         use SolutionFile
         use Storage
         use StatisticsMonitor
         use Utilities, only: toLower
         implicit none
         integer, intent(in)          :: noOutput
         integer, intent(in)          :: outputVarNames(1:noOutput)
         integer, intent(in)          :: N(3)
         class(Element_t), intent(in) :: e
         real(kind=RP), intent(out)   :: output(1:noOutput,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(in)    :: refs(NO_OF_SAVED_REFS)
         logical,       intent(in)    :: hasGradients
         logical,       intent(in)    :: hasStats
         logical,       intent(in)    :: hasSensor
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: var, i, j, k
         real(kind=RP) :: Sym, Asym
         logical       :: hasAdditionalVariables

         hasAdditionalVariables = hasUt_NS .or. hasWallY .or. hasMu_NS .or. hasStats .or. hasGradients .or. hasSensor .or. hasMu_sgs .or. hasSource

         do var = 1, noOutput
            if ( hasAdditionalVariables .or. (outputVarNames(var) .le. NO_OF_INVISCID_VARIABLES ) ) then
               associate ( Q   => e % Qout, &
                           QDot=> e % QDot_out, &
                           U_x => e % U_xout, &
                           U_y => e % U_yout, &
                           U_z => e % U_zout, &
                           mu_NS => e % mu_NSout, &
                           wallY => e % wallYout, &
                           u_tau=> e % ut_NSout, &
                           mu_sgs => e % mu_sgsout, &
                           source => e % sourceout, &
                           stats => e % statsout)

               select case (outputVarNames(var))

               case(RHO_V)
                  select case (trim(flowEq))
                  case("caa")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(ICAARHO,i,j,k)
                     end do         ; end do         ; end do
                  case default
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IRHO,i,j,k)
                     end do         ; end do         ; end do
                  end select
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * output(var,:,:,:)

               case(U_V)
                  select case (trim(flowEq))
                  case("caa")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(ICAAU,i,j,k)
                     end do         ; end do         ; end do
                  case default
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
                     end do         ; end do         ; end do
                  end select
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(V_V)
                  select case (trim(flowEq))
                  case("caa")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(ICAAV,i,j,k)
                     end do         ; end do         ; end do
                  case default
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
                     end do         ; end do         ; end do
                  end select
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(W_V)
                  select case (trim(flowEq))
                  case("caa")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(ICAAW,i,j,k)
                     end do         ; end do         ; end do
                  case default
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
                     end do         ; end do         ; end do
                  end select
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(P_V)
                  select case (trim(flowEq))
                  case ("ns")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                          ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k))
                     end do         ; end do         ; end do
                     if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)

                  case ("ins")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(INSP,i,j,k)
                     end do         ; end do         ; end do
                     if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)

                  case ("mu")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IMP,i,j,k)
                     end do         ; end do         ; end do
                     if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
                  
                  case ("caa")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(ICAAP,i,j,k)
                     end do         ; end do         ; end do
                     if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
                  end select

              case(P0_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) &
                                                             /(Q(IRHO,i,j,k)**2)     ! Vabs**2
                     output(var,i,j,k) =  output(var,i,j,k) / ( refs(GAMMA_REF)*  &
                                        (refs(GAMMA_REF)-1.0_RP)*(Q(IRHOE,i,j,k) /Q(IRHO,i,j,k)-0.5_RP * &
                                            output(var,i,j,k)) )  ! Mach Number**2
                     output(var,i,j,k) = (1+((refs(GAMMA_REF)-1.0_RP)*0.5_RP)*output(var,i,j,k))** &
                                        (refs(GAMMA_REF)/(refs(GAMMA_REF)-1.0_RP))

                     output(var,i,j,k) = (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) /Q(IRHO,i,j,k)) &
                                       * output(var,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,i,j,k) = refs(RHO_REF) * (refs(V_REF)**2) &
                                    * output(var,i,j,k)

               case(RHODOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(IRHO,i,j,k)
                  end do         ; end do         ; end do

               case(RHOUDOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(IRHOU,i,j,k)
                  end do         ; end do         ; end do

               case(RHOVDOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(IRHOV,i,j,k)
                  end do         ; end do         ; end do

               case(RHOWDOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(IRHOW,i,j,k)
                  end do         ; end do         ; end do

               case(RHOEDOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(IRHOE,i,j,k)
                  end do         ; end do         ; end do

               case(CDOT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = QDot(6,i,j,k)
                  end do         ; end do         ; end do

               case(RHO_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(IRHO,i,j,k)
                  end do         ; end do         ; end do

               case(RHOU_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(IRHOU,i,j,k)
                  end do         ; end do         ; end do

               case(RHOV_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(IRHOV,i,j,k)
                  end do         ; end do         ; end do

               case(RHOW_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(IRHOW,i,j,k)
                  end do         ; end do         ; end do

               case(RHOE_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(IRHOE,i,j,k)
                  end do         ; end do         ; end do

               case(C_S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = source(6,i,j,k)
                  end do         ; end do         ; end do

               case(T_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = (refs(GAMMA_REF) - 1.0_RP) * refs(GAMMA_REF) * POW2(refs(MACH_REF)) / Q(IRHO,i,j,k) * (Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(T_REF) * output(var,:,:,:)

               case(MACH_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = (POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                     output(var,i,j,k) = sqrt( output(var,i,j,k) / ( refs(GAMMA_REF)*(refs(GAMMA_REF)-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * output(var,i,j,k)) ) )
                  end do         ; end do         ; end do

               case(S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = refs(GAMMA_REF) * (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP * &
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k)) / (Q(IRHO,i,j,k)**refs(GAMMA_REF))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF)*POW2(refs(V_REF))/refs(RHO_REF)**refs(GAMMA_REF) * output(var,:,:,:)

               case(Vabs_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(Ht_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = refs(GAMMA_REF)*Q(IRHOE,i,j,k) - 0.5_RP*(refs(GAMMA_REF)-1.0_RP)*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)

               case(RHOU_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOU,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)

               case(RHOV_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOV,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)

               case(RHOW_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOW,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)

               case(RHOE_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOE,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)

               case(SQRHOU_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IMSQRHOU,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = sqrt(refs(RHO_REF)) * refs(V_REF) * output(var,:,:,:)
               
               case(SQRHOV_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IMSQRHOV,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = sqrt(refs(RHO_REF)) * refs(V_REF) * output(var,:,:,:)
               
               case(SQRHOW_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IMSQRHOW,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = sqrt(refs(RHO_REF)) * refs(V_REF) * output(var,:,:,:)

               case(LAMBx_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(1,i,j,k)
                  end do         ; end do         ; end do
               
               case(LAMBy_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(2,i,j,k)
                  end do         ; end do         ; end do
               
               case(LAMBz_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(3,i,j,k)
                  end do         ; end do         ; end do
               
               case(RHOBase_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * output(var,:,:,:)

               case(UBase_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBU,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(VBase_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBV,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(WBase_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBW,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(PBase_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBP,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(A2Base_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IBA2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(Nxi_V)
                  output(var,:,:,:) = e % Nsol(1)

               case(Neta_V)
                  output(var,:,:,:) = e % Nsol(2)

               case(Nzeta_V)
                  output(var,:,:,:) = e % Nsol(3)

               case(Nav_V)
                  output(var,:,:,:) = sum(e % Nsol)/real(NDIM)

               case(Xi_V)
                  output(var,:,:,:) = 0
                  output(var,:,0,0) = 3
                  output(var,0,:,0) = 0
                  output(var,0,0,:) = 0

               case(Eta_V)
                  output(var,:,:,:) = 0
                  output(var,:,0,0) = 0
                  output(var,0,:,0) = 3
                  output(var,0,0,:) = 0

               case(Zeta_V)
                  output(var,:,:,:) = 0
                  output(var,:,0,0) = 0
                  output(var,0,:,0) = 0
                  output(var,0,0,:) = 3

               case(ThreeAxes_V)
                  output(var,:,:,:) = 0
                  output(var,:,0,0) = 3
                  output(var,0,:,0) = 1.5
                  output(var,0,0,:) = 1.5

               case(eID_V)
                  output(var,:,:,:) = e % eID

               case(MPIRANK_V)
                  output(var,:,:,:) = e % mpi_rank
!
!              ******************
!              Statistics variables   
!              ******************
!
               case(U_Vmean)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(U,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(V_Vmean)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(V,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(W_Vmean)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(W,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(ReSTxx)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(UU,i,j,k) - POW2(stats(U,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(ReSTxy)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(UV,i,j,k) - stats(U,i,j,k)*stats(V,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(ReSTxz)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(UW,i,j,k) - stats(U,i,j,k)*stats(W,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(ReSTyy)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(VV,i,j,k) - POW2(stats(V,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(ReSTyz)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(VW,i,j,k) - stats(V,i,j,k)*stats(W,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(ReSTzz)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = stats(WW,i,j,k) - POW2(stats(W,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(Uf_Vrms)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(stats(UU,i,j,k) - POW2(stats(U,i,j,k)))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(Vf_Vrms)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(stats(VV,i,j,k) - POW2(stats(V,i,j,k)))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

               case(Wf_Vrms)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(stats(WW,i,j,k) - POW2(stats(W,i,j,k)))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)

!
!              ******************
!              Gradient variables
!              ******************
!
               case(UX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(VX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(WX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(UY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(VY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(WY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(UZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(VZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)
               
               case(WZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(OMEGAX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_y(3,i,j,k) - U_z(2,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(OMEGAY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_z(1,i,j,k) - U_x(3,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(OMEGAZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_x(2,i,j,k) - U_y(1,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(OMEGAABS_V)

                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(  POW2( U_y(3,i,j,k) - U_z(2,i,j,k) ) &
                                              + POW2( U_z(1,i,j,k) - U_x(3,i,j,k) ) &
                                              + POW2( U_x(2,i,j,k) - U_y(1,i,j,k) ) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(DIV_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(1,i,j,k) + U_y(2,i,j,k) + U_z(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) / Lreference * output(var,:,:,:)

               case(QCRIT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     Sym =   POW2( U_x(1,i,j,k) ) + POW2( U_y(2,i,j,k) ) + POW2( U_z(3,i,j,k) )  &
                           + 2.0_RP *( POW2( 0.5_RP * (U_x(2,i,j,k) + U_y(1,i,j,k)) ) +          &
                                       POW2( 0.5_RP * (U_x(3,i,j,k) + U_z(1,i,j,k)) ) +          &
                                       POW2( 0.5_RP * (U_y(3,i,j,k) + U_z(2,i,j,k)) ) )

                     Asym =   2.0_RP *( POW2( 0.5_RP * (U_x(2,i,j,k) - U_y(1,i,j,k)) ) +        &
                                        POW2( 0.5_RP * (U_x(3,i,j,k) - U_z(1,i,j,k)) ) +        &
                                        POW2( 0.5_RP * (U_y(3,i,j,k) - U_z(2,i,j,k)) ) )

                     output(var,i,j,k) = 0.5_RP*( Asym - Sym )
                  end do            ; end do            ; end do
!
!              ******************
!              Turbulent variables
!              ******************
!
               case(MU)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = mu_NS(1,i,j,k) !* refs(RE_REF)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) ! there is not reference state for the mu, it could be calculated if we have the Re reference

               case(U_TAU_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) =  abs(u_tau(1,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) * refs(V_REF)

               case(WallY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) =  wallY(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) * Lreference

               case(Tauw_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) =  Q(IRHO,i,j,k) * POW2(u_tau(1,i,j,k)) * sign(1.0_RP,u_tau(1,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)

               case(YPLUS)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) =  wallY(1,i,j,k) * abs(u_tau(1,i,j,k)) * Q(IRHO,i,j,k) / mu_NS(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) * Lreference

               case(MU_sgs_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = mu_sgs(1,i,j,k) !* refs(RE_REF)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) ! there is not reference state for the mu, it could be calculated if we have the Re reference

               case(Cf_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) =  Q(IRHO,i,j,k) * POW2(u_tau(1,i,j,k)) * sign(1.0_RP,u_tau(1,i,j,k)) * 2.0_RP !here we are assuming that u_inf and rho_inf are 1.0 as is common in a free stream
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) ! is non dimensional

               case(Cp_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     !here we are assuming that u_inf and rho_inf are 1.0 as is common in a free stream
                     output(var,i,j,k) = 2.0_RP * ( (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k)) - &
                                       1.0_RP / (refs(GAMMA_REF)*POW2(refs(MACH_REF))) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:) ! is non dimensional

               case(mutminf)

                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                  if (Q(6,i,j,k) .GE. 0.0_RP ) then
                        output(var,i,j,k) = Q(6,i,j,k) * ((Q(6,i,j,k)/(refs(RE_REF) *mu_NS(1,i,j,k)))**3.0_RP)/(((Q(6,i,j,k)/ (refs(RE_REF) * mu_NS(1,i,j,k)))**3.0_RP) + (7.1_RP)**3.0_RP)
                  else

                        output(var,i,j,k) = 0.0_RP

                  end if
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = output(var,:,:,:)


               case(C_V)
                  select case (trim(flowEq))
                  case ("nstpt")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IC,i,j,k)
                     end do         ; end do         ; end do
                  case ("ch")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(size(Q,1),i,j,k)
                     end do         ; end do         ; end do
                  case ("mu")
                     do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                        output(var,i,j,k) = Q(IMC,i,j,k)
                     end do         ; end do         ; end do
                  end select

               case(CX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(size(U_x,1),i,j,k)
                  end do         ; end do         ; end do

               case(CY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(size(U_y,1),i,j,k)
                  end do         ; end do         ; end do

               case(CZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(size(U_z,1),i,j,k)
                  end do         ; end do         ; end do
!
!              **************
!              Element sensor
!              **************
!
               case(SENSOR_V)
                  output(var,:,:,:) = e % sensor

               end select
               end associate

            else
               output(var,:,:,:) = 0.0_RP

            end if
         end do

      end subroutine ComputeOutputVariables

      subroutine getOutputVariablesList(list)
         implicit none
         character(len=STRING_CONSTANT_LENGTH), allocatable, intent(out) :: list(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: iVar

         allocate(list(no_of_outputVariables))
         do iVar = 1, no_of_outputVariables
            list(iVar) = trim(variableNames(outputVariableNames(iVar)))
         end do

      end subroutine getOutputVariablesList

      character(len=1024) function getOutputVariablesLabel()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: iVar

         getOutputVariablesLabel = ""
         do iVar = 1, no_of_outputVariables
            write(getOutputVariablesLabel,'(A,A,A,A)') trim(getOutputVariablesLabel), ',"',trim(variableNames(outputVariableNames(iVar))),'"'
         end do

      end function getOutputVariablesLabel

      integer function outputVariablesForVariable(iVar)
!
!        ************************************************
!           This subroutine specifies if a variable
!           implies more than one variable, e.g.,
!           Q = [rho,rhou,rhov,rhow,rhoe], or
!           V = [u,v,w]
!        ************************************************
!
         implicit none
         integer,    intent(in)     :: iVar

         select case(iVar)

         case(Q_V)
            outputVariablesForVariable = NVARS

         case(QDot_V)
            outputVariablesForVariable = NVARS

         case(Vvec_V)
            outputVariablesForVariable = 3

         case(gradV_V)
            outputVariablesForVariable = 9

         case(omega_V)
            outputVariablesForVariable = 3

         case(N_V)
            outputVariablesForVariable = 4

         case(Axes_V)
            outputVariablesForVariable = 4

         case(Vvec_Vmean)
            outputVariablesForVariable = 3

         case(ReST)
            outputVariablesForVariable = 6

         case(Vfvec_Vrms)
            outputVariablesForVariable = 3

         case(Source_V)
            outputVariablesForVariable = NVARS
         
         case(LAMBvec_V)
            outputVariablesForVariable = 3
         
         case(QBase_V)
            outputVariablesForVariable = NCONSB

         case default
            outputVariablesForVariable = 1

         end select

      end function outputVariablesForVariable

      subroutine OutputVariablesForPreliminarVariable(iVar, output)
         use Storage               , only: flowEq
         implicit none
         integer, intent(in)     :: iVar
         integer, intent(out)    :: output(:)

         select case(iVar)

         case(Q_V)
            select case (trim(flowEq))
              case ("ns")
                  output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V/)
              case ("nssa")
                  output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V, C_V/)
              case ("nstpt")
                  output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V, C_V/)
              case ("ins")
                  output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, P_V/)
              case ("ch")
                  output = (/C_V/)
              case ("mu")
                  output = (/C_V, SQRHOU_V, SQRHOV_V, SQRHOW_V, P_V/)
              case ("caa")
                  output = (/RHO_V, U_V, V_V, W_V, P_V/)
          end select

         case(QDot_V)

            if ( NVARS .eq. 5 ) then
               output = (/RHODOT_V, RHOUDOT_V, RHOVDOT_V, RHOWDOT_V, RHOEDOT_V/)

            elseif ( NVARS .eq. 6 ) then
               output = (/RHODOT_V, RHOUDOT_V, RHOVDOT_V, RHOWDOT_V, RHOEDOT_V, CDOT_V/)

            elseif ( NVARS .eq. 1 ) then
               output = (/CDOT_V/)

            end if

         case(Source_V)

            if ( NVARS .eq. 5 ) then
               output = (/RHO_S_V, RHOU_S_V, RHOV_S_V, RHOW_S_V, RHOE_S_V/)

            elseif ( NVARS .eq. 6 ) then
               output = (/RHO_S_V, RHOU_S_V, RHOV_S_V, RHOW_S_V, RHOE_S_V, C_S_V/)

            elseif ( NVARS .eq. 1 ) then
               output = (/C_S_V/)

            end if

         case(Vvec_V)
            output = (/U_V, V_V, W_V/)

         case(gradV_V)
            output = (/UX_V, VX_V, WX_V, UY_V, VY_V, WY_V, UZ_V, VZ_V, WZ_V/)

         case(omega_V)
            output = (/OMEGAX_V, OMEGAY_V, OMEGAZ_V/)

         case(N_V)
            output = (/Nxi_V, Neta_V, Nzeta_V, Nav_V/)

         case(Axes_V)
            output = (/Xi_V, Eta_V, Zeta_V, ThreeAxes_V/)

         case(Vvec_Vmean)
            output = (/U_Vmean, V_Vmean, W_Vmean/)

         case(ReST)
            output = (/ReSTxx, ReSTxy, ReSTxz, ReSTyy, ReSTyz, ReSTzz/)

         case(Vfvec_Vrms)
            output = (/Uf_Vrms, Vf_Vrms, Wf_Vrms/)

         case(LAMBvec_V)
            output = (/LAMBx_V, LAMBy_V, LAMBz_V/)
         
         case(QBase_V)
            output = (/RHOBase_V, UBase_V, VBase_V, WBase_V, PBase_V, A2Base_V/)

         case default
            output = iVar

         end select

      end subroutine OutputVariablesForPreliminarVariable

      integer function outputVariableForName(variableName)
         implicit none
         character(len=*),    intent(in)     :: variableName

         do outputVariableForName = 1, NO_OF_VARIABLES
            if ( trim(variableNames(outputVariableForName)) .eq. trim(variableName) ) return
         end do
!
!        Return "-1" if not found
!        ------------------------
         outputVariableForName = -1

      end function outputVariableForName

      integer function getNoOfCommas(input_line)
         implicit none
         character(len=*), intent(in)     :: input_line
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i

         getNoOfCommas = 0
         do i = 1, len_trim(input_line)
            if ( input_line(i:i) .eq. "," ) getNoOfCommas = getNoOfCommas + 1
         end do

      end function getNoOfCommas
end module OutputVariables
