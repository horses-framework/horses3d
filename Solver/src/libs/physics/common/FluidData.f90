module FluidData
#if defined(SPALARTALMARAS)
   use FluidData_NSSA
#elif defined(NAVIERSTOKES) && (!(TRANSPORT))
   use FluidData_NS
#elif defined(NAVIERSTOKES) && (TRANSPORT)
   use FluidData_NSTPT
#elif defined(INCNS)
   use FluidData_iNS
#elif defined(MULTIPHASE)
   use FluidData_MU
#elif defined(ACOUSTIC)
   use FluidData_CAA
#endif
#if defined(CAHNHILLIARD)
   use FluidData_CH
#endif
   implicit none

end module FluidData
