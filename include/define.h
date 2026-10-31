! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land subgrid type classification:
!    Select one of the following options.
#undef LULC_USGS
#undef LULC_IGBP
#undef LULC_IGBP_PFT
#define LULC_IGBP_PC

! 2.1 3D Urban model (put it temporarily here):
#undef URBAN_MODEL
!    Dependence:  only LULC_IGBP subgrid type for
!    single point URBAN_MODEL right now.
#if (defined URBAN_MODEL && defined SinglePoint)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 3. If defined, debug information is output.
#define CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#define RangeCheck
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. If defined, MPI parallelization is enabled.
#define USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef   Campbell_SOIL_MODEL
#define  vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define CatchLateralFlow
!    Conflicts :
#ifndef CATCHMENT
#undef CatchLateralFlow
#endif

! 6. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood
#if (defined SinglePoint)
#undef CaMa_Flood
#endif
#ifndef USEMPI
#undef CaMa_Flood
#endif

#define GridRiverLakeFlow
!    Conflicts :
!    The river network is built per worker, so it needs MPI (as CaMa_Flood).
#ifndef USEMPI
#undef GridRiverLakeFlow
#endif
#if (defined CATCHMENT || defined SinglePoint || defined CaMa_Flood)
#undef GridRiverLakeFlow
#endif

! 7. If defined, BGC model is used.
#undef BGC

!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#ifndef LULC_IGBP_PC
#undef BGC
#endif
#endif
! 7.1 If defined, CROP model is used
#undef CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif

! 8. If defined, open Land use and land cover change mode.
#undef LULCC

! 9. If defined, data assimilation is used.
#undef DataAssimilation
#if (defined DataAssimilation)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 10. Interface to AI model.
#undef USESplitAI

! 11. External lake models.
#undef EXTERNAL_LAKE

! 12. Hyperspectral scheme.
#undef HYPERSPECTRAL

! 13. If defined, material transport is compiled in (water isotopes,
!     solutes, river sediment, methane). Each material is switched on
!     in the namelist; with all of them off the results are unchanged.
#undef MaterialTransport
!    Conflicts : needs the van Genuchten-Mualem soil model.
#if (defined MaterialTransport && defined Campbell_SOIL_MODEL)
#error "MaterialTransport requires vanGenuchten_Mualem_SOIL_MODEL"
#endif
#if (defined MaterialTransport && defined CATCHMENT)
#error "MaterialTransport does not support the catchment (unit-catchment) grid (CATCHMENT); use GRID or UNSTRUCTURED (vector)."
#endif
!    Material transport is only meaningful with river-lake routing: enable it
!    automatically and drop the CaMa-Flood alternative (which does not route
!    land runoff materials). URBAN_MODEL is not supported together with it.
#if (defined MaterialTransport && defined URBAN_MODEL)
#error "MaterialTransport does not support URBAN_MODEL; disable URBAN_MODEL to build material transport."
#endif
#ifdef MaterialTransport
#ifdef CaMa_Flood
#undef CaMa_Flood
#endif
#if (!defined SinglePoint)
#ifndef GridRiverLakeFlow
#define GridRiverLakeFlow
#endif
#endif
#endif
