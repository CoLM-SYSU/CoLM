#include <define.h>

MODULE MOD_HtopReadin

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: HTOP_readin

CONTAINS

   SUBROUTINE HTOP_readin (dir_landdata, lc_year)

! ===========================================================
! Read in the canopy tree top height
! Revisions:
!  12/2025, Jiayi Xiang & Yuan: add crown structure data for tree PFTs.
! ===========================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_Namelist, only: DEF_RS_CROWN_STRUCTURE
   USE MOD_Vars_TimeInvariants
   USE MOD_LandPatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
   USE MOD_NetCDFVector
#ifdef SinglePoint
   USE MOD_SingleSrfdata
   USE MOD_NetCDFSerial
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of land cover data used
   character(len=256), intent(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: c
   character(len=256) :: landdir, lndname, cyear, fsrfdata
   integer :: i,j,t,p,ps,pe,m,n,npatch

   real(r8), allocatable :: htoplc  (:)
   real(r8), allocatable :: htoppft (:), hbotpft(:), cratio_pft(:)

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_landdata) // '/htop/' // trim(cyear)

! ==============USGS scheme===========================
#ifdef LULC_USGS

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

         ENDDO
      ENDIF

#endif
! ==============IGBP scheme===========================
#ifdef LULC_IGBP
#ifdef SinglePoint
      allocate (htoplc (numpatch))
      htoplc(:) = SITE_htop
#else
      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_read_vector (lndname, 'htop_patches', landpatch, htoplc)
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

            ! trees or woody savannas
            IF ( m<6 .or. m==8 ) THEN
               ! 01/06/2020, yuan: adjust htop reading
               ! 11/15/2021, yuan: adjust htop setting
               htop(npatch) = max(2., htoplc(npatch))
               hbot(npatch) = htoplc(npatch)*hbot0(m)/htop0(m)
               hbot(npatch) = max(1., hbot(npatch))
            ENDIF

         ENDDO
      ENDIF

      IF (allocated(htoplc))   deallocate ( htoplc )
#endif

! ==============PFT/PC scheme===========================
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifdef SinglePoint
      IF (numpft > 0) THEN
         allocate(htoppft(numpft))
         fsrfdata = trim(DEF_dir_landdata) // '/srfdata.nc'
         CAll ncio_read_serial (fsrfdata, 'canopy_height_pfts', htoppft)
         IF (DEF_RS_CROWN_STRUCTURE) THEN
            allocate(hbotpft(numpft))
            CAll ncio_read_serial (fsrfdata, 'canopy_bottom_height_pfts', hbotpft)
#ifdef LULC_IGBP_PC
            allocate(cratio_pft(numpft))
            CAll ncio_read_serial (fsrfdata, 'crown_aspect_ratio_pfts', cratio_pft)
#endif
         ENDIF
      ENDIF
#else
      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_read_vector (lndname, 'htop_pfts', landpft,   htoppft)
      IF (DEF_RS_CROWN_STRUCTURE) THEN
         lndname = trim(landdir)//'/hbot_pfts.nc'
         CALL ncio_read_vector (lndname, 'hbot_pfts', landpft, hbotpft)
#ifdef LULC_IGBP_PC
         lndname = trim(landdir)//'/cratio_pfts.nc'
         CALL ncio_read_vector (lndname, 'cratio_pfts', landpft, cratio_pft)
#endif
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            t = patchtype(npatch)        ! land cover type
                                         ! 0: soil, 1: urban, 2: wetland, 3: ice, 4: lake
            m = patchclass(npatch)       ! index of land cover type

            IF (t == 0) THEN             ! land cover type is soil, growing with vegetation
               ps = patch_pft_s(npatch)
               pe = patch_pft_e(npatch)

               DO p = ps, pe
                  n = pftclass(p)        ! pft index

                  htop_p(p) = htop0_p(n) ! for non-tree, use the default height
                  hbot_p(p) = hbot0_p(n)
#ifdef LULC_IGBP_PC
                  cratio_p(p) = 1.
#endif

                  ! for trees
                  ! 01/06/2020, yuan: adjust htop reading
                  ! 11/15/2021, yuan: adjust htop setting
                  ! 12/16/2025, xiangjy: adjust hbot and cratio reading
                  IF ( n>0 .and. n<9 ) THEN ! pft is tree
                     htop_p(p) = max(2., htoppft(p))

                     ! Diagnose default hbot from the PFT hbot/htop ratio.
                     hbot_p(p) = htoppft(p)*hbot0_p(n)/htop0_p(n)
                     hbot_p(p) = max(1., hbot_p(p))

                     IF (DEF_RS_CROWN_STRUCTURE) THEN
                        hbot_p(p) = hbotpft(p)
#ifdef LULC_IGBP_PC
                        cratio_p(p) = cratio_pft(p)
#endif
                     ENDIF

                  ENDIF
               ENDDO

               htop(npatch) = sum(htop_p(ps:pe)*pftfrac(ps:pe)) ! weighted average
               hbot(npatch) = sum(hbot_p(ps:pe)*pftfrac(ps:pe)) ! but in fact, only htop_p & hbot_p used in 3DCanopyRadiation
#ifdef LULC_IGBP_PC
               cratio(npatch) = sum(cratio_p(ps:pe)*pftfrac(ps:pe))
#endif
            ELSE ! if land cover type is not soil, use the default height
               htop(npatch) = htop0(m)
               hbot(npatch) = hbot0(m)
#ifdef LULC_IGBP_PC
               cratio(npatch) = 1.
#endif
            ENDIF

         ENDDO
      ENDIF

      IF (allocated(htoppft   )) deallocate(htoppft   )
      IF (allocated(hbotpft   )) deallocate(hbotpft   )
      IF (allocated(cratio_pft)) deallocate(cratio_pft)
#endif

   END SUBROUTINE HTOP_readin

END MODULE MOD_HtopReadin
