#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_VegOverride
!=======================================================================
! Per-patch aerenchyma parameter overrides for wetland (patchtype==2).
!
! Written by get_wetland_veg_proxy in MOD_MT_Methane_BgcLink based
! on 5-zone climate classification (tropical reed/papyrus, tropical
! swamp, temperate marsh, boreal fen, Sphagnum bog).
!
! Read by methane_aere / SiteOxAere in MOD_MT_Methane_Physics via
! getters that fall through to DEF_METHANE%* defaults when the
! per-patch override is inactive.
!
! Module exists separately from BgcLink so Physics can read these
! arrays without a Physics->BgcLink dependency (BgcLink compiles
! after Physics in the Makefile).
!
!=======================================================================

   USE MOD_Precision
   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), allocatable, public :: wetland_aere_poros  (:)   ! tiller porosity (-)
   real(r8), allocatable, public :: wetland_aere_radius (:)   ! tiller radius (m)
   real(r8), allocatable, public :: wetland_aere_tillerC(:)   ! gC per tiller
   real(r8), allocatable, public :: wetland_aere_scale  (:)   ! scale_factor_aere multiplier (-)
   logical,  allocatable, public :: wetland_aere_active (:)   ! .true. = use override for patch

   PUBLIC :: allocate_wetland_aere_overrides
   PUBLIC :: deallocate_wetland_aere_overrides
   PUBLIC :: get_aere_poros, get_aere_radius, get_aere_tillerC, get_aere_scale

CONTAINS

   SUBROUTINE allocate_wetland_aere_overrides(numpatch)
      integer, intent(in) :: numpatch
      IF (allocated(wetland_aere_poros)) RETURN
      ! Keep zero-length arrays allocated on ranks with numpatch==0; downstream
      ! getters guard by bounds, and collective setup expects allocation state
      ! to be consistent across MPI ranks.
      allocate(wetland_aere_poros  (numpatch))
      allocate(wetland_aere_radius (numpatch))
      allocate(wetland_aere_tillerC(numpatch))
      allocate(wetland_aere_scale  (numpatch))
      allocate(wetland_aere_active (numpatch))
      wetland_aere_poros  (:) = 0._r8
      wetland_aere_radius (:) = 0._r8
      wetland_aere_tillerC(:) = 0._r8
      wetland_aere_scale  (:) = 0._r8
      wetland_aere_active (:) = .false.
   END SUBROUTINE allocate_wetland_aere_overrides

   SUBROUTINE deallocate_wetland_aere_overrides()
      IF (allocated(wetland_aere_poros))   deallocate(wetland_aere_poros)
      IF (allocated(wetland_aere_radius))  deallocate(wetland_aere_radius)
      IF (allocated(wetland_aere_tillerC)) deallocate(wetland_aere_tillerC)
      IF (allocated(wetland_aere_scale))   deallocate(wetland_aere_scale)
      IF (allocated(wetland_aere_active))  deallocate(wetland_aere_active)
   END SUBROUTINE deallocate_wetland_aere_overrides

   real(r8) FUNCTION get_aere_poros(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_poros = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_poros = wetland_aere_poros(ipatch)
   END FUNCTION get_aere_poros

   real(r8) FUNCTION get_aere_radius(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_radius = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_radius = wetland_aere_radius(ipatch)
   END FUNCTION get_aere_radius

   real(r8) FUNCTION get_aere_tillerC(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_tillerC = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_tillerC = wetland_aere_tillerC(ipatch)
   END FUNCTION get_aere_tillerC

   real(r8) FUNCTION get_aere_scale(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_scale = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_scale = wetland_aere_scale(ipatch)
   END FUNCTION get_aere_scale

END MODULE MOD_MT_Methane_VegOverride
#endif
