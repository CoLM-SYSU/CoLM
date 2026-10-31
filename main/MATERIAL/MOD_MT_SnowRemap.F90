#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_SnowRemap

   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_snow_combine_work, mt_snow_divide_work
   PUBLIC :: mt_snow_merge, mt_snow_copy, mt_snow_clear

   ! Local work follows the host's layer events; no water state is changed here.
   TYPE mt_snow_combine_work
      PRIVATE
      real(r8), allocatable :: ice(:), liquid(:), solid(:)
   CONTAINS
      PROCEDURE :: init => combine_init
      PROCEDURE :: accumulate => combine_accumulate
      PROCEDURE :: collapse => combine_collapse
   END TYPE mt_snow_combine_work

   TYPE mt_snow_divide_work
      PRIVATE
      real(r8), allocatable :: ice(:,:), liquid(:,:), solid(:,:)
      real(r8), allocatable :: excess_ice(:), excess_liquid(:), excess_solid(:)
      logical :: has_solid = .false.
   CONTAINS
      PROCEDURE :: init => divide_init
      PROCEDURE :: load_layer
      PROCEDURE :: store_layer
      PROCEDURE :: halve
      PROCEDURE :: copy_layer
      PROCEDURE :: save_excess
      PROCEDURE :: scale
      PROCEDURE :: add_excess
   END TYPE mt_snow_divide_work

CONTAINS

   SUBROUTINE mt_snow_merge(lb, dest, source, liquid, ice, solid)
      integer, intent(in) :: lb, dest, source
      real(r8), intent(inout), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)
      integer :: imt

      IF (.not. present(liquid) .or. .not. present(ice)) RETURN
      DO imt = 1, size(liquid,1)
         liquid(imt,dest) = liquid(imt,dest) + liquid(imt,source)
         ice(imt,dest) = ice(imt,dest) + ice(imt,source)
         IF (present(solid)) solid(imt,dest) = solid(imt,dest) + solid(imt,source)
      ENDDO
   END SUBROUTINE mt_snow_merge

   SUBROUTINE mt_snow_copy(lb, dest, source, liquid, ice, solid)
      integer, intent(in) :: lb, dest, source
      real(r8), intent(inout), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)
      integer :: imt

      IF (.not. present(liquid) .or. .not. present(ice)) RETURN
      DO imt = 1, size(liquid,1)
         liquid(imt,dest) = liquid(imt,source)
         ice(imt,dest) = ice(imt,source)
         IF (present(solid)) solid(imt,dest) = solid(imt,source)
      ENDDO
   END SUBROUTINE mt_snow_copy

   SUBROUTINE mt_snow_clear(lb, first, last, liquid, ice, solid)
      integer, intent(in) :: lb, first, last
      real(r8), intent(inout), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)

      IF (.not. present(liquid) .or. .not. present(ice)) RETURN
      IF (last < first) RETURN
      liquid(:,first:last) = 0._r8
      ice(:,first:last) = 0._r8
      IF (present(solid)) solid(:,first:last) = 0._r8
   END SUBROUTINE mt_snow_clear

   SUBROUTINE combine_init(this, liquid, ice)
      class(mt_snow_combine_work), intent(out) :: this
      real(r8), intent(in), optional :: liquid(:,:), ice(:,:)
      integer :: nmt

      IF (.not. present(liquid) .or. .not. present(ice)) RETURN
      nmt = size(liquid,1)
      allocate(this%ice(nmt), this%liquid(nmt), this%solid(nmt))
      this%ice = 0._r8; this%liquid = 0._r8; this%solid = 0._r8
   END SUBROUTINE combine_init

   SUBROUTINE combine_accumulate(this, lb, layer, liquid, ice, solid)
      class(mt_snow_combine_work), intent(inout) :: this
      integer, intent(in) :: lb, layer
      real(r8), intent(in), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)
      integer :: imt

      IF (.not. allocated(this%ice)) RETURN
      DO imt = 1, size(liquid,1)
         this%ice(imt) = this%ice(imt) + ice(imt,layer)
         this%liquid(imt) = this%liquid(imt) + liquid(imt,layer)
         IF (present(solid)) this%solid(imt) = this%solid(imt) + solid(imt,layer)
      ENDDO
   END SUBROUTINE combine_accumulate

   SUBROUTINE combine_collapse(this, lb, liquid, ice, solid, scv)
      class(mt_snow_combine_work), intent(in) :: this
      integer, intent(in) :: lb
      real(r8), intent(inout), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:), scv(:)
      integer :: imt

      IF (.not. allocated(this%ice)) RETURN
      IF (present(scv)) scv(:) = this%ice(:)
      DO imt = 1, size(liquid,1)
         liquid(imt,1) = liquid(imt,1) + this%liquid(imt)
         IF (present(solid)) solid(imt,1) = solid(imt,1) + this%solid(imt)
      ENDDO
      CALL mt_snow_clear(lb, lb, 0, liquid, ice, solid)
   END SUBROUTINE combine_collapse

   SUBROUTINE divide_init(this, liquid, ice, solid)
      class(mt_snow_divide_work), intent(out) :: this
      real(r8), intent(in), optional :: liquid(:,:), ice(:,:), solid(:,:)
      integer :: nmt

      IF (.not. present(liquid) .or. .not. present(ice)) RETURN
      nmt = size(liquid,1)
      this%has_solid = present(solid)
      allocate(this%ice(nmt,5), this%liquid(nmt,5), this%solid(nmt,5))
      allocate(this%excess_ice(nmt), this%excess_liquid(nmt), this%excess_solid(nmt))
      this%ice = 0._r8; this%liquid = 0._r8; this%solid = 0._r8
      this%excess_ice = 0._r8; this%excess_liquid = 0._r8; this%excess_solid = 0._r8
   END SUBROUTINE divide_init

   SUBROUTINE load_layer(this, lb, work_layer, host_layer, liquid, ice, solid)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: lb, work_layer, host_layer
      real(r8), intent(in), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)
      integer :: imt

      IF (.not. allocated(this%ice)) RETURN
      DO imt = 1, size(liquid,1)
         this%ice(imt,work_layer) = ice(imt,host_layer)
         this%liquid(imt,work_layer) = liquid(imt,host_layer)
         IF (present(solid)) this%solid(imt,work_layer) = solid(imt,host_layer)
      ENDDO
   END SUBROUTINE load_layer

   SUBROUTINE store_layer(this, lb, work_layer, host_layer, liquid, ice, solid)
      class(mt_snow_divide_work), intent(in) :: this
      integer, intent(in) :: lb, work_layer, host_layer
      real(r8), intent(inout), optional :: liquid(:,lb:), ice(:,lb:), solid(:,lb:)
      integer :: imt

      IF (.not. allocated(this%ice)) RETURN
      DO imt = 1, size(liquid,1)
         ice(imt,host_layer) = this%ice(imt,work_layer)
         liquid(imt,host_layer) = this%liquid(imt,work_layer)
         IF (present(solid)) solid(imt,host_layer) = this%solid(imt,work_layer)
      ENDDO
   END SUBROUTINE store_layer

   SUBROUTINE halve(this, layer)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: layer

      IF (.not. allocated(this%ice)) RETURN
      this%ice(:,layer) = this%ice(:,layer) / 2._r8
      this%liquid(:,layer) = this%liquid(:,layer) / 2._r8
      IF (this%has_solid) this%solid(:,layer) = this%solid(:,layer) / 2._r8
   END SUBROUTINE halve

   SUBROUTINE copy_layer(this, dest, source)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: dest, source

      IF (.not. allocated(this%ice)) RETURN
      this%ice(:,dest) = this%ice(:,source)
      this%liquid(:,dest) = this%liquid(:,source)
      IF (this%has_solid) this%solid(:,dest) = this%solid(:,source)
   END SUBROUTINE copy_layer

   SUBROUTINE save_excess(this, layer, fraction)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: layer
      real(r8), intent(in) :: fraction

      IF (.not. allocated(this%ice)) RETURN
      this%excess_ice(:) = fraction * this%ice(:,layer)
      this%excess_liquid(:) = fraction * this%liquid(:,layer)
      IF (this%has_solid) this%excess_solid(:) = fraction * this%solid(:,layer)
   END SUBROUTINE save_excess

   SUBROUTINE scale(this, layer, fraction)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: layer
      real(r8), intent(in) :: fraction

      IF (.not. allocated(this%ice)) RETURN
      this%ice(:,layer) = fraction * this%ice(:,layer)
      this%liquid(:,layer) = fraction * this%liquid(:,layer)
      IF (this%has_solid) this%solid(:,layer) = fraction * this%solid(:,layer)
   END SUBROUTINE scale

   SUBROUTINE add_excess(this, layer)
      class(mt_snow_divide_work), intent(inout) :: this
      integer, intent(in) :: layer

      IF (.not. allocated(this%ice)) RETURN
      this%ice(:,layer) = this%ice(:,layer) + this%excess_ice(:)
      this%liquid(:,layer) = this%liquid(:,layer) + this%excess_liquid(:)
      IF (this%has_solid) this%solid(:,layer) = this%solid(:,layer) + this%excess_solid(:)
   END SUBROUTINE add_excess

END MODULE MOD_MT_SnowRemap
#endif
