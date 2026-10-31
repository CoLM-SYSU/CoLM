#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_Reservoir
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reservoir module in gridded mesh.
!
! Created by Shupeng Zhang, Oct 2025
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DataType

   integer :: totalnumresv
   integer :: numresv
   integer,  allocatable :: ucat2resv   (:)
   integer,  allocatable :: resv_global_id(:)
   integer,  allocatable :: resv_ucid     (:)
   type(pointer_int32_1d), allocatable :: resv_data_address (:)


   ! parameters
   integer,  allocatable :: dam_GRAND_ID  (:)  ! GRAND dam ID

   integer,  allocatable :: dam_build_year(:)  ! year in which the dam/barrier was built

   real(r8), allocatable :: volresv_total (:)  ! total reservoir volume      [m^3]
   real(r8), allocatable :: volresv_emerg (:)  ! emergency reservoir volume  [m^3]
   real(r8), allocatable :: volresv_adjust(:)  ! adjustment reservoir volume [m^3]
   real(r8), allocatable :: volresv_normal(:)  ! normal reservoir volume     [m^3]

   ! H22 static parameters, only allocated for method 2.
   real(r8), allocatable :: volresv_flood (:)  ! flood-control onset C [m^3]
   real(r8), allocatable :: resv_h22_coeff(:)  ! inflow release coefficient k [-]

   real(r8), allocatable :: qresv_flood   (:)  ! flood reservoir outflow      [m^3/s]
   real(r8), allocatable :: qresv_adjust  (:)  ! adjustment reservoir outflow [m^3/s]
   real(r8), allocatable :: qresv_normal  (:)  ! normal reservoir outflow     [m^3/s]

   ! fluxes
   real(r8), allocatable :: qresv_in      (:)  ! reservoir inflow  [m^3/s]
   real(r8), allocatable :: qresv_out     (:)  ! reservoir outflow [m^3/s]

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: reservoir_init
   PUBLIC :: reservoir_operation
   PUBLIC :: reservoir_final

CONTAINS

   ! -------
   SUBROUTINE reservoir_init ( )

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Utils
   USE MOD_Namelist,              only: DEF_ReservoirPara_file, DEF_Reservoir_Method, &
      DEF_UnitCatchment_regional, regional_unitcatchment_file
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_ucid, lake_type
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: parafile

   integer,  allocatable :: dam_seq(:), order(:), loc2all(:), ordinal_count(:)
   integer,  allocatable :: catalogue_to_active(:)
   real(r8), allocatable :: rcache (:)
   integer,  allocatable :: icache (:)

   integer :: i, j, iloc, irsv, nresv, iworker, nresv_catalogue, ndropped
   integer,  allocatable :: src_index(:), regional_index(:)
   real(r8) :: flood_capacity, upstream_area


      parafile = DEF_ReservoirPara_file

      ! 1 = Funato et al.; 2 = Hanazaki et al. as implemented in CaMa-Flood PR #124.
      IF (DEF_Reservoir_Method /= 1 .and. DEF_Reservoir_Method /= 2) &
         CALL CoLM_stop ('unsupported reservoir operation method')

      CALL ncio_read_bcast_serial (parafile, 'dam_GRAND_ID', dam_GRAND_ID)
      CALL ncio_read_bcast_serial (parafile, 'dam_seq', dam_seq)

      nresv_catalogue = size(dam_seq)
      IF (size(dam_GRAND_ID) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_GRAND_ID and dam_seq lengths differ')

      IF (DEF_UnitCatchment_regional) THEN
         IF (p_is_master) THEN
            IF (ncio_var_exist(parafile, 'seq_src_index', readflag = .false.)) &
               CALL CoLM_stop ('DEF_ReservoirPara_file must use the source unit catchment numbering')
         ENDIF
         CALL ncio_read_bcast_serial (regional_unitcatchment_file (), 'seq_src_index', src_index)
         allocate (regional_index (max(maxval(src_index), 1)))
         regional_index = 0
         DO i = 1, size(src_index)
            regional_index(src_index(i)) = i
         ENDDO
         DO i = 1, nresv_catalogue
            iloc = 0
            IF (dam_seq(i) >= 1 .and. dam_seq(i) <= size(regional_index)) iloc = regional_index(dam_seq(i))
            IF (iloc > 0) THEN
               dam_seq(i) = iloc
            ELSE
               dam_seq(i) = -i
            ENDIF
         ENDDO
         deallocate (src_index, regional_index)
      ENDIF

      numresv = 0

      allocate (order (nresv_catalogue))
      order = (/(i, i = 1, nresv_catalogue)/)

      CALL quicksort (nresv_catalogue, dam_seq, order)

      ndropped = 0

      ! One reservoir per unit catchment: when several rows share a dam_seq, keep the
      ! first row in the file and drop the others (their seq becomes -row, never matched).
      IF (nresv_catalogue > 1) THEN
         IF (any(dam_seq(2:nresv_catalogue) == dam_seq(1:nresv_catalogue-1))) THEN
            allocate (icache (nresv_catalogue))
            icache(order) = dam_seq
            i = 1
            DO WHILE (i < nresv_catalogue)
               iloc = i
               DO WHILE (iloc < nresv_catalogue)
                  IF (dam_seq(iloc+1) /= dam_seq(i)) EXIT
                  iloc = iloc + 1
               ENDDO
               IF (iloc > i) THEN
                  irsv = minval(order(i:iloc))
                  IF (p_is_master) write(*,'(A,I0,A,I0,A)') 'Warning: ', iloc-i+1, &
                     ' dams share unit catchment ', dam_seq(i), ', only the first in the file is kept:'
                  DO j = i, iloc
                     IF (p_is_master) write(*,'(A,I0,A)') '   GRAND_ID ', dam_GRAND_ID(order(j)), &
                        trim(merge(' (kept)   ', ' (dropped)', order(j) == irsv))
                     IF (order(j) /= irsv) THEN
                        icache(order(j)) = -order(j)
                        ndropped = ndropped + 1
                     ENDIF
                  ENDDO
               ENDIF
               i = iloc + 1
            ENDDO
            dam_seq = icache
            order   = (/(i, i = 1, nresv_catalogue)/)
            CALL quicksort (nresv_catalogue, dam_seq, order)
            deallocate (icache)
         ENDIF
      ENDIF

      IF (p_is_worker) THEN

         allocate (ucat2resv (numucat))
         allocate (loc2all   (numucat))
         ucat2resv = 0

         numresv = 0
         DO i = 1, numucat
            iloc = find_in_sorted_list1 (ucat_ucid(i), nresv_catalogue, dam_seq)
            IF (iloc > 0) THEN
               numresv = numresv + 1
               lake_type(i) = 2
               ucat2resv(i) = numresv
               loc2all  (numresv) = order(iloc)
            ENDIF
         ENDDO

         allocate (resv_ucid     (numresv))
         IF (numresv > 0) THEN
            DO i = 1, numucat
               IF (ucat2resv(i) > 0) resv_ucid(ucat2resv(i)) = ucat_ucid(i)
            ENDDO
         ENDIF

      ENDIF

      allocate (ordinal_count(nresv_catalogue))
      ordinal_count = 0
      IF (p_is_worker) THEN
         DO irsv = 1, numresv
            IF (loc2all(irsv) < 1 .or. loc2all(irsv) > nresv_catalogue) &
               CALL CoLM_stop ('invalid reservoir parameter-table ordinal')
            ordinal_count(loc2all(irsv)) = ordinal_count(loc2all(irsv)) + 1
         ENDDO
      ENDIF
#ifdef USEMPI
      IF (nresv_catalogue > 0) CALL mpi_allreduce (MPI_IN_PLACE, ordinal_count, nresv_catalogue, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (nresv_catalogue > 0) THEN
         IF (any(ordinal_count > 1)) THEN
            IF (p_is_master) write(*,'(A,I0)') &
               'ERROR: reservoir rows mapped to multiple active ucatch owners=', count(ordinal_count > 1)
            CALL CoLM_stop ('active reservoir ucatch IDs must map one-to-one to parameter rows')
         ENDIF
         IF (p_is_master .and. count(ordinal_count == 0) > ndropped) write(*,'(A,I0,A,I0)') &
            'Reservoir catalogue rows outside active domain: ', count(ordinal_count == 0) - ndropped, &
            ' of ', nresv_catalogue
      ENDIF

      allocate (catalogue_to_active(nresv_catalogue))
      catalogue_to_active = 0
      totalnumresv = 0
      DO i = 1, nresv_catalogue
         IF (ordinal_count(i) > 0) THEN
            totalnumresv = totalnumresv + 1
            catalogue_to_active(i) = totalnumresv
         ENDIF
      ENDDO
      dam_GRAND_ID = pack(dam_GRAND_ID, ordinal_count > 0)
      IF (p_is_worker) THEN
         allocate (resv_global_id(numresv))
         IF (numresv > 0) resv_global_id = catalogue_to_active(loc2all(1:numresv))
      ENDIF
      deallocate (ordinal_count)
      deallocate (catalogue_to_active)

#ifdef USEMPI
      IF (p_is_master) THEN

         allocate (resv_data_address (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (nresv, 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (nresv > 0) THEN
               allocate (resv_data_address(iworker)%val (nresv))
               CALL mpi_recv (resv_data_address(iworker)%val, nresv, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_send (numresv, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numresv > 0) THEN
            CALL mpi_send (resv_global_id, numresv, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
         ENDIF

      ENDIF

      IF (.not. p_is_master .and. .not. allocated(resv_data_address)) THEN
         allocate (resv_data_address (0:-1))
      ENDIF
#else
      IF (numresv > 0) THEN
         allocate (resv_data_address (0:0))
         allocate (resv_data_address(0)%val (numresv))
         resv_data_address(0)%val = resv_global_id
      ENDIF
#endif

      IF (p_is_worker) THEN

         IF (numresv > 0) THEN

            allocate (dam_build_year (numresv))

            allocate (volresv_total  (numresv))
            allocate (volresv_emerg  (numresv))
            allocate (volresv_adjust (numresv))
            allocate (volresv_normal (numresv))
            IF (DEF_Reservoir_Method == 2) THEN
               allocate (volresv_flood (numresv), resv_h22_coeff (numresv))
            ENDIF

            allocate (qresv_flood    (numresv))
            allocate (qresv_adjust   (numresv))
            allocate (qresv_normal   (numresv))

            allocate (qresv_in       (numresv))
            allocate (qresv_out      (numresv))

         ENDIF

      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_year', icache)
      IF (size(icache) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_year and dam_seq lengths differ')
      IF (p_is_worker .and. (numresv > 0)) THEN
         dam_build_year = icache(loc2all(1:numresv))
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_TotalVol_mcm', rcache)
      IF (size(rcache) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_TotalVol_mcm and dam_seq lengths differ')
      IF (p_is_worker .and. (numresv > 0)) THEN
         volresv_total = rcache(loc2all(1:numresv))*1.e6
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_ConVol_mcm', rcache)
      IF (size(rcache) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_ConVol_mcm and dam_seq lengths differ')
      IF (p_is_worker .and. (numresv > 0)) THEN
         volresv_normal = rcache(loc2all(1:numresv))*1.e6
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_Qn', rcache)
      IF (size(rcache) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_Qn and dam_seq lengths differ')
      IF (p_is_worker .and. (numresv > 0)) THEN
         qresv_normal = rcache(loc2all(1:numresv))
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_Qf', rcache)
      IF (size(rcache) /= nresv_catalogue) &
         CALL CoLM_stop ('reservoir dam_Qf and dam_seq lengths differ')
      IF (p_is_worker .and. (numresv > 0)) THEN
         qresv_flood = rcache(loc2all(1:numresv))
      ENDIF


      ! Read collectively using the original catalogue order, including on empty workers.
      ! A globally inactive catalogue does not require H22-specific parameters.
      IF (DEF_Reservoir_Method == 2 .and. totalnumresv > 0) THEN
         CALL ncio_read_bcast_serial (parafile, 'dam_upreal', rcache)
         IF (size(rcache) /= nresv_catalogue) &
            CALL CoLM_stop ('reservoir dam_upreal and dam_seq lengths differ')
      ENDIF

      IF (p_is_worker) THEN
         DO irsv = 1, numresv
            IF (.not. ieee_is_finite(volresv_total(irsv)) .or. volresv_total(irsv) <= 0._r8 .or. &
                .not. ieee_is_finite(volresv_normal(irsv)) .or. volresv_normal(irsv) <= 0._r8 .or. &
                .not. ieee_is_finite(qresv_normal(irsv)) .or. qresv_normal(irsv) < 0._r8 .or. &
                .not. ieee_is_finite(qresv_flood(irsv)) .or. qresv_flood(irsv) < 0._r8) &
               CALL CoLM_stop ('invalid active reservoir volume or outflow parameter')
            IF (DEF_Reservoir_Method == 1) THEN
               volresv_emerg (irsv) = volresv_total(irsv) * 0.94
               volresv_adjust(irsv) = volresv_total(irsv) * 0.77
               volresv_normal(irsv) = min(volresv_total(irsv)*0.7, volresv_normal(irsv))
               qresv_adjust  (irsv) = (qresv_normal(irsv) + qresv_flood(irsv)) * 0.5
            ELSE
               upstream_area = rcache(loc2all(irsv))  ! km^2
               IF (.not. ieee_is_finite(upstream_area) .or. upstream_area <= 0._r8 .or. &
                   volresv_normal(irsv) > volresv_total(irsv) .or. &
                   qresv_flood(irsv) < qresv_normal(irsv)) &
                  CALL CoLM_stop ('invalid active H22 reservoir storage, discharge or dam_upreal')
               volresv_flood(irsv) = volresv_normal(irsv)
               flood_capacity = volresv_total(irsv) - volresv_flood(irsv)
               volresv_normal(irsv) = 0.5_r8 * volresv_flood(irsv)
               ! Match the legacy LDAMH22 branch of CaMa-Flood PR #124: 0.95, not 0.8.
               volresv_emerg(irsv) = volresv_flood(irsv) + 0.95_r8 * flood_capacity
               resv_h22_coeff(irsv) = max(1._r8 - flood_capacity * 1.e-6_r8 / (upstream_area * 0.2_r8), 0._r8)
               ! Adjustment parameters are unused by H22.
               volresv_adjust(irsv) = 0._r8
               qresv_adjust(irsv) = 0._r8
            ENDIF
         ENDDO
      ENDIF

      IF (allocated(dam_seq)) deallocate(dam_seq)
      IF (allocated(order  )) deallocate(order  )
      IF (allocated(loc2all)) deallocate(loc2all)
      IF (allocated(rcache )) deallocate(rcache )
      IF (allocated(icache )) deallocate(icache )

   END SUBROUTINE reservoir_init


   SUBROUTINE reservoir_operation (method, irsv, qin, vol, qout)

   USE MOD_SPMD_Task, only: CoLM_stop
   IMPLICIT NONE
   integer,  intent(in)  :: method
   integer,  intent(in)  :: irsv
   real(r8), intent(in)  :: qin, vol
   real(r8), intent(out) :: qout

   ! local variables
   real(r8) :: q1

      SELECT CASE (method)
      CASE (1)
         ! *** Reference ***
         ! [1] Mizuki Funato, Dai Yamazaki, Dung Trung Vu.
         ! Development of an Improved Reservoir Operation Scheme for Global Flood Modeling (CaMa-Flood v4.20).
         ! ESS Open Archive . October 24, 2024.

         IF (vol > volresv_emerg(irsv)) THEN
            qout = max(qin, qresv_flood(irsv))
         ELSEIF (vol > volresv_adjust(irsv)) THEN
            qout = qresv_adjust(irsv) + (qresv_flood(irsv)-qresv_adjust(irsv)) &
               * ((vol-volresv_adjust(irsv))/(volresv_emerg(irsv)-volresv_adjust(irsv)))**0.1
            IF (qin > qresv_flood(irsv)) THEN
               q1 = qresv_normal(irsv) + (qin-qresv_normal(irsv)) &
                  * (vol-volresv_normal(irsv))/(volresv_emerg(irsv)-volresv_normal(irsv))
               qout = max(q1, qout)
            ENDIF
         ELSEIF (vol > volresv_normal(irsv)) THEN
            qout = qresv_normal(irsv) + (qresv_adjust(irsv)-qresv_normal(irsv)) &
               * ((vol-volresv_normal(irsv))/(volresv_adjust(irsv)-volresv_normal(irsv)))**3.
         ELSE
            qout = (vol/volresv_normal(irsv))**0.5 * qresv_normal(irsv)
         ENDIF

      CASE (2)
         ! Hanazaki et al. (2022), legacy LDAMH22 in CaMa-Flood PR #124.
         ! Preserve its exact interval boundaries, including the high-inflow jump at C.
         IF (vol <= volresv_normal(irsv)) THEN
            qout = qresv_normal(irsv) * (vol / volresv_flood(irsv))
         ELSEIF (vol <= volresv_flood(irsv)) THEN
            IF (qin >= qresv_flood(irsv)) THEN
               qout = 0.5_r8 * qresv_normal(irsv) &
                  + (vol-volresv_normal(irsv)) / (volresv_flood(irsv)-volresv_normal(irsv)) &
                  * (qresv_flood(irsv)-qresv_normal(irsv))
            ELSE
               qout = 0.5_r8 * qresv_normal(irsv) &
                  + ((vol-volresv_normal(irsv)) / (volresv_emerg(irsv)-volresv_normal(irsv)))**2 &
                  * (qresv_flood(irsv)-qresv_normal(irsv))
            ENDIF
         ELSEIF (vol < volresv_emerg(irsv)) THEN
            IF (qin >= qresv_flood(irsv)) THEN
               qout = qresv_flood(irsv) + resv_h22_coeff(irsv) &
                  * (vol-volresv_flood(irsv)) / (volresv_emerg(irsv)-volresv_flood(irsv)) &
                  * (qin-qresv_flood(irsv))
            ELSE
               qout = 0.5_r8 * qresv_normal(irsv) &
                  + ((vol-volresv_normal(irsv)) / (volresv_emerg(irsv)-volresv_normal(irsv)))**2 &
                  * (qresv_flood(irsv)-qresv_normal(irsv))
            ENDIF
         ELSE
            qout = max(qin, qresv_flood(irsv))
         ENDIF

      CASE DEFAULT
         CALL CoLM_stop ('unsupported reservoir operation method')
      END SELECT

   END SUBROUTINE reservoir_operation


   SUBROUTINE reservoir_final ()

   IMPLICIT NONE

      IF (allocated(ucat2resv        )) deallocate (ucat2resv        )
      IF (allocated(resv_global_id   )) deallocate (resv_global_id   )
      IF (allocated(resv_ucid        )) deallocate (resv_ucid        )
      IF (allocated(resv_data_address)) deallocate (resv_data_address)

      IF (allocated(dam_GRAND_ID     )) deallocate (dam_GRAND_ID     )
      IF (allocated(dam_build_year   )) deallocate (dam_build_year   )

      IF (allocated(volresv_total    )) deallocate (volresv_total    )
      IF (allocated(volresv_emerg    )) deallocate (volresv_emerg    )
      IF (allocated(volresv_adjust   )) deallocate (volresv_adjust   )
      IF (allocated(volresv_normal   )) deallocate (volresv_normal   )
      IF (allocated(volresv_flood    )) deallocate (volresv_flood    )
      IF (allocated(resv_h22_coeff   )) deallocate (resv_h22_coeff   )

      IF (allocated(qresv_flood      )) deallocate (qresv_flood      )
      IF (allocated(qresv_adjust     )) deallocate (qresv_adjust     )
      IF (allocated(qresv_normal     )) deallocate (qresv_normal     )

      IF (allocated(qresv_in         )) deallocate (qresv_in         )
      IF (allocated(qresv_out        )) deallocate (qresv_out        )

   END SUBROUTINE reservoir_final

END MODULE MOD_Grid_Reservoir
#endif
