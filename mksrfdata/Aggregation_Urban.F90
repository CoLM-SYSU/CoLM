#include <define.h>

!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  Aggregate/screen high-resolution urban dataset to a lower
!  resolution/subset data, suitable for running regional or point
!  cases.
!
!  Original authors: Hua Yuan and Wenzong Dong, 2021, OpenMP version.
!
! !REVISIONS:
!  05/2023, Wenzong Dong, Hua Yuan, Shupeng Zhang: porting codes to MPI
!           parallel version.
!
!  08/2025, Wenzong Dong, Hua Yuan: unifying the urban surface data
!           code for different urban type schemes.
!-----------------------------------------------------------------------

#ifdef URBAN_MODEL
SUBROUTINE Aggregation_Urban (grid_roof, grid_fgper, grid_pctt, grid_lsai, grid_pctw, grid_pop, grid_lucy, grid_alb, &
                              dir_rawdata, dir_model_landdata, lc_year)
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_AggregationRequestData
   USE MOD_5x5DataReadin
   USE MOD_DataType
   USE MOD_Utils, only: num_max_frequency
   USE MOD_LandUrban
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_Vars_Global, only: N_URB, nl_roof, nl_wall, nl_soil
   USE MOD_Urban_Const_LCZ
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: dir_rawdata
   character(len=256), intent(in) :: dir_model_landdata

   integer , intent(in) :: lc_year

   type(grid_type), intent(in) :: grid_lucy
   type(grid_type), intent(in) :: grid_lsai
   type(grid_type), intent(in) :: grid_roof
   type(grid_type), intent(in) :: grid_fgper
   type(grid_type), intent(in) :: grid_pctt, grid_pctw
   type(grid_type), intent(in) :: grid_pop, grid_alb

   ! dimensions
   integer, parameter :: rid  = 33
   integer, parameter :: ns   = 2
   integer, parameter :: nr   = 2
   integer, parameter :: ulev = 10

   ! input variables
   type(block_data_int32_2d) :: LUCY_reg
   type(block_data_real8_2d) :: pop
   type(block_data_real8_2d) :: fvegu
   type(block_data_real8_2d) :: htopu
   type(block_data_real8_2d) :: flakeu
   type(block_data_real8_2d) :: wtroof
   type(block_data_real8_2d) :: hlr
   type(block_data_real8_2d) :: htroof
   type(block_data_real8_2d) :: fgper
   type(block_data_real8_2d) :: ulai
   type(block_data_real8_2d) :: usai
   type(block_data_real8_2d) :: albroof
   type(block_data_int32_2d) :: reg_typid

   ! output variables
   integer , allocatable, dimension(:) :: LUCY_rid
   real(r8), allocatable, dimension(:) :: pop_den
   real(r8), allocatable, dimension(:) :: pct_tree
   real(r8), allocatable, dimension(:) :: htop_urb
   real(r8), allocatable, dimension(:) :: pct_water
   real(r8), allocatable, dimension(:) :: pct_roof
   real(r8), allocatable, dimension(:) :: ht_roof
   real(r8), allocatable, dimension(:) :: lai_urb
   real(r8), allocatable, dimension(:) :: sai_urb

   integer , allocatable, dimension(:) :: reg_typid_one
   integer , allocatable, dimension(:) :: LUCY_reg_one
   real(r8), allocatable, dimension(:) :: area_one
   real(r8), allocatable, dimension(:) :: pop_one
   real(r8), allocatable, dimension(:) :: fvegu_one
   real(r8), allocatable, dimension(:) :: htopu_one
   real(r8), allocatable, dimension(:) :: flakeu_one
   real(r8), allocatable, dimension(:) :: wt_roof_one
   real(r8), allocatable, dimension(:) :: ht_roof_one
   real(r8), allocatable, dimension(:) :: hlr_bld_one
   real(r8), allocatable, dimension(:) :: f_gper_one
   real(r8), allocatable, dimension(:) :: ulai_one
   real(r8), allocatable, dimension(:) :: usai_one
   real(r8), allocatable, dimension(:) :: alb_roof_one_

   real(r8), allocatable, dimension(:,:,:) :: alb_roof_one, alb_wall_one, &
                                              alb_gimp_one, alb_gper_one

   ! urban morphological and thermal paras of NCAR data
   ! input variables, NCAR look-up-table data
   real(r8), allocatable, dimension(:,:)     :: hwrbld_ncar , fgper_ncar  , &
                                                htroof_ncar , wtroof_ncar
   real(r8), allocatable, dimension(:,:)     :: emroof_ncar , emwall_ncar , &
                                                emgimp_ncar , emgper_ncar
   real(r8), allocatable, dimension(:,:)     :: thkroof_ncar, thkwall_ncar, &
                                                tbldmin_ncar, tbldmax_ncar
   real(r8), allocatable, dimension(:,:,:)   :: cvroof_ncar , cvwall_ncar , cvgimp_ncar , &
                                                tkroof_ncar , tkwall_ncar , tkgimp_ncar
   real(r8), allocatable, dimension(:,:,:,:) :: albroof_ncar, albwall_ncar, &
                                                albgimp_ncar, albgper_ncar

   ! output variables, vector data
   real(r8), allocatable, dimension(:)     :: area_urb
   real(r8), allocatable, dimension(:)     :: sarea_urb
   real(r8), allocatable, dimension(:)     :: urb_pct

   real(r8), allocatable, dimension(:)     :: hlr_bld
   real(r8), allocatable, dimension(:)     :: f_gper
   real(r8), allocatable, dimension(:)     :: f_gimp
   real(r8), allocatable, dimension(:)     :: f_isa
   real(r8), allocatable, dimension(:)     :: em_roof
   real(r8), allocatable, dimension(:)     :: em_wall
   real(r8), allocatable, dimension(:)     :: em_gimp
   real(r8), allocatable, dimension(:)     :: em_gper
   real(r8), allocatable, dimension(:)     :: thk_roof
   real(r8), allocatable, dimension(:)     :: thk_wall
   real(r8), allocatable, dimension(:)     :: tbld_min
   real(r8), allocatable, dimension(:)     :: tbld_max

   real(r8), allocatable, dimension(:,:)   :: cv_wgt
   real(r8), allocatable, dimension(:,:)   :: tk_wgt
   real(r8), allocatable, dimension(:,:)   :: cv_roof
   real(r8), allocatable, dimension(:,:)   :: cv_wall
   real(r8), allocatable, dimension(:,:)   :: cv_gimp
   real(r8), allocatable, dimension(:,:)   :: tk_roof
   real(r8), allocatable, dimension(:,:)   :: tk_wall
   real(r8), allocatable, dimension(:,:)   :: tk_gimp

   real(r8), allocatable, dimension(:,:,:) :: alb_roof
   real(r8), allocatable, dimension(:,:,:) :: alb_wall
   real(r8), allocatable, dimension(:,:,:) :: alb_gimp
   real(r8), allocatable, dimension(:,:,:) :: alb_gper

   integer , allocatable, dimension(:)     :: locpth

   ! landfile variables
   character(len=256) dir, fname
   character(len=256) landdir, landname
   character(len=4  ) cyear, c5year, cmonth, clay

   ! local vars
   real(r8) :: sumarea
   real(r8) :: thkroof, thkwall, tbldmin, tbldmax
   real(r8) :: emroof, emwall, emgper , emgimp
   real(r8) :: cvroof(nl_roof), cvwall(nl_wall), cvgimp(nl_soil)
   real(r8) :: tkroof(nl_roof), tkwall(nl_wall), tkgimp(nl_soil)

   ! index
   integer :: iurban, urb_typidx, urb_regidx
   integer :: pop_i, imonth, start_year, end_year
   integer :: ipth, ipxl, il, iy, ielm, numpth, urb_s, urb_e

   ! for surface data diag
#ifdef SrfdataDiag
   integer  :: ityp
   integer , allocatable, dimension(:) :: typindex
   real(r8), allocatable :: LUCY_rid_r8 (:)
#endif
   logical  :: first_call_LAI_urban, first_call_SAI_urban


#ifdef SrfdataDiag
      allocate( typindex(N_URB) )
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/urban/' // trim(cyear)

      first_call_LAI_urban = .true.
      first_call_SAI_urban = .true.

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Making urban data ('//trim(cyear)//') ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      write(c5year, '(i4.4)') int(lc_year/5)*5

      ! ******* Building morphology: PCT_ROOF, HT_ROOF, HL and Fgper *******
      ! if building data is missing, how to use look-up-table?
      ! a new array with urban type was used for look-up-table
IF (DEF_URBAN_type_scheme == 1) THEN
      ! only used when urban patch have nan data of building height, building fraction and HL
      landname = trim(dir_rawdata)//'/urban_physical/NCAR_urban_properties.nc'

      CALL ncio_read_bcast_serial (landname,  "WTLUNIT_ROOF", wtroof_ncar )
      CALL ncio_read_bcast_serial (landname,  "HT_ROOF"     , htroof_ncar )
      CALL ncio_read_bcast_serial (landname,  "CANYON_HWR"  , hwrbld_ncar )
      CALL ncio_read_bcast_serial (landname,  "WTROAD_PERV" , fgper_ncar  )
ENDIF

      ! allocate and read grided building height and fraction raw data
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_urban, reg_typid)
         CALL allocate_block_data (grid_roof , wtroof   )
         CALL allocate_block_data (grid_roof , htroof   )
         CALL allocate_block_data (grid_roof , hlr      )
         CALL allocate_block_data (grid_fgper, fgper    )

         CALL flush_block_data (wtroof, -999.)
         CALL flush_block_data (htroof, -999.)
         CALL flush_block_data (hlr   , -999.)
         CALL flush_block_data (fgper , -999.)

IF (DEF_URBAN_type_scheme == 1) THEN
         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_type%dir)
         fname = trim(DEF_rawdata%urban_type%fname)
         CALL read_5x5_data (dir, fname, grid_urban, "REGION_ID", reg_typid)
ENDIF

         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_roof%dir)

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            fname = trim(DEF_rawdata%urban_roof%fname)//trim(c5year)

            IF (index(DEF_rawdata%urban_roof%vname, 'GHSL')>0) THEN
               CALL read_5x5_data (dir, fname, grid_roof, "PCT_ROOF_GHSL", wtroof)
               CALL read_5x5_data (dir, fname, grid_roof, "HT_ROOF_GHSL" , htroof)
               CALL read_5x5_data (dir, fname, grid_roof, "HL_GHSL"      , hlr   )
            ELSE
               CALL read_5x5_data (dir, fname, grid_roof, "PCT_ROOF_Li", wtroof)
               CALL read_5x5_data (dir, fname, grid_roof, "HT_ROOF_Li" , htroof)
               CALL read_5x5_data (dir, fname, grid_roof, "HL_Li"      , hlr   )
            ENDIF

         ELSE

            IF (index(DEF_rawdata%urban_roof%fname, 'GHSL')>0) THEN
               fname = trim(DEF_rawdata%urban_roof%fname)//'.'//trim(c5year)
            ELSE
               fname = trim(DEF_rawdata%urban_roof%fname)
            ENDIF

            CALL read_5x5_data (dir, fname, grid_roof, "PCT_ROOF", wtroof)
            CALL read_5x5_data (dir, fname, grid_roof, "HT_ROOF" , htroof)

            dir= trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_hl%dir)
            fname = trim(DEF_rawdata%urban_hl%fname)
            CALL read_5x5_data (dir, fname, grid_roof, "HL_BLD"  , hlr   )

            dir= trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_fgper%dir)
            fname = trim(DEF_rawdata%urban_fgper%fname)
            CALL read_5x5_data (dir, fname, grid_fgper, "WTROAD_PERV", fgper)
         ENDIF

#ifdef USEMPI
IF (DEF_URBAN_type_scheme == 1) THEN
         CALL aggregation_data_daemon_multigrid (grid_in1 = grid_urban, data_i4_2d_in1 = reg_typid, &
            grid_in2 = grid_roof , data_r8_2d_in2 = wtroof, &
            grid_in3 = grid_roof , data_r8_2d_in3 = htroof, &
            grid_in4 = grid_roof , data_r8_2d_in4 = hlr   , &
            grid_in5 = grid_fgper, data_r8_2d_in5 = fgper)
ELSE
         CALL aggregation_data_daemon_multigrid (grid_in1 = grid_roof , data_r8_2d_in1 = wtroof, &
            grid_in2 = grid_roof , data_r8_2d_in2 = htroof, &
            grid_in3 = grid_roof , data_r8_2d_in3 = hlr   , &
            grid_in4 = grid_fgper, data_r8_2d_in4 = fgper)
ENDIF
#endif
      ENDIF

      IF (p_is_worker) THEN
         allocate (pct_roof (numurban))
         allocate (ht_roof  (numurban))
         allocate (hlr_bld  (numurban))
         allocate (f_gper   (numurban))
         allocate (f_gimp   (numurban))
         allocate (f_isa    (numurban))

         f_gimp (:) = 0.
         f_isa  (:) = 0.

         ! loop for urban patch to aggregate building height and fraction data with area-weighted average
         DO iurban = 1, numurban

            ! when urban patch has no data, use table data to fill gap
            ! urban type and region id for look-up-table
            urb_typidx = landurban%settyp(iurban)

IF (DEF_URBAN_type_scheme == 1) THEN
            CALL aggregation_request_data_multigrid (landurban, iurban, &
               grid_in1=grid_urban, area = area_one, data_i4_2d_in1 = reg_typid, data_i4_2d_out1 = reg_typid_one, &
               grid_in2=grid_roof , data_r8_2d_in2 = wtroof   , data_r8_2d_out2 = wt_roof_one  , &
               grid_in3=grid_roof , data_r8_2d_in3 = htroof   , data_r8_2d_out3 = ht_roof_one  , &
               grid_in4=grid_roof , data_r8_2d_in4 = hlr      , data_r8_2d_out4 = hlr_bld_one  , &
               grid_in5=grid_fgper, data_r8_2d_in5 = fgper    , data_r8_2d_out5 = f_gper_one)

            ! RG_-45_65_-50_70 of NCAR has no urban data,
            ! all urban patches of this area are assigned to region 30
            IF (all(reg_typid_one==0)) THEN
               reg_typid_one(:) = 30
            ENDIF

            IF (any(reg_typid_one==0)) THEN
               WHERE(reg_typid_one==0) reg_typid_one = num_max_frequency(reg_typid_one)
            ENDIF

            IF (all(wt_roof_one <= 0) .or. all(ht_roof_one <= 0)) THEN
               wt_roof_one = wtroof_ncar(urb_typidx,reg_typid_one)
               ht_roof_one = htroof_ncar(urb_typidx,reg_typid_one)
            ENDIF

            IF (all(f_gper_one <= 0)) THEN
               f_gper_one = fgper_ncar(urb_typidx,reg_typid_one)
            ENDIF

IF(DEF_USE_CANYON_HWR) THEN
            hlr_bld_one = -999.
ENDIF

            WHERE (hlr_bld_one<=0 .or. hlr_bld_one>=5)
               hlr_bld_one = hwrbld_ncar(urb_typidx,reg_typid_one) &
                           *(1-sqrt(wtroof_ncar(urb_typidx,reg_typid_one))) &
                           /sqrt(wtroof_ncar(urb_typidx,reg_typid_one))
            END WHERE

            WHERE (wt_roof_one<0 .or. ht_roof_one<0 .or. f_gper_one<0)
               area_one = 0
            END WHERE

ELSE IF (DEF_URBAN_type_scheme == 2) THEN
            CALL aggregation_request_data_multigrid (landurban, iurban, &
               grid_in1 = grid_roof , area = area_one, data_r8_2d_in1 = wtroof, data_r8_2d_out1 = wt_roof_one, &
               grid_in2 = grid_roof , data_r8_2d_in2 = htroof, data_r8_2d_out2 = ht_roof_one, &
               grid_in3 = grid_roof , data_r8_2d_in3 = hlr   , data_r8_2d_out3 = hlr_bld_one, &
               grid_in4 = grid_fgper, data_r8_2d_in4 = fgper , data_r8_2d_out4 = f_gper_one)

            IF (all(wt_roof_one <= 0) .or. all(ht_roof_one <= 0)) THEN
               wt_roof_one = wtroof_lcz(urb_typidx)
               ht_roof_one = htroof_lcz(urb_typidx)
            ENDIF

            WHERE (f_gper_one <= 0)
               f_gper_one = fgper_lcz(urb_typidx) / (1-wtroof_lcz(urb_typidx))
            END WHERE

IF(DEF_USE_CANYON_HWR) THEN
            hlr_bld_one = -999.
ENDIF

            WHERE (hlr_bld_one<=0 .or. hlr_bld_one>=5)
               hlr_bld_one = hwrbld_lcz(urb_typidx) &
                           *(1-sqrt(wtroof_lcz(urb_typidx))) &
                           /sqrt(wtroof_lcz(urb_typidx))
            END WHERE

            WHERE (wt_roof_one<0 .or. ht_roof_one<0)
               area_one = 0
            END WHERE

ENDIF
            ! area-weight average
            pct_roof(iurban) = sum(wt_roof_one * area_one) / sum(area_one)
            ht_roof (iurban) = sum(ht_roof_one * area_one) / sum(area_one)
            hlr_bld (iurban) = sum(hlr_bld_one * area_one) / sum(area_one)
            f_gper  (iurban) = sum(f_gper_one  * area_one) / sum(area_one)

            f_gimp  (iurban) = 1. - f_gper(iurban)
            f_isa   (iurban) = f_gimp(iurban)*(1.-pct_roof(iurban)) + pct_roof(iurban)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done_multigrid ()
#endif
      ENDIF

      ! output
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/PCT_ROOF.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'PCT_ROOF', 'urban', landurban, pct_roof, DEF_Srfdata_CompressLevel)

      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/HT_ROOF.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'HT_ROOF', 'urban', landurban, ht_roof, DEF_Srfdata_CompressLevel)

      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/HLR_BLD.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'BUILDING_HLR'  , 'urban', landurban, hlr_bld, DEF_Srfdata_CompressLevel)

      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/WTROAD_PERV.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'WTROAD_PERV', 'urban', landurban, f_gper, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_morphology_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (ht_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'HT_ROOF', compress = 6, write_mode = 'one', create_mode = .true.)

      CALL srfdata_map_and_write (pct_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'PCT_ROOF', compress = 6, write_mode = 'one', defval = 0._r8, create_mode = .false.)

      CALL srfdata_map_and_write (hlr_bld, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'HL', compress = 6, write_mode = 'one', create_mode = .false.)

      CALL srfdata_map_and_write (f_gimp, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'WT_IMPG', compress = 6, write_mode = 'one', &
         defval = 0._r8, create_mode = .false.)

      CALL srfdata_map_and_write (f_isa, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'PCT_ISA', compress = 6, write_mode = 'one', &
         defval = 0._r8, create_mode = .false.)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('Urban Roof Fraction ', pct_roof)
      CALL check_vector_data ('Urban Roof Height '  , ht_roof )
      CALL check_vector_data ('Urban Building HLR ' , hlr_bld )
      CALL check_vector_data ('WTROAD_PERV '        , f_gper  )
#endif


      ! ******* Tree : PCT_Tree, HTOP *******
      ! allocate and read the grided tree cover and tree height raw data(500m)
      ! NOTE, tree cover raw data is available every five years,
      ! tree height raw data is same from year to year
      IF (p_is_io) THEN

         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_fveg%dir)

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            fname = trim(DEF_rawdata%urban_fveg%fname)//trim(c5year)
         ELSE
            fname = trim(DEF_rawdata%urban_fveg%fname)//'.'//trim(c5year)
         ENDIF

         CALL allocate_block_data (grid_pctt, fvegu)
         CALL read_5x5_data (dir, fname, grid_pctt, "PCT_Tree", fvegu)

         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_htop%dir)

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            fname = trim(DEF_rawdata%urban_fveg%fname)//trim(c5year)
         ELSE
            fname = trim(DEF_rawdata%urban_htop%fname)
         ENDIF

         CALL allocate_block_data (grid_pctt, htopu)
         CALL read_5x5_data (dir, fname, grid_pctt, trim(DEF_rawdata%urban_htop%vname), htopu)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_pctt, &
            data_r8_2d_in1 = fvegu, data_r8_2d_in2 = htopu)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (pct_tree(numurban))
         allocate (htop_urb(numurban))

         pct_tree(:) = 0.
         htop_urb(:) = 0.

         ! loop for urban patch to aggregate tree cover and height data with area-weighted average
         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid_pctt, zip = USE_zip_for_aggregation, area = area_one, &
               data_r8_2d_in1 = fvegu, data_r8_2d_out1 = fvegu_one, &
               data_r8_2d_in2 = htopu, data_r8_2d_out2 = htopu_one)

            ! missing tree cover and tree height data (-999) were filtered
            WHERE (fvegu_one < 0)
               area_one = 0
            END WHERE

            WHERE (htopu_one<0 .or. fvegu_one<0)
               area_one = 0
            END WHERE

            ! area-weighted average
            IF (sum(area_one) > 0._r8) THEN
               pct_tree(iurban) = sum(fvegu_one * area_one) / sum(area_one)
               htop_urb(iurban) = sum(htopu_one * area_one * fvegu_one) / sum(area_one * fvegu_one)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/PCT_Tree.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'PCT_Tree', 'urban', landurban, pct_tree, DEF_Srfdata_CompressLevel)

      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/htop_urb.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'URBAN_TREE_TOP', 'urban', landurban, htop_urb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_ecology_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (pct_tree, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'PCT_TREE', compress = 6, write_mode = 'one', defval = 0._r8, create_mode = .true.)

      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_ecology_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (htop_urb, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'HTOP', compress = 6, write_mode = 'one', defval = 0._r8, create_mode = .false.)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('Urban Tree Cover ', pct_tree)
      CALL check_vector_data ('Urban Tree Top '  , htop_urb)
#endif


      ! ******* LAI, SAI *******
#ifndef LULCC
IF (DEF_LAI_CHANGE_YEARLY) THEN
      start_year = max(DEF_LAI_START_YEAR, DEF_simulation_time%start_year)
      end_year   = min(DEF_LAI_END_YEAR,   DEF_simulation_time%end_year  )
ELSE
      start_year = lc_year
      end_year   = lc_year
ENDIF
#else
      start_year = lc_year
      end_year   = lc_year
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_lsai, ulai)
      ENDIF

      IF (p_is_worker) THEN
         allocate (lai_urb (numurban))

         lai_urb(:) = 0.
      ENDIF

      DO iy = start_year, end_year

         IF (iy < 2000) THEN
            write(cyear,'(i4.4)') 2000
         ELSE
            write(cyear,'(i4.4)') iy
         ENDIF

         landdir = trim(dir_model_landdata) // '/urban/' // trim(cyear) // '/LAI'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))

         ! allocate and read grided LAI raw data
         dir = trim(dir_rawdata) // trim(DEF_rawdata%urban_lsai%dir)
         fname = trim(DEF_rawdata%urban_lsai%fname)//'.'//trim(cyear)

         ! loop for month
         DO imonth = 1, 12

            write(cmonth, '(i2.2)') imonth

            IF (p_is_master) THEN
               write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate LAI :', iy, ':', imonth, '/', 12
            ENDIF

            IF (p_is_io) THEN

               CALL read_5x5_data_time (dir, fname, grid_lsai, "URBAN_TREE_LAI", imonth, ulai)

#ifdef USEMPI
               CALL aggregation_data_daemon_multigrid (grid_in1 = grid_pctt, data_r8_2d_in1 = fvegu, &
                  grid_in2 = grid_lsai, data_r8_2d_in2 = ulai)
#endif
            ENDIF

            IF (p_is_worker) THEN

               ! loop for urban patch to aggregate LAI data
               DO iurban = 1, numurban
                  CALL aggregation_request_data_multigrid (landurban, iurban, &
                     grid_in1 = grid_pctt, area = area_one, data_r8_2d_in1 = fvegu, data_r8_2d_out1 = fvegu_one, &
                     grid_in2 = grid_lsai, data_r8_2d_in2 = ulai, data_r8_2d_out2 = ulai_one)

                  WHERE (fvegu_one<0 .or. ulai_one<0)
                     area_one = 0
                  END WHERE

                  ! area-weight average
                  IF (sum(area_one) > 0) THEN
                     lai_urb(iurban) = sum(ulai_one * area_one * fvegu_one) / &
                                       sum(area_one * fvegu_one)
                  ENDIF
               ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done_multigrid ()
#endif
            ENDIF

            ! output
            landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/LAI/urban_LAI_'//trim(cmonth)//'.nc'
            CALL ncio_create_file_vector (landname, landurban)
            CALL ncio_define_dimension_vector (landname, landurban, 'urban')
            CALL ncio_write_vector (landname, 'TREE_LAI', 'urban', landurban, lai_urb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
            typindex = (/(ityp, ityp = 1, N_URB)/)
            landname  = trim(dir_model_landdata) // '/diag/urban_ecology_'//trim(cyear)//'.nc'
            CALL srfdata_map_and_write (lai_urb, landurban%settyp, typindex, m_urb2diag, &
                  -1.0e36_r8, landname, 'Urban_Tree_LAI', compress = 6, write_mode = 'one',  &
                  lastdimname = 'Itime', lastdimvalue = imonth, defval = 0._r8, create_mode = .false.)!first_call_LAI_urban)

            IF (first_call_LAI_urban) first_call_LAI_urban = .false.
#endif

#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
            CALL check_vector_data ('Urban Tree LAI '//trim(cmonth), lai_urb)
#endif
         ENDDO
         first_call_LAI_urban = .true.
      ENDDO

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_lsai, usai)
      ENDIF

      IF (p_is_worker) THEN
         allocate (sai_urb (numurban))

         sai_urb(:) = 0.
      ENDIF

      DO iy = start_year, end_year

         IF (iy < 2000) THEN
            write(cyear,'(i4.4)') 2000
         ELSE
            write(cyear,'(i4.4)') iy
         ENDIF

         ! allocate and read grided SAI raw data
         dir = trim(dir_rawdata) // trim(DEF_rawdata%urban_lsai%dir)
         fname = trim(DEF_rawdata%urban_lsai%fname)//'.'//trim(cyear)

         ! loop for month
         DO imonth = 1, 12

            write(cmonth, '(i2.2)') imonth

            IF (p_is_master) THEN
               write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate SAI :', iy, ':', imonth, '/', 12
            ENDIF

            IF (p_is_io) THEN

               CALL read_5x5_data_time (dir, fname, grid_lsai, "URBAN_TREE_SAI", imonth, usai)

#ifdef USEMPI
               CALL aggregation_data_daemon_multigrid (grid_in1 = grid_pctt, data_r8_2d_in1 = fvegu, &
                  grid_in2 = grid_lsai, data_r8_2d_in2 = usai)
#endif
            ENDIF

            IF (p_is_worker) THEN
               ! loop for urban patch to aggregate SAI data
               DO iurban = 1, numurban
                  CALL aggregation_request_data_multigrid (landurban, iurban, &
                     grid_in1 = grid_pctt, area = area_one, data_r8_2d_in1 = fvegu, data_r8_2d_out1 = fvegu_one, &
                     grid_in2 = grid_lsai, data_r8_2d_in2 = usai , data_r8_2d_out2 = usai_one)

                  WHERE (fvegu_one<0 .or. usai_one<0)
                     area_one = 0
                  END WHERE

                  ! area-weight average
                  IF (sum(area_one) > 0) THEN
                     sai_urb(iurban) = sum(usai_one * area_one * fvegu_one) / &
                                       sum(area_one * fvegu_one)
                  ENDIF
               ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done_multigrid ()
#endif
            ENDIF

            ! output
            landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/LAI/urban_SAI_'//trim(cmonth)//'.nc'
            CALL ncio_create_file_vector (landname, landurban)
            CALL ncio_define_dimension_vector (landname, landurban, 'urban')
            CALL ncio_write_vector (landname, 'TREE_SAI', 'urban', landurban, sai_urb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
            landname  = trim(dir_model_landdata) // '/diag/urban_ecology_'//trim(cyear)//'.nc'
            CALL srfdata_map_and_write (sai_urb, landurban%settyp, typindex, m_urb2diag, &
                  -1.0e36_r8, landname, 'TREE_SAI', compress = 6, write_mode = 'one',  &
                  lastdimname = 'Itime', lastdimvalue = imonth, defval = 0._r8, create_mode = .false.)!first_call_SAI_urban)

            IF (first_call_SAI_urban) first_call_SAI_urban = .false.
#endif

#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
            CALL check_vector_data ('Urban Tree SAI '//trim(cmonth), sai_urb)
#endif
         ENDDO
         first_call_SAI_urban = .true.
      ENDDO


      ! ******* PCT_Water *******
      ! allocate and read grided water cover raw data
      IF (p_is_io) THEN

         CALL allocate_block_data (grid_pctw, flakeu)

         dir = trim(dir_rawdata) // trim(DEF_rawdata%urban_flake%dir)

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            fname = trim(DEF_rawdata%urban_flake%fname)//trim(c5year)
         ELSE
            fname = trim(DEF_rawdata%urban_flake%fname)//'.'//trim(c5year)
         ENDIF

         CALL read_5x5_data (dir, fname, grid_pctw, "PCT_Water", flakeu)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_pctw, flakeu)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (pct_water (numurban))

         pct_water (:) = 0.
         ! loop for urban patch to aggregate water cover data with area-weighted average
         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid_pctw, zip = USE_zip_for_aggregation, area = area_one, &
               data_r8_2d_in1 = flakeu, data_r8_2d_out1 = flakeu_one)

            WHERE (flakeu_one < 0)
               area_one = 0
            END WHERE
            ! only calculate when urban patch have water cover
            IF (sum(area_one) > 0) THEN
               pct_water(iurban) = sum(flakeu_one * area_one) / sum(area_one)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/PCT_Water.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'PCT_Water', 'urban', landurban, pct_water, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_ecology_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (pct_water, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'PCT_WATER', compress = 6, write_mode = 'one', defval = 0._r8, create_mode = .false.)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('Urban Water Cover ', pct_water)
#endif


      ! ******* LUCY_id *******
      ! allocate and read the LUCY id
      IF (p_is_io) THEN

         landname = TRIM(dir_rawdata)//trim(DEF_rawdata%urban_lucy%dir)//'LUCY_regionid.nc'

         CALL allocate_block_data (grid_lucy, LUCY_reg)
         CALL flush_block_data    (LUCY_reg , 0       )

         CALL ncio_read_block (landname, 'LUCY_REGION_ID', grid_lucy, LUCY_reg)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_lucy, data_i4_2d_in1 = LUCY_reg)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate ( LUCY_rid (numurban))

         LUCY_rid (:) = 0

         ! loop for each urban patch to get the LUCY id of all fine grid
         ! of iurban patch, then assign the most frequency id to this urban patch
         DO iurban = 1, numurban

            CALL aggregation_request_data (landurban, iurban, grid_lucy, &
               zip = USE_zip_for_aggregation, &
               data_i4_2d_in1 = LUCY_reg, data_i4_2d_out1 = LUCY_reg_one)
            ! the most frequency id to this urban patch
            LUCY_rid(iurban) = num_max_frequency (LUCY_reg_one, 0)
         ENDDO
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/LUCY_region_id.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'LUCY_id', 'urban', landurban, LUCY_rid, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_human_'//trim(cyear)//'.nc'

      IF (allocated(LUCY_rid)) THEN
         allocate (LUCY_rid_r8 (size(LUCY_rid)))
         LUCY_rid_r8 = real(LUCY_rid, r8)
      ENDIF

      CALL srfdata_map_and_write (LUCY_rid_r8, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'LUCY_REGION_ID', compress = 6, write_mode = 'one', create_mode = .true.)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('LUCY_ID ', LUCY_rid)
#endif


      ! ******* POP_DEN *******
      ! allocate and read the grided population raw data(500m)
      ! NOTE, the population is year-by-year
      IF (p_is_io) THEN

         CALL allocate_block_data (grid_pop, pop)

         dir = TRIM(dir_rawdata) // trim(DEF_rawdata%urban_pop%dir)

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            write(c5year, '(i4.4)') int(lc_year/5)*5
            fname = trim(DEF_rawdata%urban_pop%fname)//trim(c5year)

            ! population data is year by year,
            ! so pop_i is calculated to determine the dimension of POP data reads
            IF (mod(lc_year,5) == 0) THEN
               pop_i = 1
            ELSE
               pop_i = 5 - (ceiling(lc_year*1./5.)*5 - lc_year) + 1
            ENDIF

            ! read the population data of total 5x5 region
            CALL read_5x5_data_time (dir, fname, grid_pop, "POP_DEN", pop_i, pop)
         ELSE
            fname = trim(DEF_rawdata%urban_pop%fname)//'.'//trim(cyear)
            CALL read_5x5_data (dir, fname, grid_pop, "POP_density", pop)
         ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_pop, data_r8_2d_in1 = pop)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (pop_den (numurban))

         pop_den (:) = 0.

         ! loop for urban patch to aggregate population data with area-weighted average
         DO iurban = 1, numurban
            ! request all fine grid data and area of the iurban urban patch
            ! a one dimension vector will be returned
            CALL aggregation_request_data (landurban, iurban, grid_pop, zip = USE_zip_for_aggregation, &
               area = area_one, data_r8_2d_in1 = pop, data_r8_2d_out1 = pop_one)

            WHERE (pop_one <= 0)
               area_one = 0
            END WHERE
            ! area-weighted average
            IF (sum(area_one) > 0._r8) THEN
               pop_den(iurban) = sum(pop_one * area_one) / sum(area_one)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/POP.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'POP_DEN', 'urban', landurban, pop_den, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_human_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (pop_den, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'POPULATION_DEN', compress = 6, write_mode = 'one', defval = 0._r8, create_mode = .false.)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('POP_DEN ', pop_den)
#endif


IF (DEF_URBAN_type_scheme == 1) THEN
      ! look up table of NCAR urban properties (using look-up tables)
      landname = TRIM(dir_rawdata)//'/urban_physical/NCAR_urban_properties.nc'

      CALL ncio_read_bcast_serial (landname,  "EM_ROOF"       , emroof_ncar )
      CALL ncio_read_bcast_serial (landname,  "EM_WALL"       , emwall_ncar )
      CALL ncio_read_bcast_serial (landname,  "EM_IMPROAD"    , emgimp_ncar )
      CALL ncio_read_bcast_serial (landname,  "EM_PERROAD"    , emgper_ncar )
      CALL ncio_read_bcast_serial (landname,  "ALB_ROOF"      , albroof_ncar)
      CALL ncio_read_bcast_serial (landname,  "ALB_WALL"      , albwall_ncar)
      CALL ncio_read_bcast_serial (landname,  "ALB_IMPROAD"   , albgimp_ncar)
      CALL ncio_read_bcast_serial (landname,  "ALB_PERROAD"   , albgper_ncar)
      CALL ncio_read_bcast_serial (landname,  "TK_ROOF"       , tkroof_ncar )
      CALL ncio_read_bcast_serial (landname,  "TK_WALL"       , tkwall_ncar )
      CALL ncio_read_bcast_serial (landname,  "TK_IMPROAD"    , tkgimp_ncar )
      CALL ncio_read_bcast_serial (landname,  "CV_ROOF"       , cvroof_ncar )
      CALL ncio_read_bcast_serial (landname,  "CV_WALL"       , cvwall_ncar )
      CALL ncio_read_bcast_serial (landname,  "CV_IMPROAD"    , cvgimp_ncar )
      CALL ncio_read_bcast_serial (landname,  "THICK_ROOF"    , thkroof_ncar)
      CALL ncio_read_bcast_serial (landname,  "THICK_WALL"    , thkwall_ncar)
      CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MIN", tbldmin_ncar)
      CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MAX", tbldmax_ncar)
ENDIF

      IF (p_is_io) THEN

         CALL allocate_block_data (grid_alb, albroof)

         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%urban_alb%dir)
         fname = trim(DEF_rawdata%urban_alb%fname)
         CALL read_5x5_data (dir, fname, grid_alb, "ALB_ROOF", albroof)

#ifdef USEMPI
IF (DEF_URBAN_type_scheme == 1) THEN
         CALL aggregation_data_daemon_multigrid (grid_in1 = grid_urban, data_i4_2d_in1 = reg_typid, &
            grid_in2 = grid_alb, data_r8_2d_in1 = albroof)
ELSE
         CALL aggregation_data_daemon_multigrid (grid_in1 = grid_patch, &
            grid_in2 = grid_alb, data_r8_2d_in2 = albroof)
ENDIF
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (em_roof          (numurban))
         allocate (em_wall          (numurban))
         allocate (em_gimp          (numurban))
         allocate (em_gper          (numurban))
         allocate (thk_roof         (numurban))
         allocate (thk_wall         (numurban))
         allocate (tbld_min         (numurban))
         allocate (tbld_max         (numurban))
         allocate (tk_wgt     (ulev, numurban))
         allocate (cv_wgt     (ulev, numurban))
         allocate (cv_roof    (ulev, numurban))
         allocate (cv_wall    (ulev, numurban))
         allocate (cv_gimp    (ulev, numurban))
         allocate (tk_roof    (ulev, numurban))
         allocate (tk_wall    (ulev, numurban))
         allocate (tk_gimp    (ulev, numurban))
         allocate (alb_roof (nr, ns, numurban))
         allocate (alb_wall (nr, ns, numurban))
         allocate (alb_gimp (nr, ns, numurban))
         allocate (alb_gper (nr, ns, numurban))

         allocate (urb_pct          (numurban))
         allocate (sarea_urb        (numurban))
         allocate (area_urb         (numurban))

         ! initialization
         urb_pct  (:)     = 0.
         sarea_urb(:)     = 0.
         area_urb (:)     = 0.

         em_roof  (:)     = 0.
         em_wall  (:)     = 0.
         em_gimp  (:)     = 0.
         em_gper  (:)     = 0.
         thk_roof (:)     = 0.
         thk_wall (:)     = 0.
         tbld_min (:)     = 0.
         tbld_max (:)     = 0.
         tk_wgt   (:,:)   = 0.
         cv_wgt   (:,:)   = 0.
         cv_roof  (:,:)   = 0.
         cv_wall  (:,:)   = 0.
         cv_gimp  (:,:)   = 0.
         tk_roof  (:,:)   = 0.
         tk_wall  (:,:)   = 0.
         tk_gimp  (:,:)   = 0.
         alb_roof (:,:,:) = 0.
         alb_wall (:,:,:) = 0.
         alb_gimp (:,:,:) = 0.
         alb_gper (:,:,:) = 0.

         ! loop for each urban patch to aggregate urban morphological and thermal paras with area-weighted average
         DO iurban = 1, numurban
            ! urban region and type id for look-up-table
            urb_typidx = landurban%settyp(iurban)

IF (DEF_URBAN_type_scheme == 1) THEN
            CALL aggregation_request_data_multigrid (landurban, iurban, &
               grid_in1 = grid_urban, area = area_one, data_i4_2d_in1 = reg_typid, data_i4_2d_out1 = reg_typid_one, &
               grid_in2 = grid_alb  , data_r8_2d_in2 = albroof, data_r8_2d_out2 = alb_roof_one_)
ELSE
            CALL aggregation_request_data_multigrid (landurban, iurban, grid_in1 = grid_patch, area = area_one, &
               grid_in2 = grid_alb, data_r8_2d_in2 = albroof, data_r8_2d_out2 = alb_roof_one_)
ENDIF

            sumarea          = sum(area_one)
            area_urb(iurban) = sumarea

            IF (.not. allocated (alb_roof_one)) THEN
               allocate (alb_roof_one (nr, ns, size(area_one)))
               allocate (alb_wall_one (nr, ns, size(area_one)))
               allocate (alb_gimp_one (nr, ns, size(area_one)))
               allocate (alb_gper_one (nr, ns, size(area_one)))
            ELSE
               deallocate (alb_roof_one)
               deallocate (alb_wall_one)
               deallocate (alb_gimp_one)
               deallocate (alb_gper_one)

               allocate (alb_roof_one (nr, ns, size(area_one)))
               allocate (alb_wall_one (nr, ns, size(area_one)))
               allocate (alb_gimp_one (nr, ns, size(area_one)))
               allocate (alb_gper_one (nr, ns, size(area_one)))
            ENDIF

            ! loop for each finer grid to aggregate data
            DO ipxl = 1, size(area_one)

IF (DEF_URBAN_type_scheme == 1) THEN
               ! same for above, assign reg id for RG_-45_65_-50_70
               IF (all(reg_typid_one==0)) THEN
                  reg_typid_one(:) = 30
               ENDIF

               IF (any(reg_typid_one==0)) THEN
                  WHERE(reg_typid_one==0) reg_typid_one = num_max_frequency(reg_typid_one)
               ENDIF

               urb_regidx = reg_typid_one(ipxl)

               emroof  = emroof_ncar (urb_typidx,urb_regidx)
               emwall  = emwall_ncar (urb_typidx,urb_regidx)
               emgper  = emgper_ncar (urb_typidx,urb_regidx)
               emgimp  = emgimp_ncar (urb_typidx,urb_regidx)

               thkroof = thkroof_ncar(urb_typidx,urb_regidx)
               thkwall = thkwall_ncar(urb_typidx,urb_regidx)
               tbldmax = tbldmax_ncar(urb_typidx,urb_regidx)
               tbldmin = tbldmin_ncar(urb_typidx,urb_regidx)

               cvroof   (:) = cvroof_ncar (urb_typidx,urb_regidx,:)
               cvwall   (:) = cvwall_ncar (urb_typidx,urb_regidx,:)
               cvgimp   (:) = cvgimp_ncar (urb_typidx,urb_regidx,:)

               tkroof   (:) = tkroof_ncar (urb_typidx,urb_regidx,:)
               tkwall   (:) = tkwall_ncar (urb_typidx,urb_regidx,:)
               tkgimp   (:) = tkgimp_ncar (urb_typidx,urb_regidx,:)

               IF (alb_roof_one_(ipxl) > 0) THEN
                  alb_roof_one(:,:,ipxl) = alb_roof_one_(ipxl)
               ELSE
                  alb_roof_one(:,:,ipxl) = albroof_ncar(urb_typidx,urb_regidx,:,:)
               ENDIF

               alb_wall_one(:,:,ipxl) = albwall_ncar(urb_typidx,urb_regidx,:,:)
               alb_gper_one(:,:,ipxl) = albgper_ncar(urb_typidx,urb_regidx,:,:)
               alb_gimp_one(:,:,ipxl) = albgimp_ncar(urb_typidx,urb_regidx,:,:)

ELSE IF (DEF_URBAN_type_scheme == 2) THEN

               emroof  = emroof_lcz (urb_typidx)
               emwall  = emwall_lcz (urb_typidx)
               emgper  = emgper_lcz (urb_typidx)
               emgimp  = emgimp_lcz (urb_typidx)

               thkroof = thkroof_lcz(urb_typidx)
               thkwall = thkwall_lcz(urb_typidx)
               tbldmax = tbldmax_lcz(urb_typidx)
               tbldmin = tbldmin_lcz(urb_typidx)

               cvroof   (:) = cvroof_lcz (urb_typidx)
               cvwall   (:) = cvwall_lcz (urb_typidx)
               cvgimp   (:) = cvgimp_lcz (urb_typidx)

               tkroof   (:) = tkroof_lcz (urb_typidx)
               tkwall   (:) = tkwall_lcz (urb_typidx)
               tkgimp   (:) = tkgimp_lcz (urb_typidx)

               IF (alb_roof_one_(ipxl) > 0) THEN
                  alb_roof_one(:,:,ipxl) = alb_roof_one_(ipxl)
               ELSE
                  alb_roof_one(:,:,ipxl) = albroof_lcz(urb_typidx)
               ENDIF

               alb_wall_one(:,:,ipxl) = albwall_lcz(urb_typidx)
               alb_gper_one(:,:,ipxl) = albgper_lcz(urb_typidx)
               alb_gimp_one(:,:,ipxl) = albgimp_lcz(urb_typidx)
ENDIF

               em_roof (iurban) = em_roof  (iurban) + emroof  * area_one(ipxl)
               em_wall (iurban) = em_wall  (iurban) + emwall  * area_one(ipxl)
               em_gimp (iurban) = em_gimp  (iurban) + emgimp  * area_one(ipxl)
               em_gper (iurban) = em_gper  (iurban) + emgper  * area_one(ipxl)

               thk_roof(iurban) = thk_roof (iurban) + thkroof * area_one(ipxl)
               thk_wall(iurban) = thk_wall (iurban) + thkwall * area_one(ipxl)
               tbld_min(iurban) = tbld_min (iurban) + tbldmin * area_one(ipxl)
               tbld_max(iurban) = tbld_max (iurban) + tbldmax * area_one(ipxl)

               ! tkgimp and cvgimp may have nan-values, and need to be calculated separately
               DO il = 1, 10
                  IF (tkgimp(il) .ne. -999.) THEN
                     tk_gimp(il,iurban) = tk_gimp(il,iurban) + tkgimp(il) * area_one(ipxl)
                     tk_wgt (il,iurban) = tk_wgt (il,iurban) + area_one(ipxl)
                  ENDIF

                  IF (cvgimp(il) .ne. -999.) THEN
                     cv_gimp(il,iurban) = cv_gimp(il,iurban) + cvgimp(il) * area_one(ipxl)
                     cv_wgt (il,iurban) = cv_wgt (il,iurban) + area_one(ipxl)
                  ENDIF
               ENDDO

               cv_roof (:,iurban) = cv_roof (:,iurban) + cvroof(:) * area_one(ipxl)
               cv_wall (:,iurban) = cv_wall (:,iurban) + cvwall(:) * area_one(ipxl)
               tk_roof (:,iurban) = tk_roof (:,iurban) + tkroof(:) * area_one(ipxl)
               tk_wall (:,iurban) = tk_wall (:,iurban) + tkwall(:) * area_one(ipxl)

               alb_roof(:,:,iurban) = alb_roof(:,:,iurban) + alb_roof_one(:,:,ipxl) * area_one(ipxl)
               alb_wall(:,:,iurban) = alb_wall(:,:,iurban) + alb_wall_one(:,:,ipxl) * area_one(ipxl)
               alb_gimp(:,:,iurban) = alb_gimp(:,:,iurban) + alb_gimp_one(:,:,ipxl) * area_one(ipxl)
               alb_gper(:,:,iurban) = alb_gper(:,:,iurban) + alb_gper_one(:,:,ipxl) * area_one(ipxl)

            ENDDO

            em_roof  (iurban) = em_roof   (iurban) / sumarea
            em_wall  (iurban) = em_wall   (iurban) / sumarea
            em_gimp  (iurban) = em_gimp   (iurban) / sumarea
            em_gper  (iurban) = em_gper   (iurban) / sumarea
            thk_roof (iurban) = thk_roof  (iurban) / sumarea
            thk_wall (iurban) = thk_wall  (iurban) / sumarea
            tbld_min (iurban) = tbld_min  (iurban) / sumarea
            tbld_max (iurban) = tbld_max  (iurban) / sumarea

            cv_roof(:,iurban) = cv_roof (:,iurban) / sumarea
            cv_wall(:,iurban) = cv_wall (:,iurban) / sumarea
            tk_roof(:,iurban) = tk_roof (:,iurban) / sumarea
            tk_wall(:,iurban) = tk_wall (:,iurban) / sumarea

            DO il = 1, 10
               IF (tk_wgt(il,iurban) > 0.) THEN
                  tk_gimp(il,iurban) = tk_gimp(il,iurban) / tk_wgt(il,iurban)
               ENDIF

               IF (cv_wgt(il,iurban) > 0.) THEN
                  cv_gimp(il,iurban) = cv_gimp(il,iurban) / cv_wgt(il,iurban)
               ENDIF
            ENDDO

            alb_roof(:,:,iurban) = alb_roof(:,:,iurban) / sumarea
            alb_wall(:,:,iurban) = alb_wall(:,:,iurban) / sumarea
            alb_gimp(:,:,iurban) = alb_gimp(:,:,iurban) / sumarea
            alb_gper(:,:,iurban) = alb_gper(:,:,iurban) / sumarea

         ENDDO

         IF (numelm>0 .and. numurban>0) THEN
            DO ielm = 1, numelm
               numpth = count(landurban%eindex==landelm%eindex(ielm))

               IF (allocated(locpth)) deallocate(locpth)
               allocate(locpth(numpth))

               locpth = pack([(ipth, ipth=1, numurban)], &
                        landurban%eindex==landelm%eindex(ielm))

               urb_s = minval(locpth)
               urb_e = maxval(locpth)

               DO iurban = urb_s, urb_e
                  sarea_urb(urb_s:urb_e) = sarea_urb(urb_s:urb_e) + area_urb(iurban)
               ENDDO
            ENDDO

            urb_pct(:) = area_urb(:)/sarea_urb(:)
         ENDIF

#ifdef USEMPI
         CALL aggregation_worker_done_multigrid ()
#endif
      ENDIF

      !output
      write(cyear,'(i4.4)') lc_year
      landname = trim(dir_model_landdata) // '/urban/'//trim(cyear)//'/urban.nc'
      CALL ncio_create_file_vector (landname, landurban)

      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_define_dimension_vector (landname, landurban, 'numsolar', ns)
      CALL ncio_define_dimension_vector (landname, landurban, 'numrad'  , nr)
      CALL ncio_define_dimension_vector (landname, landurban, 'ulev'    , ulev)

      CALL ncio_write_vector (landname, 'EM_ROOF'       , 'urban', landurban, em_roof , DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'EM_WALL'       , 'urban', landurban, em_wall , DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'EM_IMPROAD'    , 'urban', landurban, em_gimp , DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'EM_PERROAD'    , 'urban', landurban, em_gper , DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'THICK_ROOF'    , 'urban', landurban, thk_roof, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'THICK_WALL'    , 'urban', landurban, thk_wall, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'T_BUILDING_MIN', 'urban', landurban, tbld_min, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'T_BUILDING_MAX', 'urban', landurban, tbld_max, DEF_Srfdata_CompressLevel)

      CALL ncio_write_vector (landname, 'CV_ROOF'   , 'ulev', ulev, 'urban', landurban, cv_roof, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'CV_WALL'   , 'ulev', ulev, 'urban', landurban, cv_wall, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'TK_ROOF'   , 'ulev', ulev, 'urban', landurban, tk_roof, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'TK_WALL'   , 'ulev', ulev, 'urban', landurban, tk_wall, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'TK_IMPROAD', 'ulev', ulev, 'urban', landurban, tk_gimp, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'CV_IMPROAD', 'ulev', ulev, 'urban', landurban, cv_gimp, DEF_Srfdata_CompressLevel)

      CALL ncio_write_vector (landname, 'ALB_ROOF'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_roof, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'ALB_WALL'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_wall, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'ALB_IMPROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_gimp, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (landname, 'ALB_PERROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_gper, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_type_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (urb_pct, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'PCT_TYPE', compress = 6, write_mode = 'one', &
         stat_mode = 'fraction', defval = 0._r8, pctshared = urb_pct, create_mode = .true.)

      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_model_landdata) // '/diag/urban_phyical_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (em_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'EM_ROOF', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (em_wall, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'EM_WALL', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (em_gper, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'EM_PERROAD', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (em_gimp, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'EM_IMPROAD', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (thk_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'THICK_ROOF', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (thk_wall, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'THICK_WALL', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (tbld_min, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'T_BUILDING_MIN', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (tbld_max, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'T_BUILDING_MAX', compress = 6, write_mode = 'one')

      DO il = 1, nl_roof
         CALL srfdata_map_and_write (cv_roof(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'CV_ROOF', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)

         CALL srfdata_map_and_write (tk_roof(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'TK_ROOF', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)
      ENDDO

      DO il = 1, nl_wall
         CALL srfdata_map_and_write (cv_wall(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'CV_WALL', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)

         CALL srfdata_map_and_write (tk_wall(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'TK_WALL', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)
      ENDDO

      DO il = 1, nl_soil
         CALL srfdata_map_and_write (cv_gimp(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'CV_IMPROAD', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)

         CALL srfdata_map_and_write (tk_gimp(il,:), landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'TK_IMPROAD', compress = 6, write_mode = 'one', &
            lastdimname = 'ulev', lastdimvalue = il)
      ENDDO

      CALL srfdata_map_and_write (alb_roof(1,1,:), landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'ALB_ROOF', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (alb_wall(1,1,:), landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'ALB_WALL', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (alb_gper(1,1,:), landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'ALB_PERROAD', compress = 6, write_mode = 'one')

      CALL srfdata_map_and_write (alb_gimp(1,1,:), landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'ALB_IMPROAD', compress = 6, write_mode = 'one')

      deallocate (typindex)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('EM_ROOF '       , em_roof )
      CALL check_vector_data ('EM_WALL '       , em_wall )
      CALL check_vector_data ('EM_IMPROAD '    , em_gimp )
      CALL check_vector_data ('EM_PERROAD '    , em_gper )
      CALL check_vector_data ('ALB_ROOF '      , alb_roof)
      CALL check_vector_data ('ALB_WALL '      , alb_wall)
      CALL check_vector_data ('ALB_IMPROAD '   , alb_gimp)
      CALL check_vector_data ('ALB_PERROAD '   , alb_gper)
      CALL check_vector_data ('TK_ROOF '       , tk_roof )
      CALL check_vector_data ('TK_WALL '       , tk_wall )
      CALL check_vector_data ('TK_IMPROAD '    , tk_gimp )
      CALL check_vector_data ('CV_ROOF '       , cv_roof )
      CALL check_vector_data ('CV_WALL '       , cv_wall )
      CALL check_vector_data ('CV_IMPROAD '    , cv_gimp )
      CALL check_vector_data ('THICK_ROOF '    , thk_roof)
      CALL check_vector_data ('THICK_WALL '    , thk_wall)
      CALL check_vector_data ('T_BUILDING_MIN ', tbld_min)
      CALL check_vector_data ('T_BUILDING_MAX ', tbld_max)
#endif

      IF (p_is_worker) THEN

         IF ( allocated (LUCY_rid ) ) deallocate (LUCY_rid  )
         IF ( allocated (pop_den  ) ) deallocate (pop_den   )
         IF ( allocated (pct_tree ) ) deallocate (pct_tree  )
         IF ( allocated (htop_urb ) ) deallocate (htop_urb  )
         IF ( allocated (pct_water) ) deallocate (pct_water )
         IF ( allocated (pct_roof ) ) deallocate (pct_roof  )
         IF ( allocated (ht_roof  ) ) deallocate (ht_roof   )
         IF ( allocated (lai_urb  ) ) deallocate (lai_urb   )
         IF ( allocated (sai_urb  ) ) deallocate (sai_urb   )
         IF ( allocated (area_urb ) ) deallocate (area_urb  )
         IF ( allocated (sarea_urb) ) deallocate (sarea_urb )
         IF ( allocated (urb_pct  ) ) deallocate (urb_pct   )

IF (DEF_URBAN_type_scheme == 1) THEN
         IF ( allocated (htroof_ncar ) ) deallocate (htroof_ncar )
         IF ( allocated (wtroof_ncar ) ) deallocate (wtroof_ncar )
         IF ( allocated (hwrbld_ncar ) ) deallocate (hwrbld_ncar )
         IF ( allocated (fgper_ncar  ) ) deallocate (fgper_ncar  )
         IF ( allocated (emroof_ncar ) ) deallocate (emroof_ncar )
         IF ( allocated (emwall_ncar ) ) deallocate (emwall_ncar )
         IF ( allocated (emgper_ncar ) ) deallocate (emgper_ncar )
         IF ( allocated (emgimp_ncar ) ) deallocate (emgimp_ncar )
         IF ( allocated (albroof_ncar) ) deallocate (albroof_ncar)
         IF ( allocated (albwall_ncar) ) deallocate (albwall_ncar)
         IF ( allocated (albgper_ncar) ) deallocate (albgper_ncar)
         IF ( allocated (albgimp_ncar) ) deallocate (albgimp_ncar)
         IF ( allocated (thkroof_ncar) ) deallocate (thkroof_ncar)
         IF ( allocated (thkwall_ncar) ) deallocate (thkwall_ncar)
         IF ( allocated (tbldmin_ncar) ) deallocate (tbldmin_ncar)
         IF ( allocated (tbldmax_ncar) ) deallocate (tbldmax_ncar)
         IF ( allocated (cvroof_ncar ) ) deallocate (cvroof_ncar )
         IF ( allocated (cvwall_ncar ) ) deallocate (cvwall_ncar )
         IF ( allocated (cvgimp_ncar ) ) deallocate (cvgimp_ncar )
         IF ( allocated (tkroof_ncar ) ) deallocate (tkroof_ncar )
         IF ( allocated (tkwall_ncar ) ) deallocate (tkwall_ncar )
         IF ( allocated (tkgimp_ncar ) ) deallocate (tkgimp_ncar )
ENDIF

         IF ( allocated (hlr_bld ) ) deallocate (hlr_bld )
         IF ( allocated (f_gper  ) ) deallocate (f_gper  )
         IF ( allocated (f_gimp  ) ) deallocate (f_gimp  )
         IF ( allocated (f_isa   ) ) deallocate (f_isa   )
         IF ( allocated (em_roof ) ) deallocate (em_roof )
         IF ( allocated (em_wall ) ) deallocate (em_wall )
         IF ( allocated (em_gimp ) ) deallocate (em_gimp )
         IF ( allocated (em_gper ) ) deallocate (em_gper )
         IF ( allocated (thk_roof) ) deallocate (thk_roof)
         IF ( allocated (thk_wall) ) deallocate (thk_wall)
         IF ( allocated (tbld_min) ) deallocate (tbld_min)
         IF ( allocated (tbld_max) ) deallocate (tbld_max)
         IF ( allocated (tk_wgt  ) ) deallocate (tk_wgt  )
         IF ( allocated (cv_wgt  ) ) deallocate (cv_wgt  )
         IF ( allocated (cv_roof ) ) deallocate (cv_roof )
         IF ( allocated (cv_wall ) ) deallocate (cv_wall )
         IF ( allocated (cv_gimp ) ) deallocate (cv_gimp )
         IF ( allocated (cv_wgt  ) ) deallocate (cv_wgt  )
         IF ( allocated (tk_roof ) ) deallocate (tk_roof )
         IF ( allocated (tk_wall ) ) deallocate (tk_wall )
         IF ( allocated (tk_gimp ) ) deallocate (tk_gimp )
         IF ( allocated (tk_wgt  ) ) deallocate (tk_wgt  )
         IF ( allocated (alb_roof) ) deallocate (alb_roof)
         IF ( allocated (alb_wall) ) deallocate (alb_wall)
         IF ( allocated (alb_gimp) ) deallocate (alb_gimp)
         IF ( allocated (alb_gper) ) deallocate (alb_gper)

         IF ( allocated (area_one     ) ) deallocate (area_one     )
         IF ( allocated (LUCY_reg_one ) ) deallocate (LUCY_reg_one )
         IF ( allocated (pop_one      ) ) deallocate (pop_one      )
         IF ( allocated (fvegu_one    ) ) deallocate (fvegu_one    )
         IF ( allocated (htopu_one    ) ) deallocate (htopu_one    )
         IF ( allocated (flakeu_one   ) ) deallocate (flakeu_one   )
         IF ( allocated (wt_roof_one  ) ) deallocate (wt_roof_one  )
         IF ( allocated (ht_roof_one  ) ) deallocate (ht_roof_one  )
         IF ( allocated (hlr_bld_one  ) ) deallocate (hlr_bld_one  )
         IF ( allocated (f_gper_one   ) ) deallocate (f_gper_one   )
         IF ( allocated (ulai_one     ) ) deallocate (ulai_one     )
         IF ( allocated (usai_one     ) ) deallocate (usai_one     )

      ENDIF

END SUBROUTINE Aggregation_Urban
#endif
