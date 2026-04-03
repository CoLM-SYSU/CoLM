#include <define.h>

#ifdef SinglePoint
MODULE MOD_HistSingle

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out model results at sites to history files.
!
!  Created by Shupeng Zhang, July 2023
!
!  TODO...(need complement)
!----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_Namelist, only: USE_SITE_HistWriteBack
   USE MOD_SPMD_Task

   logical :: memory_to_disk

   integer :: ntime_mem
   integer :: itime_mem
   integer, allocatable :: time_memory (:)

   ! -- data type --
   type hist_memory_type
      character(len=256) :: varname
      real(r8), allocatable :: v2d (:,:)
      real(r8), allocatable :: v3d (:,:,:)
      real(r8), allocatable :: v4d (:,:,:,:)
      type(hist_memory_type), pointer :: next
   END type hist_memory_type

   type(hist_memory_type), target  :: hist_memory
   type(hist_memory_type), pointer :: thisvar, nextvar

CONTAINS

   ! -- initialize history IO --
   SUBROUTINE hist_single_init ()

      USE MOD_Namelist
      IMPLICIT NONE

      ! Local Variables
      real(r8) :: secs_group, secs_write

      IF (USE_SITE_HistWriteBack) THEN

         IF ( trim(DEF_HIST_groupby) == 'YEAR' ) THEN
            secs_group = 366*24*3600
         ELSEIF ( trim(DEF_HIST_groupby) == 'MONTH' ) THEN
            secs_group = 31*24*3600
         ELSEIF ( trim(DEF_HIST_groupby) == 'DAY' ) THEN
            secs_group = 24*3600
         ENDIF

         select CASE (trim(adjustl(DEF_HIST_FREQ)))
         CASE ('TIMESTEP')
            secs_write = DEF_simulation_time%timestep
         CASE ('HOURLY')
            secs_write = 3600
         CASE ('DAILY')
            secs_write = 24*3600
         CASE ('MONTHLY')
            secs_write = 31*24*3600
         CASE ('YEARLY')
            secs_write = 366*24*3600
         END select

         ntime_mem = ceiling(secs_group / secs_write) + 2

         allocate (time_memory (ntime_mem))

         itime_mem = 0

         hist_memory%next => null()
         hist_memory%varname = ''

         thisvar => hist_memory

      ENDIF

   END SUBROUTINE hist_single_init

   ! -- finalize history IO --
   SUBROUTINE hist_single_final ()

      IMPLICIT NONE

      IF (USE_SITE_HistWriteBack) THEN

         thisvar => hist_memory%next
         DO WHILE (associated(thisvar))
            nextvar => thisvar%next
            IF (allocated(thisvar%v2d)) deallocate(thisvar%v2d)
            IF (allocated(thisvar%v3d)) deallocate(thisvar%v3d)
            IF (allocated(thisvar%v4d)) deallocate(thisvar%v4d)
            deallocate(thisvar)
            thisvar => nextvar
         ENDDO

         deallocate (time_memory)

      ENDIF

   END SUBROUTINE hist_single_final

   ! -- write history time --
   SUBROUTINE hist_single_write_time (filename, filelast, dataname, time, itime)

      USE MOD_Namelist
      USE MOD_TimeManager
      USE MOD_SingleSrfData
      USE MOD_NetCDFSerial
      USE MOD_Landpatch, only: numpatch
#ifdef URBAN_MODEL
      USE MOD_Landurban, only: numurban
#endif
      IMPLICIT NONE

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: filelast
      character (len=*), intent(in) :: dataname

      integer, intent(in)  :: time(3)
      integer, intent(out) :: itime

      ! Local variables
      integer :: minutes
      logical :: fexists

      inquire (file=filename, exist=fexists)
      IF ((.not. fexists) .or. (trim(filename) /= trim(filelast))) THEN
         CALL ncio_create_file (trim(filename))
         CALL ncio_define_dimension(filename, 'patch', numpatch)

         CALL ncio_write_serial (filename, 'lat', SITE_lat_location)
         CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
         CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

         CALL ncio_write_serial (filename, 'lon', SITE_lon_location)
         CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
         CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

         CALL ncio_write_colm_dimension (filename)

         IF (.not. USE_SITE_HistWriteBack) THEN
            CALL ncio_define_dimension(filename, 'time', 0)
         ENDIF

      ENDIF

      IF (USE_SITE_HistWriteBack) THEN

         minutes = minutes_since_1900 (time(1), time(2), time(3))

         select CASE (trim(adjustl(DEF_HIST_FREQ)))
         CASE ('HOURLY')
            minutes = minutes - 30
         CASE ('DAILY')
            minutes = minutes - 720
         CASE ('MONTHLY')
            minutes = minutes - 21600
         CASE ('YEARLY')
            minutes = minutes - 262800
         END select

         itime_mem = itime_mem + 1
         time_memory(itime_mem) = minutes

         IF (memory_to_disk) THEN
            CALL ncio_define_dimension(filename, 'time', itime_mem)
            CALL ncio_write_serial (filename, dataname, time_memory(1:itime_mem), 'time')
            CALL ncio_put_attr (filename, dataname, 'long_name', 'time')
            CALL ncio_put_attr (filename, dataname, 'units', 'minutes since 1900-1-1 0:0:0')
         ENDIF

         thisvar => hist_memory

      ELSE
         CALL ncio_write_time (filename, dataname, time, itime, DEF_HIST_FREQ)
      ENDIF

   END SUBROUTINE hist_single_write_time

   ! -- write 2D data --
   SUBROUTINE single_write_2d ( &
         acc_vec, file_hist, varname, itime_in_file, longname, units)

      USE MOD_Vars_Global,      only: spval
      IMPLICIT NONE

      real(r8),         intent(inout) :: acc_vec(:)
      character(len=*), intent(in)    :: file_hist
      character(len=*), intent(in)    :: varname
      integer,          intent(in)    :: itime_in_file
      character(len=*), intent(in)    :: longname
      character(len=*), intent(in)    :: units

      IF (USE_SITE_HistWriteBack) THEN

         IF (.not. associated(thisvar%next)) THEN
            allocate (thisvar%next)
            thisvar => thisvar%next

            thisvar%next    => null()
            thisvar%varname = varname
            allocate(thisvar%v2d (size(acc_vec),ntime_mem))
         ELSE
            thisvar => thisvar%next
         ENDIF

         IF (thisvar%varname /= varname) THEN
            write(*,*) 'Warning: history variable in memory is wrong: ' &
               // trim(thisvar%varname) // ' should be ' // trim(varname)
            CALL CoLM_stop ()
         ENDIF

         thisvar%v2d(:,itime_mem) = acc_vec(:)

         IF (memory_to_disk) THEN
            CALL ncio_write_serial (file_hist, varname, &
               thisvar%v2d(:,1:itime_mem), 'patch', 'time')
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ELSE
         CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            'patch', 'time')
         IF (itime_in_file == 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF
      ENDIF

   END SUBROUTINE single_write_2d

   ! -- write 3D data --
   SUBROUTINE single_write_3d ( &
         acc_vec, file_hist, varname, itime_in_file, dim1name, ndim1, longname, units)

      USE MOD_Vars_1DAccFluxes, only: nac
      USE MOD_Vars_Global,      only: spval
      IMPLICIT NONE

      real(r8),         intent(inout) :: acc_vec(:,:)
      character(len=*), intent(in)    :: file_hist
      character(len=*), intent(in)    :: varname
      integer,          intent(in)    :: itime_in_file
      character(len=*), intent(in)    :: dim1name
      integer,          intent(in)    :: ndim1
      character(len=*), intent(in)    :: longname
      character(len=*), intent(in)    :: units

      WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac

      IF (USE_SITE_HistWriteBack) THEN

         IF (.not. associated(thisvar%next)) THEN
            allocate (thisvar%next)
            thisvar => thisvar%next

            thisvar%next    => null()
            thisvar%varname = varname
            allocate(thisvar%v3d (ndim1,size(acc_vec,2),ntime_mem))
         ELSE
            thisvar => thisvar%next
         ENDIF

         IF (thisvar%varname /= varname) THEN
            write(*,*) 'Warning: history variable in memory is wrong: ' &
               // trim(thisvar%varname) // ' should be ' // trim(varname)
            CALL CoLM_stop ()
         ENDIF

         thisvar%v3d(:,:,itime_mem) = acc_vec

         IF (memory_to_disk) THEN
            CALL ncio_define_dimension (file_hist, dim1name, ndim1)
            CALL ncio_write_serial (file_hist, varname, thisvar%v3d(:,:,1:itime_mem), &
               dim1name, 'patch', 'time')
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ELSE
         CALL ncio_define_dimension (file_hist, dim1name, ndim1)
         CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, 'patch', 'time')
         IF (itime_in_file == 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF
      ENDIF

   END SUBROUTINE single_write_3d

   ! -- write 4D data --
   SUBROUTINE single_write_4d ( &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, ndim1, dim2name, ndim2, longname, units)

      USE MOD_Vars_Global,      only: spval
      IMPLICIT NONE

      real(r8),         intent(inout) :: acc_vec(:,:,:)
      character(len=*), intent(in)    :: file_hist
      character(len=*), intent(in)    :: varname
      integer,          intent(in)    :: itime_in_file
      character(len=*), intent(in)    :: dim1name
      integer,          intent(in)    :: ndim1
      character(len=*), intent(in)    :: dim2name
      integer,          intent(in)    :: ndim2
      character(len=*), intent(in)    :: longname
      character(len=*), intent(in)    :: units

      IF (USE_SITE_HistWriteBack) THEN

         IF (.not. associated(thisvar%next)) THEN
            allocate (thisvar%next)
            thisvar => thisvar%next

            thisvar%next    => null()
            thisvar%varname = varname
            allocate(thisvar%v4d (ndim1,ndim2,size(acc_vec,3),ntime_mem))
         ELSE
            thisvar => thisvar%next
         ENDIF

         IF (thisvar%varname /= varname) THEN
            write(*,*) 'Warning: history variable in memory is wrong: ' &
               // trim(thisvar%varname) // ' should be ' // trim(varname)
            CALL CoLM_stop ()
         ENDIF

         thisvar%v4d(:,:,:,itime_mem) = acc_vec

         IF (memory_to_disk) THEN
            CALL ncio_define_dimension (file_hist, dim1name, ndim1)
            CALL ncio_define_dimension (file_hist, dim2name, ndim2)
            CALL ncio_write_serial (file_hist, varname, thisvar%v4d(:,:,:,1:itime_mem), &
               dim1name, dim2name, 'patch', 'time')
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ELSE
         CALL ncio_define_dimension (file_hist, dim1name, ndim1)
         CALL ncio_define_dimension (file_hist, dim2name, ndim2)
         CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, dim2name, 'patch', 'time')
         IF (itime_in_file == 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF
      ENDIF

   END SUBROUTINE single_write_4d

END MODULE MOD_HistSingle
#endif
