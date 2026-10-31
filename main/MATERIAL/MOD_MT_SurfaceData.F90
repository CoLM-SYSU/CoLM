#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_SurfaceData
   USE MOD_Block, only: get_filename_block
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_clip_methane_surface_data
CONTAINS
   SUBROUTINE mt_clip_methane_surface_data(dir_landdata_in, dir_landdata_out, input_year, lc_cyear, &
      iblk, jblk, patchmask, requires_lake_soilc, requires_spatial_ph, clip_vector, clip_error)
      character(len=*), intent(in) :: dir_landdata_in, dir_landdata_out, input_year, lc_cyear
      integer, intent(in) :: iblk, jblk
      logical, intent(in) :: patchmask(:), requires_lake_soilc, requires_spatial_ph
      character(len=256) :: file_in, file_out, fileblock
      logical :: fexists
      interface
         SUBROUTINE clip_vector(file_in, file_out, iblk, jblk, varname, vecmask)
            character(len=*), intent(in) :: file_in, file_out, varname
            integer, intent(in) :: iblk, jblk
            logical, intent(in) :: vecmask(:)
         END SUBROUTINE clip_vector
         SUBROUTINE clip_error(message)
            character(len=*), intent(in) :: message
         END SUBROUTINE clip_error
      end interface

               ! Methane surface data (lake-bottom soil carbon, soil pH), made by
               ! mksrfdata when methane needs them; required inputs must exist.

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/lake_soilc_patches.nc'
               CALL get_filename_block (file_in, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (.not. fexists .and. requires_lake_soilc) CALL clip_error('Missing required lake carbon: '//trim(fileblock))
               IF (fexists) THEN
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/soil/' // lc_cyear)
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/lake_soilc_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'lake_soilc_patches', patchmask)
               ENDIF

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/methane_ph_patches.nc'
               CALL get_filename_block (file_in, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (.not. fexists .and. requires_spatial_ph) CALL clip_error('Missing required methane pH: '//trim(fileblock))
               IF (fexists) THEN
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/soil/' // lc_cyear)
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/methane_ph_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'methane_ph_patches', patchmask)
               ENDIF
   END SUBROUTINE mt_clip_methane_surface_data
END MODULE MOD_MT_SurfaceData
#endif
