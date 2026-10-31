#!/usr/bin/env python3
"""Check caller gross interception and production canopy isotope mixing."""
from pathlib import Path
import re
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
caller = (root / 'main/MATERIAL/MOD_MT_LandPhase.F90').read_text()
lines = re.findall(r'^\s*work%gross_(?:rain|snow) = .*$', caller, re.M)
assert len(lines) == 2, 'production gross interception expressions missing'
gross = '\n'.join(lines).replace('work%gross_rain','gross_intr_rain').replace('work%gross_snow','gross_intr_snow').replace('canopy_drain_rain','mt_leaf_qdrp_rain').replace('canopy_drain_snow','mt_leaf_qdrp_snow')
text = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Precip.F90').read_text()
precip = re.search(r'^   SUBROUTINE mt_precip \(.*?^   END SUBROUTINE mt_precip', text, re.M | re.S)[0]
source = '''module fixture
implicit none
integer, parameter :: r8=kind(1d0), nmts=1
real(r8), parameter :: mt_tiny=1d-20, mt_water_min_for_ratio=1d-12
logical :: DEF_VEG_SNOW=.true.
real(r8) :: mt_ldew_rain(1,1),mt_ldew_snow(1,1),a_mt_precip(1,1)
real(r8) :: mt_pg_rain_ground(1,1),mt_pg_snow_ground(1,1)
real(r8) :: mt_canopy_solid(1,1),mt_surface_solid(1,1),mt_waterstorage_solid(1,1)
real(r8),allocatable :: mt_waterstorage(:,:)
contains
logical function mt_uses_generic_water(i)
integer :: i
mt_uses_generic_water=.true.
end function
real(r8) function material_forcing_rain_value(i,p)
integer :: i,p
material_forcing_rain_value=5d0
end function
real(r8) function material_forcing_snow_value(i,p)
integer :: i,p
material_forcing_snow_value=7d0
end function
subroutine mt_equilibrate_dissolved(i,w,m,s)
integer :: i
real(r8) :: w,m,s
! This fixture is an isotope: it has no solubility limit.
end subroutine
''' + precip + '''
subroutine handoff(qintr_rain,qintr_snow,mt_leaf_qdrp_rain,mt_leaf_qdrp_snow,gross_intr_rain,gross_intr_snow)
real(r8),intent(in) :: qintr_rain,qintr_snow,mt_leaf_qdrp_rain,mt_leaf_qdrp_snow
real(r8),intent(out) :: gross_intr_rain,gross_intr_snow
''' + gross + '''
end subroutine
end module
program test
use fixture
implicit none
real(r8) :: gr,gs,old_snow,input_snow,drip,expected_ratio
integer :: c
do c=1,3
 old_snow=10d0
 drip=2d0
 input_snow=real(c-1,r8)
 call handoff(0d0,input_snow-drip,0d0,drip,gr,gs)
 call expect(gr,0d0,'rain gross')
 call expect(gs,input_snow,'signed snow gross')
 mt_ldew_rain=0;mt_ldew_snow=20
 mt_canopy_solid=0;mt_surface_solid=0;mt_waterstorage_solid=0;a_mt_precip=0
 expected_ratio=(20d0+7d0*input_snow)/(old_snow+input_snow)
 ! Rain passes through and triggers interception while old snow unloads.
 call mt_precip(1,1d0,1d0,input_snow,input_snow-drip,0d0,input_snow-drip, &
  1d0,drip,0d0,old_snow+input_snow-drip,0d0,old_snow, &
  0d0,gr,gs,0d0,0d0,0d0,0d0)
 call expect(mt_pg_rain_ground(1,1),5d0,'rain signature')
 call expect(mt_pg_snow_ground(1,1),drip*expected_ratio,'unloaded snow signature')
 call expect(mt_ldew_snow(1,1),(old_snow+input_snow-drip)*expected_ratio,'remaining snow')
 call expect(mt_ldew_snow(1,1)+mt_pg_snow_ground(1,1),20d0+7d0*input_snow,'snow budget')
enddo
call handoff(0.3d0,0.2d0,0.1d0,0.1d0,gr,gs)
call expect(gr,0.4d0,'positive net rain')
call expect(gs,0.3d0,'positive net snow')
print *, 'PASS canopy negative/zero/positive net interception and old-snow unloading'
contains
subroutine expect(actual,want,label)
real(r8) :: actual,want
character(*) :: label
if(abs(actual-want)>1d-12)then
 print *, 'FAIL ',label,actual,want
 error stop 1
endif
end subroutine
end program
'''
with tempfile.TemporaryDirectory(prefix='colm-canopy-') as d:
    p = Path(d)
    (p / 'probe.f90').write_text(source)
    subprocess.run(['gfortran', '-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                    '-ffree-line-length-none', 'probe.f90', '-o', 'probe'], cwd=p, check=True)
    subprocess.run(['./probe'], cwd=p, check=True, timeout=20)
