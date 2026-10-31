#!/usr/bin/env python3
"""Exercise the actual extracted thin-snow carry routine, including family isolation."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT/'main/MATERIAL/GENERIC_WATER/MOD_MT_Snow.F90').read_text()
match = re.search(r'^\s*SUBROUTINE mt_snowmelt_carry\b.*?END SUBROUTINE mt_snowmelt_carry', source, re.M|re.S|re.I)
assert match, 'production mt_snowmelt_carry missing'
fixture = '''module fixture
use iso_fortran_env, only: real64
implicit none
integer, parameter :: r8=real64
integer :: nmts=2
real(r8), parameter :: mt_tiny=1.e-30_r8
real(r8) :: mt_scv(2,2),mt_sm_carry(2,2)
contains
logical function mt_uses_generic_water(i)
integer,intent(in)::i
mt_uses_generic_water=i==1
end function
''' + match.group() + '''
end module
program check_snow
use fixture
implicit none
integer :: k
real(r8) :: scv,before,expect_store,expect_carry
integer :: snl
! Complete melt, partial melt, unchanged, growing, layered and zero material cases.
do k=1,6
 mt_scv=8._r8; mt_sm_carry=11._r8
 before=4._r8; scv=4._r8; snl=0
 expect_store=8._r8; expect_carry=0._r8
 select case(k)
 case(1)
  scv=0._r8; expect_store=0._r8; expect_carry=8._r8
 case(2)
  scv=1._r8; expect_store=2._r8; expect_carry=6._r8
 case(4)
  scv=8._r8
 case(5)
  snl=-2; scv=0._r8
 case(6)
  nmts=0; expect_carry=11._r8
 end select
 call mt_snowmelt_carry(1,snl,scv,before)
 if(mt_scv(1,1)/=expect_store.or.mt_sm_carry(1,1)/=expect_carry) error stop 'thin-snow partition'
 if(k<=5.and.mt_scv(1,1)+mt_sm_carry(1,1)/=8._r8) error stop 'material conservation'
 if(mt_scv(2,1)/=8._r8.or.mt_sm_carry(2,1)/=11._r8) error stop 'provider-owned material changed'
 if(any(mt_scv(:,2)/=8._r8).or.any(mt_sm_carry(:,2)/=11._r8)) error stop 'other patch changed'
enddo
print *, 'PASS 6 snowmelt cases: partition, conservation, provider and patch isolation'
end program
'''
compiler = shutil.which('gfortran')
if not compiler:
    raise SystemExit('gfortran is required')
with tempfile.TemporaryDirectory(prefix='colm-snow-carry-') as tmp:
    path=Path(tmp)
    (path/'check.f90').write_text(fixture)
    subprocess.run([compiler,'-O0','-g','-ffree-line-length-none','-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow','check.f90','-o','check'],cwd=path,check=True)
    subprocess.run([str(path/'check')],cwd=path,check=True,timeout=10)
