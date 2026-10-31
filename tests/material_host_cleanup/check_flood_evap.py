#!/usr/bin/env python3
"""Check production flood evaporation bookkeeping, optionally against saved CoLMMAIN."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=Path, help='pre-refactor CoLMMAIN.F90')
args = parser.parse_args()
source = (ROOT / 'main/MATERIAL/MOD_MT_RiverLakeCoupling.F90').read_text()
match = re.search(r'^   SUBROUTINE mt_rlc_land_flood_evap\b.*?'
                  r'^   END SUBROUTINE mt_rlc_land_flood_evap', source, re.M | re.S)
assert match, 'Missing production flood-evaporation helper'
reference = ''
if args.baseline:
    original = args.baseline.read_text()
    start = original.index('         IF (allocated(flood_mt_evap_patch)) THEN')
    block = original[start:original.index('\n#endif', start)]
    reference = '''subroutine reference_evap(ipatch,patchtype,flood_credit_patch,fevpg_fld, &
    deltim,temperature,forc_q,forc_psrf,forc_us,forc_vs)
use MOD_Namelist, only: LWINFILT => DEF_GridRiverLake_FloodFeedback
use MOD_MT_Evapo, only: mt_flood_evap_loss
integer, intent(in) :: ipatch,patchtype
real(r8), allocatable, intent(in) :: flood_credit_patch(:)
real(r8), intent(in) :: fevpg_fld,deltim,temperature,forc_q,forc_psrf,forc_us,forc_vs
type work_type
 real(r8) :: flood_evap_temperature
end type
type(work_type) :: mt_work
mt_work%flood_evap_temperature=temperature
''' + block + '\nend subroutine\n'
fixture = '''module MOD_Namelist
implicit none
logical :: DEF_GridRiverLake_FloodFeedback
end module
module MOD_MT_Evapo
use iso_fortran_env, only: real64
implicit none
integer, parameter :: r8=real64
integer :: calls=0, recorded_patch=0
real(r8) :: recorded(10)=0._r8, before(3)=0._r8
contains
subroutine mt_flood_evap_loss(ipatch,water_credit,material_credit,water_evap, &
    temperature,forc_q,forc_psrf,forc_us,forc_vs,material_loss)
integer, intent(in) :: ipatch
real(r8), intent(in) :: water_credit,material_credit(:),water_evap
real(r8), intent(in) :: temperature,forc_q,forc_psrf,forc_us,forc_vs
! Instrumented inout allows observing host clearing before the real out service.
real(r8), intent(inout) :: material_loss(:)
calls=calls+1
recorded_patch=ipatch
recorded=[water_credit,water_evap,temperature,forc_q,forc_psrf,forc_us,forc_vs,material_credit]
before=material_loss
material_loss=material_credit*water_evap/water_credit
end subroutine
end module
module fixture
use MOD_Namelist
use MOD_MT_Evapo, only: r8
implicit none
real(r8), allocatable :: flood_mt_credit_patch(:,:),flood_mt_evap_patch(:,:)
contains
'''
program = '''program check
use fixture
use MOD_MT_Evapo, only: calls, recorded_patch, recorded, before
use iso_fortran_env, only: int64
implicit none
real(r8), allocatable :: water(:)
real(r8) :: evaporation, actual(3,2), actual_recorded(10), actual_before(3)
integer :: scenario,patchtype,implementation,actual_calls,actual_patch,nimpl
nimpl=1
#ifdef HAVE_BASELINE
nimpl=2
#endif
do scenario=1,7
 do implementation=1,nimpl
  if(allocated(water))deallocate(water)
  if(allocated(flood_mt_credit_patch))deallocate(flood_mt_credit_patch)
  if(allocated(flood_mt_evap_patch))deallocate(flood_mt_evap_patch)
  allocate(water(2),flood_mt_credit_patch(3,2),flood_mt_evap_patch(3,2))
  water=[0.7_r8,0.125_r8]
  flood_mt_credit_patch(:,1)=999._r8
  flood_mt_credit_patch(:,2)=[125._r8,-250._r8,375._r8]
  flood_mt_evap_patch=-999._r8
  DEF_GridRiverLake_FloodFeedback=.true.
  patchtype=0
  evaporation=0.5_r8
  calls=0
  recorded_patch=0
  recorded=-333._r8
  before=-444._r8
  select case(scenario)
  case(2)
   DEF_GridRiverLake_FloodFeedback=.false.
   deallocate(water,flood_mt_credit_patch)
  case(3)
   patchtype=2
   deallocate(water,flood_mt_credit_patch)
  case(4)
   evaporation=0._r8
   deallocate(water,flood_mt_credit_patch)
  case(5)
   deallocate(water,flood_mt_credit_patch,flood_mt_evap_patch)
  case(6)
   evaporation=-0.5_r8
   deallocate(water,flood_mt_credit_patch)
  case(7)
   flood_mt_credit_patch(:,2)=0._r8
  end select
  if(implementation==1)then
   call mt_rlc_land_flood_evap(2,patchtype,water,evaporation,10._r8, &
      285._r8,0.004_r8,90000._r8,3._r8,-4._r8)
#ifdef HAVE_BASELINE
  else
   call reference_evap(2,patchtype,water,evaporation,10._r8, &
      285._r8,0.004_r8,90000._r8,3._r8,-4._r8)
#endif
  endif
  if(scenario==1.or.scenario==7)then
   if(calls/=1.or.recorded_patch/=2)error stop 'wrong callback count or patch'
   if(any(before/=0._r8))error stop 'ledger not cleared before callback'
   if(any(recorded(:7)/=[125._r8,5._r8,285._r8,0.004_r8,90000._r8,3._r8,-4._r8])) &
      error stop 'wrong conversion or atmospheric inputs'
   if(any(recorded(8:)/=flood_mt_credit_patch(:,2)))error stop 'wrong material credit'
   if(scenario==1)then
    if(any(flood_mt_evap_patch(:,2)/=[5._r8,-10._r8,15._r8]))error stop 'wrong output ledger'
   else
    if(any(flood_mt_evap_patch(:,2)/=0._r8))error stop 'zero credit ledger'
   endif
  else
   if(calls/=0)error stop 'inactive callback'
   if(allocated(flood_mt_evap_patch))then
    if(any(flood_mt_evap_patch(:,2)/=0._r8))error stop 'inactive ledger not cleared'
   endif
  endif
  if(allocated(flood_mt_evap_patch))then
   if(any(flood_mt_evap_patch(:,1)/=-999._r8))error stop 'other patch changed'
  endif
  if(implementation==1)then
   actual_calls=calls
   actual_patch=recorded_patch
   actual_recorded=recorded
   actual_before=before
   if(allocated(flood_mt_evap_patch))actual=flood_mt_evap_patch
  else
   if(actual_calls/=calls.or.actual_patch/=recorded_patch)error stop 'baseline callback mismatch'
   if(any(transfer(actual_recorded,[0_int64],10)/=transfer(recorded,[0_int64],10))) &
     error stop 'baseline input bits differ'
   if(any(transfer(actual_before,[0_int64],3)/=transfer(before,[0_int64],3))) &
     error stop 'baseline callback order differs'
   if(allocated(flood_mt_evap_patch))then
    if(any(transfer(actual,[0_int64],6)/=transfer(flood_mt_evap_patch,[0_int64],6))) &
      error stop 'baseline ledger bits differ'
   endif
  endif
 enddo
enddo
print *, 'PASS 7 flood-evaporation scenarios: guards, units, callback order and patch isolation'
end program
'''
compiler = shutil.which('gfortran')
assert compiler, 'gfortran is required'
with tempfile.TemporaryDirectory(prefix='colm-flood-evap-') as directory:
    tmp = Path(directory)
    (tmp / 'fixture.f90').write_text(fixture + match.group() + '\n' + reference + 'end module\n' + program)
    binary = tmp / 'check'
    subprocess.run([compiler, '-cpp', *(['-DHAVE_BASELINE'] if args.baseline else []),
                    '-ffree-line-length-none', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-O0',
                    str(tmp / 'fixture.f90'), '-o', str(binary)], cwd=tmp, check=True)
    subprocess.run([str(binary)], check=True)
    if args.baseline:
        print('PASS 7 bitwise baseline comparisons of callback inputs and ledger outputs')
