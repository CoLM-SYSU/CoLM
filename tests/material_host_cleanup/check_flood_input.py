#!/usr/bin/env python3
"""Check flood-input masses; optionally compare a saved pre-refactor CoLMMAIN."""
from pathlib import Path
import argparse
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=Path, help='pre-refactor CoLMMAIN.F90 for bitwise comparison')
args = parser.parse_args()
source = (ROOT / 'main/MATERIAL/MOD_MT_RiverLakeCoupling.F90').read_text()
match = re.search(r'^   SUBROUTINE mt_rlc_land_flood_input\b.*?'
                  r'^   END SUBROUTINE mt_rlc_land_flood_input', source, re.M | re.S)
assert match, 'Missing production flood-input helper'
block = ''
if args.baseline:
    original = args.baseline.read_text()
    start = original.index('         flood_input_material = 0._r8')
    block = original[start:original.index('\n#endif', start)]

fixture = '''module MOD_Namelist
implicit none
logical :: DEF_GridRiverLake_FloodFeedback
end module
module fixture
use iso_fortran_env, only: real64, int64
use MOD_Namelist
implicit none
integer, parameter :: r8=real64, nmts=3
type material_def
 real(r8) :: max_dissolved_conc
end type
type(material_def) :: materials(nmts)
real(r8), allocatable :: flood_mt_credit_patch(:,:), flood_mt_evap_patch(:,:)
contains
logical function mt_has_dissolved_limit(imt)
 integer, intent(in) :: imt
 mt_has_dissolved_limit=imt/=1
end function
subroutine CoLM_stop(message)
 character(len=*), intent(in) :: message
 print *, message
 error stop 17
end subroutine
'''
reference = '''subroutine reference_input(ipatch,patchtype,flood_credit_patch,fevpg_fld,qinfl_fld,deltim,flood_input_material)
use MOD_Namelist, only: LWINFILT => DEF_GridRiverLake_FloodFeedback
integer, intent(in) :: ipatch,patchtype
real(r8), intent(in) :: fevpg_fld,qinfl_fld,deltim
real(r8), intent(out) :: flood_input_material(nmts)
real(r8), allocatable, intent(in) :: flood_credit_patch(:)
integer :: imt_loc
'''+block+'''
end subroutine
''' if args.baseline else ''
program = '''program check
use fixture
implicit none
real(r8) :: actual(nmts), previous(nmts), evaporation, infiltration
real(r8), allocatable :: water_credit(:)
integer :: scenario,patchtype
character(len=32) :: mode
call get_command_argument(1,mode)
materials%max_dissolved_conc=[0._r8,0.25_r8,100._r8]
do scenario=1,9
 if(allocated(water_credit)) deallocate(water_credit)
 allocate(water_credit(2))
 if(allocated(flood_mt_credit_patch)) deallocate(flood_mt_credit_patch)
 if(allocated(flood_mt_evap_patch)) deallocate(flood_mt_evap_patch)
 allocate(flood_mt_credit_patch(nmts,2),flood_mt_evap_patch(nmts,2))
 flood_mt_credit_patch(:,1)=999._r8
 flood_mt_credit_patch(:,2)=[190._r8,95._r8,285._r8]
 flood_mt_evap_patch=0._r8
 DEF_GridRiverLake_FloodFeedback=.true.
 patchtype=0
 water_credit=0.1_r8
 evaporation=0.5_r8
 infiltration=2._r8
 if(len_trim(mode)>0)then
   ! More infiltration than the tolerated roundoff when no water is left.
   water_credit=0._r8
   evaporation=0._r8
   infiltration=2.e-13_r8
#ifdef HAVE_BASELINE
   if(trim(mode)=='reference-fatal')then
     call reference_input(2,patchtype,water_credit,evaporation,infiltration,10._r8,actual)
   else
#endif
     call mt_rlc_land_flood_input(2,patchtype,water_credit,evaporation,infiltration,10._r8,actual)
#ifdef HAVE_BASELINE
   endif
#endif
   error stop 'expected missing-credit rejection'
 endif
 select case(scenario)
 case(2)
   DEF_GridRiverLake_FloodFeedback=.false.
   deallocate(water_credit)
 case(3)
   patchtype=2
   deallocate(water_credit)
 case(4)
   infiltration=0._r8
 case(5)
   deallocate(flood_mt_credit_patch,flood_mt_evap_patch,water_credit)
 case(6)
   flood_mt_credit_patch(:,2)=[-19._r8,-9.5_r8,285._r8]
 case(7)
   water_credit=0._r8
   evaporation=0._r8
   infiltration=0.5e-13_r8
 case(8)
   evaporation=10._r8
   infiltration=0._r8
 case(9)
   flood_mt_evap_patch(:,2)=[9.5_r8,4.75_r8,0._r8]
 end select
 actual=-999._r8
 previous=-999._r8
 call mt_rlc_land_flood_input(2,patchtype,water_credit,evaporation,infiltration,10._r8,actual)
#ifdef HAVE_BASELINE
 call reference_input(2,patchtype,water_credit,evaporation,infiltration,10._r8,previous)
 if(any(transfer(actual,[0_int64],nmts)/=transfer(previous,[0_int64],nmts))) &
   error stop 'refactor changed flood input bits'
#endif
 select case(scenario)
 case(1)
   if(any(abs(actual-[40._r8,5._r8,60._r8])>1.e-12_r8))error stop 'fraction or solubility cap'
 case(6)
   if(any(abs(actual-[-4._r8,0._r8,60._r8])>1.e-12_r8))error stop 'family-specific negative clamp'
 case(9)
   if(any(abs(actual-[38._r8,5._r8,60._r8])>1.e-12_r8))error stop 'evaporation credit deduction'
 case default
   if(any(actual/=0._r8))error stop 'inactive or roundoff branch did not clear output'
 end select
enddo
print *, 'PASS 9 flood-input scenarios: expected mass and inactive branches'
end program
'''
compiler = shutil.which('gfortran')
assert compiler, 'gfortran is required'
with tempfile.TemporaryDirectory(prefix='colm-flood-input-') as directory:
    tmp = Path(directory)
    (tmp / 'fixture.f90').write_text(fixture + match.group() + '\n' + reference + 'end module\n' + program)
    binary = tmp / 'check'
    subprocess.run([compiler, '-cpp', *(['-DHAVE_BASELINE'] if args.baseline else []),
                    '-ffree-line-length-none', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-O0',
                    str(tmp / 'fixture.f90'), '-o', str(binary)], cwd=tmp, check=True)
    subprocess.run([str(binary)], check=True)
    for mode in ['fatal'] + (['reference-fatal'] if args.baseline else []):
        result = subprocess.run([str(binary), mode], capture_output=True, text=True)
        assert result.returncode != 0, 'Missing-credit branch must reject impossible infiltration'
        assert 'infiltration without published material credit' in result.stdout + result.stderr
    print('PASS production rejects infiltration without remaining water credit')
    if args.baseline:
        print('PASS 9 bitwise baseline comparisons and baseline missing-credit rejection')
