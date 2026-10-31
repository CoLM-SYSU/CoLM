#!/usr/bin/env python3
"""Record real material driver/history helpers; optionally compare frozen host blocks.

--baseline names the complete pre-phase-4 working-tree snapshot, not Git HEAD.
The services record every argument and assert callback order; no model data needed.
"""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def routine(source, name):
    found = re.search(r'\bSUBROUTINE\s+' + name + r'\b.*?END SUBROUTINE\s+' + name + r'\b',
                      source, re.I | re.S)
    if not found:
        raise ValueError(f'missing routine {name}')
    return found.group()


def call(source, name, optional=False):
    lines = source.splitlines()
    found = [index for index, line in enumerate(lines)
             if re.match(r'\s*CALL\s+' + name + r'\s*\(', line, re.I)]
    if optional and not found:
        return ''
    if len(found) != 1:
        raise ValueError(f'expected one call to {name}, found {len(found)}')
    start = end = found[0]
    while lines[end].rstrip().endswith('&'):
        end += 1
    return '\n'.join(lines[start:end + 1])


def cpp_block(source, start):
    lines = source[start:].splitlines(keepends=True)
    depth = 0
    for index, line in enumerate(lines):
        if re.match(r'\s*#\s*(if|ifdef|ifndef)\b', line):
            depth += 1
        elif re.match(r'\s*#\s*endif\b', line):
            depth -= 1
            if depth == 0:
                return ''.join(lines[:index + 1])
    raise ValueError('unclosed CPP block')


def material_block(source, call_name):
    position = source.index(call(source, call_name).strip())
    start = source.rfind('#ifdef MaterialTransport', 0, position)
    if start < 0:
        raise ValueError(f'{call_name} is not in a material block')
    return cpp_block(source, start)


SERVICES = '''#include <define.h>
module records
use MOD_Precision
implicit none
integer :: events(16), nevents=0, expected_step
contains
subroutine record(code)
integer,intent(in)::code
nevents=nevents+1
if(nevents>size(events)) error stop 'too many callbacks'
events(nevents)=code
write(20) code
end subroutine
subroutine expect(codes)
integer,intent(in)::codes(:)
if(nevents/=size(codes)) error stop 'callback count changed'
if(any(events(:nevents)/=codes)) error stop 'callback order changed'
nevents=0
end subroutine
subroutine bgc_driver(ipatch,idate,deltim,latitude,longitude)
integer,intent(in)::ipatch,idate(3)
real(r8),intent(in)::deltim,latitude,longitude
call record(12)
write(20) ipatch,idate,deltim,latitude,longitude
end subroutine
end module
module MOD_MT_Lifecycle
use records
implicit none
contains
subroutine material_lifecycle_lake_step(step,ipatch,idate,deltim,isub,nsub)
integer,intent(in)::step,ipatch,idate(3),isub,nsub
real(r8),intent(in)::deltim
if(step/=expected_step) error stop 'lake optional step changed'
call record(11)
write(20) step,ipatch,idate,deltim,isub,nsub
end subroutine
subroutine material_lifecycle_wetland_decomp(ipatch,deltim)
integer,intent(in)::ipatch
real(r8),intent(in)::deltim
call record(13)
write(20) ipatch,deltim
end subroutine
subroutine material_lifecycle_soil_step(step,ipatch,idate,deltim)
integer,intent(in)::step,ipatch,idate(3)
real(r8),intent(in)::deltim
if(step/=expected_step) error stop 'soil optional step changed'
call record(14)
write(20) step,ipatch,idate,deltim
end subroutine
subroutine mt_lifecycle_route_write_history(filename,itime)
character(*),intent(in)::filename
integer,intent(in)::itime
call record(21)
write(20) len(filename),filename,itime
end subroutine
subroutine mt_lifecycle_route_flush_history()
call record(23)
end subroutine
end module
#ifdef BGC
module MOD_MT_NitrogenLeaching
use records
implicit none
contains
subroutine mt_bgc_n_leaching(ipatch,patchtype,deltim)
integer,intent(in)::ipatch,patchtype
real(r8),intent(in)::deltim
call record(15)
write(20) ipatch,patchtype,deltim
end subroutine
end module
#endif
#ifdef GridRiverLakeFlow
module MOD_MT_RiverLake
use records
implicit none
contains
subroutine mt_riverlake_write_history(filename,itime,acctime)
character(*),intent(in)::filename
integer,intent(in)::itime
real(r8),intent(in)::acctime(:)
call record(20)
write(20) len(filename),filename,itime,size(acctime),acctime
end subroutine
subroutine mt_riverlake_flush_acc()
call record(22)
end subroutine
end module
#endif
'''


def fixture(root, baseline):
    land = (root / 'main/MATERIAL/MOD_MT_LandPhase.F90').read_text()
    driver = (root / 'main/CoLMDRIVER.F90').read_text()
    history = (root / 'main/HYDRO/MOD_Grid_RiverLakeHist.F90').read_text()
    names = ['mt_resolve_step', 'mt_lake_step', 'mt_wetland_decomp', 'mt_soil_step']
    if not baseline:
        names.append('mt_land_after_bgc')
    helpers = '\n'.join(routine(land, name) for name in names)
    bgc_position = driver.index(call(driver, 'bgc_driver').strip())
    start = driver.rfind('#if (defined BGC)', 0, bgc_position)
    if start < 0:
        raise ValueError('post-BGC dispatch lost its BGC guard')
    bgc = cpp_block(driver, start)
    write_name = 'mt_riverlake_write_history' if baseline else 'mt_rlc_write_history'
    flush_name = 'mt_riverlake_flush_acc' if baseline else 'mt_rlc_flush_history'
    history_write = material_block(history, write_name)
    history_flush = material_block(history, flush_name)
    main = (root / 'main/CoLMMAIN.F90').read_text()
    if baseline:
        start = main.index('IF (.not. (patchtype == 2 .and. .not. DEF_USE_Dynamic_Wetland)) THEN')
        publication = main[start:main.index('ENDIF', start) + len('ENDIF')]
    else:
        publication = call(land, 'mt_rlc_land_publish_input')
    if not baseline:
        coupling = (root / 'main/MATERIAL/MOD_MT_RiverLakeCoupling.F90').read_text()
        helpers += '\n#ifdef GridRiverLakeFlow\n' + '\n'.join(
            routine(coupling, name) for name in ['mt_rlc_write_history', 'mt_rlc_flush_history',
                                                'mt_rlc_land_publish_input']) + '\n#endif\n'
    return '''#include <define.h>
module fixture
use records
use MOD_MT_Lifecycle
#ifdef BGC
use MOD_MT_NitrogenLeaching
#endif
#ifdef GridRiverLakeFlow
use MOD_MT_RiverLake
#endif
implicit none
#ifdef GridRiverLakeFlow
real(r8),allocatable::flood_mt_land_patch(:,:)
#endif
contains
''' + helpers + '''
subroutine driver_case(ptype,istep_in)
integer,intent(in)::ptype
integer,intent(in),optional::istep_in
integer::istep_local,i,k,steps_in_one_deltim,idate(3),patchtype(2)
real(r8)::deltim,deltim_phy,patchlatr(2),patchlonr(2),PI
i=2; k=3; steps_in_one_deltim=4; idate=[2005,152,12345]
patchtype=ptype; deltim=1800._r8; deltim_phy=450._r8
patchlatr=[0.3_r8,0.51_r8]; patchlonr=[1.3_r8,1.71_r8]; PI=acos(-1._r8)
expected_step=1
if(present(istep_in)) expected_step=istep_in
write(20) ptype,present(istep_in),expected_step
''' + call(driver, 'mt_resolve_step', optional=True) + '\n' + call(driver, 'mt_lake_step') + '\n' + bgc + '''
#ifdef BGC
if(ptype==0)then
 call expect([11,12,13,14,15])
else
 call expect([11,13,14,15])
endif
#else
call expect([11])
#endif
end subroutine
#ifdef GridRiverLakeFlow
subroutine history_case(p_is_worker,n)
logical,intent(in)::p_is_worker
integer,intent(in)::n
integer::itime_in_file_ucat,i
real(r8)::acctime_ucat(n)
character(len=32)::file_hist_ucat
file_hist_ucat='case_hist_unitcat_2005-06.nc'
itime_in_file_ucat=7
acctime_ucat=[(17.25_r8*i,i=1,n)]
write(20) p_is_worker,n
''' + history_write + '''
if(p_is_worker)then
''' + history_flush + '''
endif
if(p_is_worker)then
 call expect([20,21,22,23])
else
 call expect([20,21])
endif
end subroutine
subroutine publication_case(ptype,dynamic,active,ns)
integer,intent(in)::ptype,ns
logical,intent(in)::dynamic,active
integer::ipatch,patchtype,j
logical::DEF_USE_Dynamic_Wetland
real(r8)::flood_input_material(ns)
patchtype=ptype; DEF_USE_Dynamic_Wetland=dynamic; ipatch=2
flood_input_material=[(1.25_r8*j-0.5_r8,j=1,ns)]
if(active)then
 allocate(flood_mt_land_patch(ns,3))
 flood_mt_land_patch=-7.5_r8
endif
! Inactive paths must not index the unpublished ledger.
if(.not.active.or.(ptype==2.and..not.dynamic)) ipatch=99
write(20) ptype,dynamic,active,ns,ipatch
''' + publication + '''
if(allocated(flood_mt_land_patch).neqv.active) error stop 'ledger allocation changed'
write(20) allocated(flood_mt_land_patch)
if(active)then
 if(any(flood_mt_land_patch(:,1)/=-7.5_r8).or.any(flood_mt_land_patch(:,3)/=-7.5_r8)) &
  error stop 'other patch ledger changed'
 if(ptype==2.and..not.dynamic)then
  if(any(flood_mt_land_patch(:,2)/=-7.5_r8)) error stop 'static wetland ledger changed'
 else
  if(any(flood_mt_land_patch(:,2)/=flood_input_material)) error stop 'ledger input not published'
 endif
 write(20) flood_mt_land_patch
 deallocate(flood_mt_land_patch)
endif
end subroutine
#endif
end module
program check
use fixture
implicit none
integer::ptype,dynamic,active,ns
open(20,file='trace.bin',access='stream',form='unformatted',status='replace')
do ptype=0,3,3
 call driver_case(ptype)
 call driver_case(ptype,-3)
 call driver_case(ptype,0)
 call driver_case(ptype,7)
enddo
#ifdef GridRiverLakeFlow
call history_case(.false.,0)
call history_case(.true.,0)
call history_case(.false.,3)
call history_case(.true.,3)
do ptype=0,3
 do dynamic=0,1
  do active=0,1
   do ns=0,3,3
    call publication_case(ptype,dynamic==1,active==1,ns)
   enddo
  enddo
 enddo
enddo
#endif
close(20)
end program
'''


def run(root, directory, compiler, bgc, grid, baseline=False):
    directory.mkdir()
    (directory / 'define.h').write_text('#define MaterialTransport\n' +
                                      ('#define BGC\n' if bgc else '') +
                                      ('#define GridRiverLakeFlow\n' if grid else ''))
    (directory / 'services.F90').write_text(SERVICES)
    (directory / 'fixture.F90').write_text(fixture(root, baseline))
    subprocess.run([compiler, '-O0', '-g', '-cpp', '-ffree-line-length-none',
                    '-fdefault-real-8', '-fdefault-double-8', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-I' + str(directory),
                    str(ROOT / 'share/MOD_Precision.F90'), 'services.F90',
                    'fixture.F90', '-o', 'check'], cwd=directory, check=True)
    subprocess.run([str(directory / 'check')], cwd=directory, check=True, timeout=15)
    return (directory / 'trace.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='colm-batch-routing-') as tmp:
        for bgc in [False, True]:
            for grid in [False, True]:
                label = f'bgc{int(bgc)}-grid{int(grid)}'
                current = run(ROOT, Path(tmp) / label, compiler, bgc, grid)
                if args.baseline:
                    old = run(args.baseline.resolve(), Path(tmp) / (label + '-baseline'),
                              compiler, bgc, grid, baseline=True)
                    if current != old:
                        raise SystemExit(f'FAIL {label}: callback arguments/order differ from frozen baseline')
                print(f'PASS {label}: 8 driver cases' + (' + 4 history + 32 publication cases' if grid else '') +
                      (f'; {len(current)} baseline-identical trace bytes' if args.baseline else ''))


if __name__ == '__main__':
    main()
