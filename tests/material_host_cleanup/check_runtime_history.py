#!/usr/bin/env python3
"""Exercise material history state and compare moved restart-count logic to frozen host."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def routine(source, name):
    return re.search(r'\bSUBROUTINE\s+' + name + r'\b.*?END SUBROUTINE\s+' + name + r'\b',
                     source, re.I | re.S).group()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True, help='frozen pre-phase-4 source tree')
    args = parser.parse_args()
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    baseline = args.baseline.resolve()
    old_count = routine((baseline / 'main/MOD_HistAccRestart.F90').read_text(), 'hacc_frcsat_count')
    old_count = re.sub(r'\bhacc_frcsat_count\b', 'old_count', old_count)
    old_acc = routine((baseline / 'main/MOD_Vars_1DAccFluxes.F90').read_text(), 'acc1d')
    registry = (ROOT / 'main/MATERIAL/MOD_MT_HistoryRestart.F90').read_text()
    fields = re.findall(r"CALL field([123])d\s*\(file, mode, '[^']+',\s*(\w+)", registry)
    declarations = '\n'.join(f'real(r8),allocatable::{name}({",".join(":" for _ in range(int(rank)))})'
                             for rank, name in fields)
    stubs = '''module MOD_Vars_Global
use MOD_Precision
real(r8),parameter::spval=-9999._r8
end module
module MOD_MT_Namelist
logical::enabled=.true.
contains
logical function mt_any_enabled()
mt_any_enabled=enabled
end function
end module
module MOD_SPMD_Task
logical::p_is_io=.true.,p_is_master=.true.
contains
subroutine CoLM_stop(message)
character(*),optional::message
if(present(message)) print '(A)',trim(message)
stop 7
end subroutine
end module
module MOD_Pixelset
type pixelset_type
integer::nblkgrp=1
integer::xblkgrp(1)=1,yblkgrp(1)=1
end type
end module
module MOD_LandPatch
use MOD_Pixelset
type(pixelset_type)::landpatch
end module
module MOD_Block
contains
subroutine get_filename_block(file,x,y,out)
character(*),intent(in)::file
integer,intent(in)::x,y
character(*),intent(out)::out
out=file
end subroutine
end module
module netcdf
integer,parameter::NF90_NOWRITE=0,NF90_NOERR=0,NF90_ENOTVAR=-49
logical::version_present=.true.,count_present=.true.,open_bad=.false.,close_bad=.false.,inquire_bad=.false.
contains
integer function nf90_open(file,mode,ncid)
character(*),intent(in)::file
integer,intent(in)::mode
integer,intent(out)::ncid
ncid=1
nf90_open=0
if(open_bad) nf90_open=-1
end function
integer function nf90_inq_varid(ncid,name,varid)
integer,intent(in)::ncid
character(*),intent(in)::name
integer,intent(out)::varid
varid=1
nf90_inq_varid=0
if(name=='hacc_frcsat_count_version'.and..not.version_present) nf90_inq_varid=NF90_ENOTVAR
if(name=='nac_frcsat'.and..not.count_present) nf90_inq_varid=NF90_ENOTVAR
if(inquire_bad) nf90_inq_varid=-1
end function
integer function nf90_close(ncid)
integer,intent(in)::ncid
nf90_close=0
if(close_bad) nf90_close=-1
end function
end module
module MOD_NetCDFVector
use MOD_Pixelset
contains
logical function ncio_vector_exist(file,name,pixelset)
character(*),intent(in)::file,name
type(pixelset_type),intent(in)::pixelset
ncio_vector_exist=.true.
end function
end module
module MOD_MT_Vars
use MOD_Precision
logical::mt_land_restart_loaded=.true.
''' + declarations + '\nend module\n'
    fixture = '''module fixture
use MOD_Precision
use MOD_Vars_Global
use MOD_MT_HostHistory
use MOD_MT_HistoryRestart
use MOD_MT_Namelist
use MOD_SPMD_Task
use MOD_LandPatch
use netcdf
use,intrinsic::ieee_arithmetic
implicit none
integer,parameter::HACC_WRITE=1
real(r8)::nac=3._r8,disk_count(4),version_value=1._r8
real(r8),allocatable::a_frcsat(:)
contains
subroutine require(ok,message)
logical,intent(in)::ok
character(*),intent(in)::message
if(.not.ok) then
print *,message
stop 8
endif
end subroutine
subroutine write_marker(file,name,val)
character(*),intent(in)::file,name
real(r8),intent(in)::val
print *, 'write_marker ',trim(name),val
end subroutine
subroutine read_marker(file,name,vmin,vmax)
character(*),intent(in)::file,name
real(r8),intent(out)::vmin,vmax
print *, 'read_marker ',trim(name)
vmin=version_value;vmax=version_value
end subroutine
subroutine require_var(file,name,pixelset,vecname,rank,n1,n2)
character(*),intent(in)::file,name,vecname
type(pixelset_type),intent(in)::pixelset
integer,intent(in)::rank,n1,n2
print *, 'require_var ',trim(name),rank,n1,n2
end subroutine
subroutine global_max(n)
integer,intent(inout)::n(:)
print *, 'global_max ',n
end subroutine
subroutine hacc_field1d(file,mode,name,var,pixelset,vecname)
character(*),intent(in)::file,name,vecname
integer,intent(in)::mode
real(r8),allocatable,intent(inout)::var(:)
type(pixelset_type),intent(in)::pixelset
print *, 'field1d ',trim(name),mode
if(mode==2.and.allocated(var)) var=disk_count
end subroutine
''' + old_acc + '\n' + old_count + '''
subroutine state_check()
real(r8)::values(4)
enabled=.false.
call mt_allocate_frcsat_valid(4)
call require(all(frcsat_hist_valid),'valid allocation')
call mt_allocate_frcsat_count(4)
call require(.not.allocated(nac_frcsat),'disabled count allocation')
enabled=.true.
call mt_allocate_frcsat_count(4)
call mt_flush_frcsat_count()
call require(all(nac_frcsat==0._r8),'count reset')
allocate(a_frcsat(4));a_frcsat=[spval,spval,spval,5._r8]
values=[1._r8,spval,2._r8,3._r8]
frcsat_hist_valid=[.true.,.true.,.false.,.true.]
call mt_accumulate_frcsat(values,a_frcsat,acc1d)
call require(all(a_frcsat==[1._r8,spval,spval,8._r8]),'filtered sums')
call require(all(nac_frcsat==[1._r8,0._r8,0._r8,1._r8]),'filtered counts')
nac_frcsat(1)=-1._r8
call mt_accumulate_frcsat(values,a_frcsat,acc1d)
call require(all(a_frcsat==[2._r8,spval,spval,11._r8]),'repeat sums')
call require(all(nac_frcsat==[-1._r8,0._r8,0._r8,2._r8]),'legacy count preserved')
enabled=.false.
call mt_accumulate_frcsat(values,a_frcsat,acc1d)
call require(all(a_frcsat==[3._r8,spval,2._r8,14._r8]),'disabled unfiltered sums')
call require(all(nac_frcsat==[-1._r8,0._r8,0._r8,2._r8]),'disabled count untouched')
call mt_deallocate_frcsat_count()
call mt_deallocate_frcsat_valid()
deallocate(a_frcsat)
call mt_flush_frcsat_count()
call mt_deallocate_frcsat_count()
call require(.not.allocated(nac_frcsat).and..not.allocated(frcsat_hist_valid),'release')
print *, 'PASS state allocation/reset/filter/legacy/disabled/release'
end subroutine
end module
program check
use fixture
implicit none
character(32)::mode,arg
integer::scenario,read_write
call get_command_argument(1,mode)
if(trim(mode)=='state')then
call state_check()
stop
endif
call get_command_argument(2,arg);read(arg,*)scenario
call mt_allocate_frcsat_count(4);call mt_flush_frcsat_count()
allocate(a_frcsat(4));a_frcsat=[1._r8,spval,2._r8,3._r8]
disk_count=[1._r8,0._r8,2._r8,3._r8]
read_write=2
select case(scenario)
case(1)
read_write=1
case(2)
case(3)
version_present=.false.;count_present=.false.
case(4)
version_present=.false.
case(5)
version_value=2._r8
case(6)
open_bad=.true.
case(7)
close_bad=.true.
case(8)
inquire_bad=.true.
case(9)
disk_count(1)=ieee_value(0._r8,ieee_quiet_nan)
case(10)
disk_count(1)=4._r8
case(11)
disk_count(1)=0.5_r8
case(12)
disk_count(2)=1._r8
case(13)
disk_count(1)=-2._r8
case(14)
enabled=.false.
case(15)
deallocate(nac_frcsat);deallocate(a_frcsat)
p_is_io=.false.;p_is_master=.false.
case(16)
disk_count(1)=-1._r8
end select
if(trim(mode)=='current')then
call mt_history_frcsat_count('fixture.nc',read_write,HACC_WRITE,nac,a_frcsat, &
hacc_field1d,write_marker,read_marker,require_var,global_max)
else
call old_count('fixture.nc',read_write)
endif
if(allocated(nac_frcsat)) print *, 'counts ',nac_frcsat
end program
'''
    with tempfile.TemporaryDirectory(prefix='colm-runtime-history-') as tmp:
        work = Path(tmp)
        (work / 'define.h').write_text('#define MaterialTransport\n')
        (work / 'stubs.F90').write_text(stubs)
        (work / 'fixture.F90').write_text(fixture)
        subprocess.run([compiler, '-O0', '-g', '-cpp', '-ffree-line-length-none', '-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow', '-I'+tmp,
                        str(ROOT/'share/MOD_Precision.F90'), 'stubs.F90',
                        str(ROOT/'main/MATERIAL/MOD_MT_HostHistory.F90'),
                        str(ROOT/'main/MATERIAL/MOD_MT_HistoryRestart.F90'),
                        'fixture.F90', '-o', 'check'], cwd=work, check=True)
        subprocess.run([str(work/'check'), 'state'], cwd=work, check=True)
        expected_errors = set(range(4, 14))
        for case in range(1, 17):
            outputs = [subprocess.run([str(work/'check'), mode, str(case)], cwd=work,
                                     capture_output=True, timeout=10) for mode in ['current', 'baseline']]
            assert outputs[0].returncode == outputs[1].returncode, (case, outputs)
            assert outputs[0].stdout == outputs[1].stdout, (case, outputs)
            assert outputs[0].returncode == (7 if case in expected_errors else 0), (case, outputs)
        print('PASS 16 frozen history-count cases: callback order, legacy metadata, errors and empty ranks')


if __name__ == '__main__':
    main()
