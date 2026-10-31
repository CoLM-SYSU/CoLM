#!/usr/bin/env python3
"""Compare host boundary observations against a frozen pre-cleanup source tree.

Water solvers and history callbacks are recording stubs; thermal work/diagnostic
storage is real. No model data, MPI, or third-party Python packages are required.
"""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def call_span(text, name):
    match = re.search(r'\bCALL\s+' + name + r'\s*\(', text, re.I)
    if match is None:
        raise AssertionError(f'missing production call: {name}')
    depth = 1
    end = match.end()
    while depth:
        depth += (text[end] == '(') - (text[end] == ')')
        end += 1
    return match.start(), end


def preprocess(source, directory):
    path = directory / 'fragment.F90'
    path.write_text(source)
    return subprocess.check_output(['gfortran', '-E', '-P', '-cpp',
                                    '-I' + str(directory), str(path)], text=True)


def fixture(root, directory):
    main = (root / 'main/CoLMMAIN.F90').read_text()
    start = main.index('         IF (.not. DEF_USE_VariablySaturatedFlow) THEN',
                       main.index('CALL THERMAL'))
    flux_call = main.find('CALL mt_land_water_fluxes', start)
    end = flux_call if flux_call >= 0 else main.index('CALL mt_land_water_step', start)
    # Close the still-open MaterialTransport guard before preprocessing the slice.
    water = preprocess('#include <define.h>\n' + main[start:end] + '\n#endif\n', directory)
    water = '\n'.join(line.split('!')[0] for line in water.splitlines())
    if flux_call < 0:
        # Stage consolidation owns this observation; full call mapping/order is
        # also checked by check_stage_consolidation.py.
        land = (root / 'main/MATERIAL/MOD_MT_LandPhase.F90').read_text()
        assignment = re.search(r'^\s*work%irrigation\s*=.*$', land, re.M).group().strip()
        water += '\n' + assignment.replace('work%', 'mt_work%')
    for name, branch in [('WATER_2014', 1), ('WATER_VSF', 2)]:
        a, b = call_span(water, name)
        args = [re.search(r'\b' + name + r'\s*=\s*(\w+)', water[a:b], re.I).group(1)
                for name in ('qflx_irrig_drip', 'qflx_irrig_flood', 'qflx_irrig_paddy')]
        water = water[:a] + f'call change_water({branch},' + ','.join(args) + ')' + water[b:]

    thermal = (root / 'main/MOD_Thermal.F90').read_text()
    a, b = call_span(thermal, 'mt_thermal_reset')
    reset = thermal[a:b].replace('thermal_diag', 'mt_work%thermal')
    history = (root / 'main/MOD_HistAccRestart.F90').read_text()
    fields = history[:history.index('END SUBROUTINE hacc_fields')]
    legacy = re.search(r'\bCALL\s+hacc_frcsat_count\s*\(', fields, re.I)
    a, b = call_span(fields, 'hacc_frcsat_count' if legacy else 'mt_history_frcsat_count')
    history_call = fields[a:b]
    wrapper = ''
    if legacy:
        wrapper = re.search(r'\bSUBROUTINE hacc_frcsat_count\b.*?END SUBROUTINE hacc_frcsat_count',
                            history, re.I | re.S).group()
    callbacks = '\n'.join(f'integer function {name}()\n{name}={i}\nend function'
                          for i, name in enumerate(('hacc_field1d', 'write_marker',
                              'read_marker', 'require_var', 'global_max'), 11))
    return r'''#include <define.h>
program boundary_regression
use MOD_Precision
use MOD_MT_LandWork
use MOD_MT_HydroDiagnostics
use MOD_MT_ThermalDiagnostics
implicit none
type(mt_land_work)::mt_work
type(mt_water_diagnostics)::water_diag
integer,parameter::nl_soil=3,HACC_WRITE=1
integer::branch,scenario,k,lb,mode,nac,unit
logical::DEF_USE_VariablySaturatedFlow
real(r8)::qflx_irrig_drip,qflx_irrig_flood,qflx_irrig_paddy,rsub,rnof,rsur
real(r8)::liquid(-3:3),ice(-3:3),temperature(-3:3),a_frcsat(3)
character(len=12)::file='restart-test'
open(newunit=unit,file='observations.bin',access='stream',form='unformatted',status='replace')
do branch=1,2
DEF_USE_VariablySaturatedFlow=branch==2
do scenario=1,8
qflx_irrig_drip=real(scenario,r8);qflx_irrig_flood=2; qflx_irrig_paddy=3
mt_work%irrigation=-999;rnof=17;rsur=11;rsub=-1
''' + water + r'''
write(unit) branch,scenario,qflx_irrig_drip,qflx_irrig_flood,qflx_irrig_paddy,mt_work%irrigation,rsub
enddo
enddo
liquid=7;ice=9;temperature=271
do k=1,12
lb=1-mod(k,5)
call mt_work%prepare_snow(-4,nl_soil,min(lb-1,0),ice)
! Seed stale observations, including arrays with mismatched bounds.
call mt_thermal_reset(mt_work%thermal,-5,5)
mt_work%thermal%soil_thaw=99;mt_work%thermal%soil_freeze=88
mt_work%thermal%canopy=mt_canopy_diagnostics(77._r8,66._r8,55._r8)
call mt_work%prepare_thermal(lb,nl_soil,liquid(lb:),ice(lb:),10._r8,11._r8,12._r8, &
 water_diag,temperature(lb:),0.25_r8,mod(k,2)==0)
''' + reset + r'''
if(lbound(mt_work%thermal%soil_thaw,1)/=lb.or.ubound(mt_work%thermal%soil_thaw,1)/=nl_soil) error stop 'bounds'
if(any(mt_work%thermal%soil_thaw/=0).or.any(mt_work%thermal%soil_freeze/=0)) error stop 'stale phase'
if(mt_work%thermal%canopy%melt/=0.or.mt_work%thermal%canopy%freeze/=0.or.mt_work%thermal%canopy%raw/=0) error stop 'stale canopy'
write(unit) lb,mt_work%thermal%soil_thaw,mt_work%thermal%soil_freeze,mt_work%thermal%canopy
write(unit) mt_work%wliq_before,mt_work%wice_before,mt_work%wa_before,mt_work%wdsrf_before,mt_work%wetwat_before,mt_work%recharge
#ifdef GridRiverLakeFlow
write(unit) mt_work%flood_evap_temperature
#endif
call mt_work%release()
enddo
do mode=1,3
do nac=0,2
a_frcsat=[0._r8,0.25_r8,1._r8]
''' + history_call + r'''
enddo
enddo
close(unit)
contains
subroutine change_water(which,drip,flood,paddy)
integer,intent(in)::which
real(r8),intent(inout)::drip,flood,paddy
! Mutations detect premature sampling; cancellation detects changed sum ordering.
select case(mod(scenario,4))
case(0)
drip=1.e16_r8;flood=-1.e16_r8;paddy=real(which,r8)
case(1)
drip=0;flood=3._r8*which;paddy=7
case(2)
drip=real(which,r8)/3;flood=0;paddy=2._r8/7
case(3)
drip=-real(which,r8);flood=11;paddy=0
end select
end subroutine
subroutine mt_history_frcsat_count(path,operation,write_mode,count,values,f1,f2,f3,f4,f5)
character(len=*),intent(in)::path
integer,intent(in)::operation,write_mode,count
real(r8),intent(in)::values(:)
integer,external::f1,f2,f3,f4,f5
write(unit) path,operation,write_mode,count,values,f1(),f2(),f3(),f4(),f5()
end subroutine
''' + wrapper + '\n' + callbacks + '\nend program\n'


def run(root, work, macros, optimization):
    work.mkdir()
    (work / 'define.h').write_text('\n'.join('#define ' + macro for macro in macros) + '\n')
    (work / 'fixture.F90').write_text(fixture(root, work))
    sources = [root / 'share/MOD_Precision.F90']
    sources += [root / 'main/MATERIAL' / name for name in (
        'MOD_MT_HydroDiagnostics.F90', 'MOD_MT_ThermalDiagnostics.F90', 'MOD_MT_LandWork.F90')]
    command = ['gfortran', '-cpp', '-ffree-line-length-none', '-fcheck=all', optimization,
               '-I' + str(work), '-J' + str(work), *map(str, sources), 'fixture.F90', '-o', 'check']
    subprocess.run(command, cwd=work, check=True, capture_output=True, text=True)
    subprocess.run(['./check'], cwd=work, check=True, capture_output=True, text=True)
    return (work / 'observations.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    args = parser.parse_args()
    if shutil.which('gfortran') is None:
        raise SystemExit('gfortran is required')
    total = 0
    with tempfile.TemporaryDirectory(prefix='colm-optional-') as temporary:
        for index, (crop, grid) in enumerate(((False, False), (True, False), (False, True), (True, True))):
            macros = ['MaterialTransport'] + (['CROP'] if crop else []) + (['GridRiverLakeFlow'] if grid else [])
            for optimization in ('-O0', '-O2'):
                directory = Path(temporary) / f'{index}{optimization}'
                directory.mkdir()
                expected = run(args.baseline.resolve(), directory / 'baseline', macros, optimization)
                actual = run(ROOT, directory / 'current', macros, optimization)
                assert expected == actual, f'boundary observations differ: {macros}, {optimization}'
                total += len(actual)
    print(f'PASS: 8 configurations, 128 water branches, 96 thermal lifecycles, 72 history calls; {total} bytes identical')


if __name__ == '__main__':
    main()
