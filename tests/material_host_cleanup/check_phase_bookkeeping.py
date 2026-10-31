#!/usr/bin/env python3
"""Exercise extracted production phase helpers with real work storage and recording services.

No model data or Git history required. --baseline additionally compiles the original
host blocks and compares every callback input bitwise with the refactored helpers.
"""
import argparse
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=Path)
args = parser.parse_args()


def routine(source, name):
    return re.search(r'\bSUBROUTINE\s+' + name + r'\b.*?END SUBROUTINE\s+' + name + r'\b', source, re.I | re.S).group()


def interface(source, name):
    text = routine(source, name)
    lines = [line.split('!')[0] for line in re.split(r'\n\s*CONTAINS\b', text, flags=re.I)[0].splitlines()]
    header = text[:text.index(')') + 1]
    declarations = [line for line in lines if re.search(r'intent\s*\(', line, re.I)]
    return header + '\n' + '\n'.join(declarations)


def stub(source, name, body):
    signature = interface(source, name)
    records = [f'write(20) "{name}"']
    for line in signature.splitlines():
        if not re.search(r'intent\s*\(in\)', line, re.I):
            continue
        fields = re.split(r',\s*(?![^()]*\))', line.split('::')[1])
        for field in fields:
            field = field.strip().split('(')[0]
            if 'optional' in line.lower():
                records += [f'write(20) present({field})', f'if(present({field})) write(20) {field}']
            else:
                records += [f'write(20) {field}']
    return signature + '\n' + '\n'.join(records) + '\n' + body + '\nend subroutine ' + name


land = (ROOT / 'main/MATERIAL/MOD_MT_LandPhase.F90').read_text()
hist = (ROOT / 'main/MATERIAL/MOD_MT_Hist.F90').read_text()
precip = (ROOT / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Precip.F90').read_text()
cons = (ROOT / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Conservation.F90').read_text()
# Exercise the real THERMAL entry initializer after preparing its work object.
thermal = (ROOT / 'main/MOD_Thermal.F90').read_text()
thermal_reset = re.search(r'^\s*CALL mt_thermal_reset\(thermal_diag, lb, nl_soil\)', thermal, re.M).group().strip()
thermal_reset = thermal_reset.replace('thermal_diag', 'active%thermal').replace('nl_soil', '3')
helpers = '\n'.join([routine(land, 'mt_land_interception'), routine(land, 'mt_land_water_fluxes'), routine(hist, 'mt_land_finish_step')])
services = '\n'.join([
    stub(precip, 'mt_canopy_lump', '''call require(stage==0, 'canopy first')
stage=1
rain_eff=ldew_rain+1._r8
snow_eff=ldew_snow+2._r8'''),
    stub(cons, 'mt_save_storage', '''call require(stage==1, 'storage after canopy')
stage=2
#ifdef CROP
call require(present(waterstorage).eqv.DEF_USE_IRRIGATION, 'irrigation presence')
#else
call require(.not.present(waterstorage), 'no irrigation without CROP')
#endif
if(present(waterstorage)) call near(waterstorage,13._r8,'irrigation storage')'''),
    stub(precip, 'mt_precip', '''call require(stage==2,'precip after snapshot')
stage=3
call near(gross_intr_rain,3._r8,'negative net plus actual rain drip')
call near(gross_intr_snow,0._r8,'negative snow drip clipped')
call near(xsc_rain_out,2._r8,'old rain release')
call near(xsc_snow_out,0._r8,'snow release clipped')
call near(ldew_rain,4._r8,'effective intercepted rain')
call near(ldew_snow,6._r8,'effective intercepted snow')
call near(ldew_rain_old,7._r8,'pre-interception rain')
call near(ldew_snow_old,8._r8,'pre-interception snow')
call near(ldew_smelt_mass,0._r8,'no interception melt')
call near(ldew_frzc_mass,0._r8,'no interception freeze')'''),
    stub(cons, 'mt_apply_decay', '''call require(stage==3,'decay first at finish')
call require(allocated(active%wliq_before),'scratch alive at decay')
stage=4'''),
    stub(cons, 'mt_balance_check', '''call require(stage==4,'balance after decay')
call require(allocated(active%wliq_before),'scratch alive at balance')
stage=5
xerr_material=123._r8
call near(water_err_in,0.25_r8,'water error')
call near(water_dS_in,expected_ds,'storage change budget')
call near(water_input_in,expected_input,'water input budget')
call near(water_output_in,expected_output,'water output budget')
call near(water_evap_in,expected_output-8._r8,'evaporation budget')
call near(water_rnof_in,8._r8,'runoff budget')
#ifdef GridRiverLakeFlow
call require(present(flood_heterogeneous_in),'flood heterogeneity present')
call require(flood_heterogeneous_in.eqv.(enabled.and.ptype==0),'heterogeneity value')
#else
call require(.not.present(flood_heterogeneous_in),'heterogeneity absent')
#endif'''),
    stub(hist, 'mt_hist_accumulate', '''call require(stage==5,'history after balance')
call require(allocated(active%wliq_before).and.allocated(active%snow_ice_before),'scratch alive at history')
call near(active%balance_error,123._r8,'balance result retained')
call require(size(wliq_soisno)==nl_soil-snl,'history snow extent')
call near(wliq_soisno(snl+1),21._r8+snl+1,'history lower layer')
stage=6'''),
])
fixture = '''#include <define.h>
module fixture
use MOD_Precision
use MOD_MT_HydroDiagnostics
use MOD_MT_LandWork
use MOD_MT_ThermalDiagnostics, only: mt_thermal_reset
implicit none
real(r8),parameter::mt_tiny=1.e-30_r8
type(mt_land_work),target::active
type(mt_water_diagnostics)::diag
logical::DEF_USE_IRRIGATION,enabled
integer::stage,ptype
real(r8)::expected_ds,expected_input,expected_output
contains
subroutine require(ok,label)
logical::ok
character(*)::label
if(.not.ok)then
 print *, 'FAIL ',label
 stop 1
endif
end subroutine
subroutine near(value,expected,label)
real(r8)::value,expected
character(*)::label
call require(abs(value-expected)<1.e-12_r8,label)
end subroutine
''' + services + '\n@HELPERS@\nend module\n' + '''program check
use fixture
implicit none
integer::s,lb,i,k,ir,split,p,f
real(r8)::liquid(-3:3),ice(-3:3),temp(-3:3)
open(20,file='trace.bin',access='stream',form='unformatted',status='replace')
do i=-3,3
 liquid(i)=20._r8+i
 ice(i)=40._r8+i
 temp(i)=270._r8+i
enddo
do s=0,-3,-3
do ir=0,1
do split=0,1
do p=0,2,2
do f=0,1
ptype=p;enabled=f==1;DEF_USE_IRRIGATION=ir==1
lb=s+1
! Repeat on the same object to detect retained allocation or patch state.
do k=1,2
call active%prepare_snow(-4,3,s,ice)
call require(lbound(active%snow_ice_before,1)==-3,'snow allocation lower bound')
if(s<0) call require(all(active%snow_ice_before(s+1:0)==ice(s+1:0)),'snow snapshot')
call require(all(active%snow_ice_before(1:)==0),'soil excluded from snow snapshot')
call require(all(active%snow_ice_before(-3:s)==0),'inactive snow zero')
active%thermal%canopy%melt=1;active%thermal%canopy%freeze=2
active%thermal%canopy%raw=3;active%recharge=4
call active%prepare_thermal(lb,3,liquid(lb:),ice(lb:),10._r8,11._r8,12._r8,diag,temp(lb:),0.25_r8,split==1)
call require(lbound(active%wliq_before,1)==lb,'thermal allocation lower bound')
call require(all(active%wliq_before==liquid(lb:)).and.all(active%wice_before==ice(lb:)),'thermal snapshot')
call near(active%thermal%canopy%melt,1._r8,'work preparation preserves thermal output')
@THERMAL_RESET@
call require(all(active%thermal%soil_thaw==0).and.all(active%thermal%soil_freeze==0),'thermal scratch reset')
call near(active%wa_before,10._r8,'aquifer snapshot')
call near(active%wdsrf_before,11._r8,'surface snapshot')
call near(active%wetwat_before,12._r8,'wetland snapshot')
call near(active%thermal%canopy%melt+active%thermal%canopy%freeze+active%thermal%canopy%raw+active%recharge,0._r8,'thermal scalar reset')
call near(diag%etroot_aquifer+diag%rsub_source_surface,0._r8,'diagnostics reset between patches')
#ifdef GridRiverLakeFlow
if(split==1)then
 call near(active%flood_evap_temperature,0.25_r8*temp(lb)+0.75_r8*temp(1),'split flood temperature')
else
 call near(active%flood_evap_temperature,temp(lb),'unsplit flood temperature')
endif
#endif
call active%save_water(lb,3,liquid(lb:)+1._r8,ice(lb:)+2._r8)
call require(all(active%wliq_before==liquid(lb:)+1).and.all(active%wice_before==ice(lb:)+2),'post-thermal snapshot')
active%irrigation_storage_before=13
active%canopy_rain_before=7;active%canopy_snow_before=8
stage=0
call mt_land_interception(2,2._r8,s,3,active,7._r8,3._r8,4._r8,280._r8,2._r8,3._r8, &
 -3._r8,-2._r8,-1._r8,7._r8,3._r8,1._r8,5._r8,-4._r8)
diag%etroot_aquifer=4;diag%rsub_source_aquifer=-5
diag%etroot_surface=2;diag%rsub_source_surface=4
call mt_land_water_fluxes(active,diag,12._r8,10._r8,2._r8)
call near(active%recharge,3._r8,'negative aquifer source excluded')
call near(active%infiltration,7._r8,'surface extraction excluded')
diag%rsub_source_aquifer=6
call mt_land_water_fluxes(active,diag,12._r8,10._r8,2._r8)
call near(active%recharge,6._r8,'positive aquifer source restored')
write(20) active%recharge,active%infiltration
expected_ds=10;expected_input=6;expected_output=18
#if defined CaMa_Flood || defined GridRiverLakeFlow
if(p==0)expected_output=16
if(p==0.and.enabled)then
 expected_ds=16;expected_input=12
endif
#endif
call mt_land_finish_step(active,2,s,-4,3,2._r8,p,0.25_r8,30._r8,20._r8,1._r8,2._r8,5._r8,4._r8, &
 3._r8,4._r8,liquid(s+1:)+1._r8,ice(s+1:)+2._r8,12._r8,11._r8,12._r8,9._r8 &
#if defined CaMa_Flood || defined GridRiverLakeFlow
 ,enabled,3._r8,1._r8 &
#endif
 )
call require(stage==6,'all finish callbacks executed')
call require(.not.allocated(active%wliq_before).and..not.allocated(active%wice_before), 'water snapshots released')
call require(.not.allocated(active%thermal%soil_thaw).and..not.allocated(active%thermal%soil_freeze).and..not.allocated(active%snow_ice_before),'scratch released')
enddo
enddo
enddo
enddo
enddo
enddo
close(20)
end program
'''


def baseline_helpers(path):
    old = (path / 'main/CoLMMAIN.F90').read_text()
    interception = old[old.index('         mt_work%rain_release ='):old.index('\n#endif\n\n!----', old.index('         mt_work%rain_release ='))]
    flux_start = old.index('         mt_work%recharge = (wa')
    flux = old[flux_start:old.index('\n#ifdef GridRiverLakeFlow', flux_start)]
    finish_start = old.index('         mt_work%flood_input = 0._r8')
    finish = old[finish_start:old.index('\n#endif\n\n#if', finish_start)]
    replacements = {'mt_work': 'work', 'mt_leaf_qdrp_rain': 'canopy_drain_rain', 'mt_leaf_qdrp_snow': 'canopy_drain_snow', 'LWINFILT': 'flood_enabled', 'qinfl_fld': 'flood_input', 'fevpg_fld': 'flood_evaporation'}
    blocks = []
    for source, name, body in [(land,'mt_land_interception',interception),(land,'mt_land_water_fluxes',flux),(hist,'mt_land_finish_step',finish)]:
        for before, after in replacements.items():
            body = re.sub(r'\b'+before+r'\b',after,body)
        # Preserve preprocessor-controlled declarations in the actual interface.
        current = routine(source,name)
        if name=='mt_land_finish_step':
            declaration = current[:current.index('      work%flood_input')]
        else:
            declaration = interface(source,name)
        blocks.append(declaration+'\n'+body+'\nend subroutine '+name)
    return '\n'.join(blocks)

compiler = shutil.which('gfortran')
if not compiler:
    raise SystemExit('gfortran required')
with tempfile.TemporaryDirectory(prefix='colm-phase-bookkeeping-') as tmp:
    path = Path(tmp)
    for crop in [False, True]:
        for flood in ['', 'GridRiverLakeFlow', 'CaMa_Flood']:
            (path/'define.h').write_text('#define MaterialTransport\n'+('#define CROP\n' if crop else '')+(f'#define {flood}\n' if flood else ''))
            traces = []
            for variant in [helpers] + ([baseline_helpers(args.baseline)] if args.baseline else []):
                (path/'fixture.F90').write_text(fixture.replace('@HELPERS@',variant).replace('@THERMAL_RESET@',thermal_reset))
                subprocess.run([compiler,'-O0','-g','-cpp','-ffree-line-length-none','-fcheck=all','-I'+tmp,
                    str(ROOT/'share/MOD_Precision.F90'),str(ROOT/'main/MATERIAL/MOD_MT_HydroDiagnostics.F90'),
                    str(ROOT/'main/MATERIAL/MOD_MT_ThermalDiagnostics.F90'),str(ROOT/'main/MATERIAL/MOD_MT_LandWork.F90'),'fixture.F90','-o','check'],cwd=path,check=True)
                subprocess.run([str(path/'check')],cwd=path,check=True,timeout=10)
                traces.append((path/'trace.bin').read_bytes())
            if len(traces)==2:
                assert traces[0]==traces[1], f'callback values/order differ: CROP={crop}, flood={flood}'
            print(f'PASS 64 bookkeeping/lifecycle cases: CROP={crop}, flood={flood or "off"}'+(', baseline trace bitwise equal' if args.baseline else '')+f', sha256={hashlib.sha256(traces[0]).hexdigest()}')
    # Compile the actual CaMa guard separately; no routing model/data needed.
    guard = '''module guard_fixture
implicit none
integer::nmts=2
logical::mt_uses_generic_water(2)=[.false.,.true.]
contains
subroutine CoLM_stop(message)
character(*)::message
print '(a)',message
stop 17
end subroutine
''' + routine(land, 'mt_land_check_flood_feedback') + '''
end module
program check_guard
use guard_fixture
implicit none
character(1)::case_id
call get_command_argument(1,case_id)
select case(case_id)
case('1')
 call mt_land_check_flood_feedback(0,.false.)
case('2')
 call mt_land_check_flood_feedback(2,.true.)
case('3')
 nmts=0
 call mt_land_check_flood_feedback(0,.true.)
case('4')
 mt_uses_generic_water=.false.
 call mt_land_check_flood_feedback(0,.true.)
case('5')
 call mt_land_check_flood_feedback(0,.true.)
end select
end program
'''
    (path/'guard.F90').write_text(guard)
    subprocess.run([compiler,'-O0','-g','-fcheck=all','guard.F90','-o','guard'],cwd=path,check=True)
    for case in range(1,6):
        result = subprocess.run([str(path/'guard'),str(case)],cwd=path,capture_output=True,text=True,timeout=10)
        if case < 5:
            assert result.returncode == 0 and not result.stdout and not result.stderr, (case,result)
        else:
            assert result.returncode == 17, result
            assert result.stdout.strip() == 'Material transport: CaMa flood feedback (LWINFILT) is not supported.', result.stdout
    print('PASS 5 CaMa material guard cases: disabled/nonsoil/empty/provider-only allowed; generic rejected with original message')
