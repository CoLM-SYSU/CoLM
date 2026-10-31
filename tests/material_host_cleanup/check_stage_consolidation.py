#!/usr/bin/env python3
"""Compare actual old/new land-stage entry calls with recording physics services."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def routine(source, name):
    return re.search(r'\bSUBROUTINE\s+' + name + r'\b.*?END SUBROUTINE\s+' + name + r'\b', source, re.I | re.S).group()


def call_span(source, name):
    match = re.search(r'\bCALL\s+' + name + r'\s*\(', source, re.I)
    end = match.end()
    depth = 1
    while depth:
        depth += (source[end] == '(') - (source[end] == ')')
        end += 1
    return match.start(), end


def preprocess(path, directory):
    return subprocess.check_output(['gfortran', '-E', '-P', '-cpp', '-I' + str(directory), str(path)], text=True)


def recording_stub(source, name):
    body = '\n'.join(line.split('!')[0] for line in routine(source, name).splitlines())
    header = body[:body.index(')\n') + 1]
    declarations = [line for line in re.split(r'\n\s*CONTAINS\b', body, flags=re.I)[0].splitlines() if re.search(r'intent\s*\(in\)', line, re.I)]
    output = [header, 'implicit none', *declarations, f'write(unit) "{name}"']
    for line in declarations:
        for field in re.split(r',\s*(?![^()]*\))', line.split('::')[1]):
            var = field.strip().split('(')[0]
            if 'optional' in line.lower():
                output += [f'write(unit) present({var})', f'if(present({var})) write(unit) {var}']
            else:
                output += [f'write(unit) {var}']
    if name == 'mt_evapo':
        output += ['mt_scv(:,ipatch)=mt_scv(:,ipatch)+0.125_r8']
    return '\n'.join(output + [f'end subroutine {name}'])


def fixture(root, directory):
    main = preprocess(root / 'main/CoLMMAIN.F90', directory)
    land = preprocess(root / 'main/MATERIAL/MOD_MT_LandPhase.F90', directory)
    start, end = call_span(main, 'mt_land_water_step')
    prefix = main.rfind('mt_work%irrigation =', 0, start)
    if prefix >= 0:
        start = prefix
        publication = main.find('CALL mt_rlc_land_publish_input', end)
        if publication >= 0:
            _, tail = call_span(main[publication:], 'mt_rlc_land_publish_input')
            end = publication + tail
    water = main[start:end]
    start, end = call_span(main, 'mt_land_thermal_step')
    if 'CALL mt_snowmelt_carry' in main[end:main.index('CALL WATER_2014', end)]:
        _, tail = call_span(main[end:], r'mt_work%save_water')
        end += tail
    thermal = main[start:end]
    physics = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90').read_text()
    evap = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Evapo.F90').read_text()
    precip = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Precip.F90').read_text()
    snow = routine((root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Snow.F90').read_text(), 'mt_snowmelt_carry')
    assert snow.count('         ! snow melt') == 1, 'snow carry observation anchor changed'
    snow = snow.replace('         ! snow melt', '      write(unit) "carry",ipatch,snl,scv,scv_before,mt_work%wliq_before,mt_work%wice_before\n         ! snow melt')
    helpers = '\n'.join(routine(land, name) for name in ('mt_land_thermal_step', 'mt_land_water_step', 'mt_land_water_fluxes'))
    helpers += '\n' + routine(precip, 'mt_canopy_parts') + '\n' + snow
    helpers += '\n' + '\n'.join(recording_stub(source, name) for source, name in ((physics, 'mt_wetland'), (physics, 'mt_soil_water'), (evap, 'mt_evapo')))
    scalars = 'deltim rsur rsub qseva qsdew qsubl qfros qseva_soil qsdew_soil qsubl_soil qfros_soil qseva_snow qsdew_snow qsubl_snow qfros_snow etr sm fsno wa wdsrf wetwat pg_rain pg_snow forc_us forc_vs forc_q forc_psrf tleaf lai rst rss qinfl_fld fevpg_fld qinfl qflx_irrig_drip qflx_irrig_flood qflx_irrig_paddy ldew ldew_rain ldew_snow scv'.split()
    init = '\n'.join(f'{name}={100+i}._r8' for i, name in enumerate(scalars))
    return '''#include <define.h>
module stage_fixture
use MOD_Precision
use MOD_MT_LandWork
use MOD_MT_HydroDiagnostics
use MOD_MT_ThermalDiagnostics
implicit none
integer,parameter::nl_soil=4
real(r8),parameter::mt_tiny=1.e-30_r8,tfrz=273.15_r8
integer::nmts,ipatch=2,snl,patchtype,lb,lbsn,unit=20,scenario
logical::DEF_USE_Dynamic_Wetland,DEF_SPLIT_SOILSNOW,DEF_VEG_SNOW
real(r8)::@SCALARS@
real(r8)::qlayer(0:nl_soil),wliq_soisno(-4:nl_soil),wice_soisno(-4:nl_soil), &
 t_soisno(-4:nl_soil),porsl(nl_soil),dz_soisno(-4:nl_soil)
real(r8),allocatable::flood_credit_patch(:),flood_input_material(:),mt_scv(:,:),mt_sm_carry(:,:)
type(mt_land_work)::mt_work
type(mt_water_diagnostics)::water_diag
contains
@HELPERS@
logical function mt_uses_generic_water(i)
integer,intent(in)::i
mt_uses_generic_water=mod(i,2)==1
end function
subroutine mt_rlc_land_flood_input(ipatch,patchtype,water_credit,evaporation,infiltration,deltim,material_input)
integer,intent(in)::ipatch,patchtype
real(r8),allocatable,intent(in)::water_credit(:)
real(r8),intent(in)::evaporation,infiltration,deltim
real(r8),intent(out)::material_input(nmts)
integer::i
write(unit) 'input',ipatch,patchtype,allocated(water_credit),evaporation,infiltration,deltim
if(allocated(water_credit)) write(unit) water_credit
material_input=[(real(i,r8)+infiltration*deltim-evaporation,i=1,nmts)]
end subroutine
subroutine mt_rlc_land_publish_input(ipatch,patchtype,dynamic_wetland,material_input)
integer,intent(in)::ipatch,patchtype
logical,intent(in)::dynamic_wetland
real(r8),intent(in)::material_input(:)
write(unit) 'publish',ipatch,patchtype,dynamic_wetland,material_input
end subroutine
subroutine exercise()
@THERMAL@
write(unit) 'snapshot',mt_work%wliq_before,mt_work%wice_before,mt_scv,mt_sm_carry
@WATER@
write(unit) 'work',mt_work%irrigation,mt_work%infiltration,mt_work%recharge
end subroutine
subroutine setup()
integer::j
lb=snl+1;lbsn=min(lb,0)
@INIT@
if(scenario==0) then
qflx_irrig_drip=1.e16_r8;qflx_irrig_flood=-1.e16_r8;qflx_irrig_paddy=1
endif
scv=real(scenario,r8);tleaf=270+scenario
DEF_VEG_SNOW=mod(scenario,2)==1
qlayer=[(real(j,r8),j=1,size(qlayer))]
wliq_soisno=[(200._r8+j,j=1,size(wliq_soisno))];wice_soisno=wliq_soisno+20
t_soisno=wliq_soisno+40;porsl=0.4_r8;dz_soisno=1
allocate(mt_work%wliq_before(lb:nl_soil),mt_work%wice_before(lb:nl_soil))
mt_work%wliq_before=400;mt_work%wice_before=500
mt_work%wa_before=9;mt_work%wdsrf_before=10;mt_work%wetwat_before=11
mt_work%canopy_rain_prethermal=7;mt_work%canopy_snow_prethermal=8
mt_work%irrigation_storage_ground=12;mt_work%snow_before=2
call mt_thermal_reset(mt_work%thermal,lb,nl_soil)
mt_work%thermal%canopy=mt_canopy_diagnostics(1._r8,2._r8,3._r8)
mt_work%thermal%soil_thaw=4;mt_work%thermal%soil_freeze=5
call water_diag%reset(nl_soil,lbsn)
water_diag%etroot_aquifer=3;water_diag%rsub_source_aquifer=-2
water_diag%etroot_surface=0.125_r8;water_diag%rsub_source_surface=0.25_r8
water_diag%permeable_soil=[.true.,.false.,.true.,.false.]
water_diag%etroot_actual=6;water_diag%rsub_source_layer=7
water_diag%wblc_ice_sink=8;water_diag%snow_qout_layer=9
allocate(flood_input_material(nmts),mt_scv(nmts,2),mt_sm_carry(nmts,2))
mt_scv=5;mt_sm_carry=7
if(scenario>0) then
allocate(flood_credit_patch(2));flood_credit_patch=[0.25_r8,0.75_r8]
endif
end subroutine
subroutine cleanup()
deallocate(mt_work%wliq_before,mt_work%wice_before,mt_work%thermal%soil_thaw,mt_work%thermal%soil_freeze)
deallocate(flood_input_material,mt_scv,mt_sm_carry)
if(allocated(flood_credit_patch)) deallocate(flood_credit_patch)
end subroutine
end module
program check
use stage_fixture
implicit none
integer::s,p,d,k,n
open(unit,file='trace.bin',access='stream',form='unformatted',status='replace')
do n=0,3,3
do s=0,-3,-1
do p=0,2,2
do d=0,1
do k=0,1
do scenario=0,3
nmts=n;snl=s;patchtype=p;DEF_USE_Dynamic_Wetland=d==1;DEF_SPLIT_SOILSNOW=k==1
call setup()
call exercise()
call cleanup()
enddo
enddo
enddo
enddo
enddo
enddo
close(unit)
end program
'''.replace('@SCALARS@', ','.join(scalars)).replace('@HELPERS@', helpers).replace('@THERMAL@', thermal).replace('@WATER@', water).replace('@INIT@', init)


def run(root, directory, macros, optimization):
    directory.mkdir()
    (directory/'define.h').write_text('\n'.join('#define '+m for m in macros)+'\n')
    (directory/'fixture.F90').write_text(fixture(root, directory))
    sources = [root/'share/MOD_Precision.F90'] + [root/'main/MATERIAL'/name for name in ('MOD_MT_HydroDiagnostics.F90','MOD_MT_ThermalDiagnostics.F90','MOD_MT_LandWork.F90')]
    subprocess.run(['gfortran','-cpp','-ffree-line-length-none','-fcheck=all',optimization,'-I'+str(directory),*map(str,sources),'fixture.F90','-o','check'],cwd=directory,check=True)
    subprocess.run(['./check'],cwd=directory,check=True,timeout=20)
    return (directory/'trace.bin').read_bytes()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,required=True)
    args=parser.parse_args()
    if not shutil.which('gfortran'):raise SystemExit('gfortran required')
    count=0
    with tempfile.TemporaryDirectory(prefix='colm-stage-') as tmp:
        for grid,crop in ((False,False),(True,False),(False,True),(True,True)):
            macros=['MaterialTransport']+(['GridRiverLakeFlow'] if grid else [])+(['CROP'] if crop else [])
            for optimization in ('-O0','-O2'):
                directory=Path(tmp)/f'{grid}-{crop}-{optimization}';directory.mkdir()
                old=run(args.baseline.resolve(),directory/'old',macros,optimization)
                new=run(ROOT,directory/'new',macros,optimization)
                assert old==new, f'land stage trace differs: {macros}, {optimization}'
                count+=len(new)
    print(f'PASS 8 configurations x 256 thermal/water stage pairs: {count} bytes identical')


if __name__=='__main__':main()
