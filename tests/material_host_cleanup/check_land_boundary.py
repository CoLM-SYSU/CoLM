#!/usr/bin/env python3
"""Compare old/new land entry calls by recording every physics argument."""
from pathlib import Path
import argparse
import re
import shutil
import subprocess
import tempfile

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline-main',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[2]
old=args.baseline_main.read_text()
land=(root/'main/MATERIAL/MOD_MT_LandPhase.F90').read_text()
main=(root/'main/CoLMMAIN.F90').read_text()
physics=(root/'main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90').read_text()

def clean(text):
 return '\n'.join(line.split('!')[0] for line in text.splitlines())

def call(text,name):
 start=re.search(r'CALL\s+'+name+r'\s*\(',text,re.I).start()
 begin=text.index('(',start); depth=1; end=begin+1
 while depth:
  if text[end]=='(':depth+=1
  if text[end]==')':depth-=1
  end+=1
 return text[start:end]

def routine(text,name):
 return re.search(r'\bSUBROUTINE\s+'+name+r'\b.*?END SUBROUTINE\s+'+name+r'\b',text,re.I|re.S).group()

def split(text):
 result=[]; depth=0; start=0
 for i,c in enumerate(text):
  if c=='(':depth+=1
  if c==')':depth-=1
  if c==',' and depth==0:
   result.append(text[start:i].strip());start=i+1
 result.append(text[start:].strip())
 return result

def stub(name):
 text=clean(routine(physics,name))
 header=text[:text.index(')\n')+1]
 fields=[]; decl=[]
 for line in re.split(r'\n\s*CONTAINS\b',text,flags=re.I)[0].splitlines():
  if 'intent(in)' in line:
   decl.append(line)
   for f in split(line.split('::')[1]):
    fields.append((f.split('(')[0].strip(),'optional' in line))
 body=[header,'implicit none',*decl,'write(output_unit) "'+name+'"']
 for name,optional in fields:
  if optional: body += ['write(output_unit) present('+name+')','if(present('+name+')) then']
  body += ['write(output_unit) '+name]
  if optional: body += ['endif']
 return '\n'.join(body+['end subroutine'])

work_map={
'wliq_soisno_old_mt':'wliq_before','wice_soisno_old_mt':'wice_before',
'wa_old_mt':'wa_before','wdsrf_old_mt':'wdsrf_before','wetwat_old_mt':'wetwat_before',
'qinfl_mt':'infiltration','qcharge_mt':'recharge','irrig_ground_mt':'irrigation',
'waterstorage_mt_ground':'irrigation_storage_ground','raw_mt':'thermal%canopy%raw'}
diag_map={name+'_mt':name for name in ['wblc_ice_sink','etroot_actual','etroot_aquifer','etroot_surface',
'dew_overflow','frost_displaced','rsub_source_layer','rsub_source_surface','rsub_source_aquifer',
'imperv_evap_wdsrf','imperv_evap_soil','imperv_subl_soil','snow_qout_layer','permeable_soil']}
diag_map.update({'late_runoff_mt':'late_runoff','qgtop_mt':'qgtop'})
reference='if(patchtype==2.and..not.DEF_USE_Dynamic_Wetland)then\n'+call(clean(old),'mt_wetland')+'\nelse\n'+call(clean(old),'mt_soil_water')+'\nendif'
for source,target in work_map.items():reference=re.sub(r'\b'+source+r'\b','mt_work%'+target,reference)
for source,target in diag_map.items():reference=re.sub(r'\b'+source+r'\b','water_diag%'+target,reference)
actual=call(clean(main),'mt_land_water_step')
wrapper=routine(land,'mt_land_water_step')
scalars='deltim rsur rsub qseva qsdew qsubl qfros qseva_soil qsdew_soil qsubl_soil qfros_soil qseva_snow qsdew_snow qsubl_snow qfros_snow etr sm fsno wa wdsrf wetwat pg_rain pg_snow forc_us forc_vs forc_q forc_psrf tleaf lai rst rss qinfl_fld fevpg_fld qinfl qflx_irrig_drip qflx_irrig_flood qflx_irrig_paddy'.split()
arrays={'qlayer':'0:nl_soil','wliq_soisno':'-4:nl_soil','wice_soisno':'-4:nl_soil',
 't_soisno':'-4:nl_soil','porsl':'nl_soil','dz_soisno':'-4:nl_soil','flood_input_material':'2'}
initial=[]
for i,name in enumerate(scalars):initial.append(f'{name}={100+i}._r8')
for i,name in enumerate(arrays):initial.append(f'{name}=[({200+i*20}._r8+j,j=1,size({name}))]')
for i,name in enumerate(work_map.values()):
 if name in ['wliq_before','wice_before']:
  initial += [f'if(allocated(mt_work%{name}))deallocate(mt_work%{name})',f'allocate(mt_work%{name}(snl+1:nl_soil))',f'mt_work%{name}=[({400+i*20}._r8+j,j=1,size(mt_work%{name}))]']
 else:initial.append(f'mt_work%{name}={400+i*20}._r8')
initial.append('call water_diag%reset(nl_soil,lbsn)')
for i,name in enumerate(diag_map.values()):
 if name=='permeable_soil':initial.append('water_diag%permeable_soil=[.true.,.false.,.true.,.false.]')
 elif name in ['wblc_ice_sink','etroot_actual','rsub_source_layer','snow_qout_layer']:
  initial.append(f'water_diag%{name}=[({700+i*20}._r8+j,j=1,size(water_diag%{name}))]')
 else:initial.append(f'water_diag%{name}={700+i*20}._r8')
# The legacy reference below checks physics argument mapping; stage preparation
# and ordering are compared with frozen MAIN by check_stage_consolidation.py.
initial += ['mt_work%irrigation=qflx_irrig_drip+qflx_irrig_flood+qflx_irrig_paddy',
 'call mt_land_water_fluxes(mt_work,water_diag,wa,qinfl,deltim)',
 'if(.not.allocated(flood_credit_patch))allocate(flood_credit_patch(2))',
 'flood_credit_patch=[0.25_r8,0.75_r8]']
flood_services='''subroutine mt_rlc_land_flood_input(ipatch,patchtype,water_credit,evaporation,infiltration,deltim,material_input)
integer,intent(in)::ipatch,patchtype
real(r8),allocatable,intent(in)::water_credit(:)
real(r8),intent(in)::evaporation,infiltration,deltim
real(r8),intent(out)::material_input(nmts)
material_input=flood_input_material
end subroutine
subroutine mt_rlc_land_publish_input(ipatch,patchtype,dynamic_wetland,material_input)
integer,intent(in)::ipatch,patchtype
logical,intent(in)::dynamic_wetland
real(r8),intent(in)::material_input(:)
end subroutine'''
fixture='\n'.join(['#include <define.h>','module fixture','use MOD_Precision','use MOD_MT_HydroDiagnostics','use MOD_MT_LandWork','implicit none',
'integer,parameter::nl_soil=4,nmts=2','real(r8),parameter::mt_tiny=1.e-30_r8','real(r8),allocatable::flood_credit_patch(:)','integer::ipatch=1,snl,patchtype,lbsn,output_unit=20',
'logical::DEF_USE_Dynamic_Wetland,DEF_SPLIT_SOILSNOW',
'real(r8)::'+','.join(scalars),*['real(r8)::'+n+'('+d+')' for n,d in arrays.items()],
 'type(mt_land_work)::mt_work','type(mt_water_diagnostics)::water_diag','contains',stub('mt_wetland'),stub('mt_soil_water'),wrapper,routine(land,'mt_land_water_fluxes'),flood_services,
 'subroutine expected()',reference,'end subroutine','subroutine current()',actual,'end subroutine',
 'subroutine setup()','integer::j','lbsn=min(snl+1,0)',*initial,'end subroutine','end module',
 '''program check
use fixture
implicit none
integer::s,p,d,k
open(20,file='old.bin',access='stream',form='unformatted',status='replace')
open(21,file='new.bin',access='stream',form='unformatted',status='replace')
do s=0,-3,-1
do p=0,2,2
do d=0,1
do k=0,1
snl=s;patchtype=p;DEF_USE_Dynamic_Wetland=d==1;DEF_SPLIT_SOILSNOW=k==1
call setup()
output_unit=20
call expected()
output_unit=21
call current()
enddo
enddo
enddo
enddo
close(20);close(21)
end program'''])
compiler=shutil.which('gfortran')
if not compiler:raise SystemExit('gfortran required')
with tempfile.TemporaryDirectory(prefix='colm-land-boundary-') as tmp:
 path=Path(tmp)
 (path/'fixture.F90').write_text(fixture)
 for flood in [False,True]:
  (path/'define.h').write_text('#define MaterialTransport\n'+('#define GridRiverLakeFlow\n' if flood else ''))
  subprocess.run([compiler,'-O0','-g','-cpp','-ffree-line-length-none','-fcheck=all','-I'+tmp,
   str(root/'share/MOD_Precision.F90'),str(root/'main/MATERIAL/MOD_MT_HydroDiagnostics.F90'),str(root/'main/MATERIAL/MOD_MT_ThermalDiagnostics.F90'),str(root/'main/MATERIAL/MOD_MT_LandWork.F90'),'fixture.F90','-o','check'],cwd=path,check=True)
  subprocess.run([str(path/'check')],cwd=path,check=True,timeout=10)
  assert (path/'old.bin').read_bytes()==(path/'new.bin').read_bytes(),f'land call arguments changed, flood={flood}'
  print(f'PASS 32 land boundary cases, flood={flood}: all values and optional presence bitwise identical')
