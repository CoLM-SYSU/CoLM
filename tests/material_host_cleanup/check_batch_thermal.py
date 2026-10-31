#!/usr/bin/env python3
"""Compare actual phase-change routines and thermal observations to frozen phase 3."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = r'''#include <define.h>
program thermal_regression
use MOD_Precision
use MOD_Namelist
use MOD_PhaseChange
use MOD_Hydro_SoilFunction
use MOD_Const_Physical, only: tfrz
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifndef LEGACY
use MOD_MT_ThermalDiagnostics
#endif
implicit none
integer,parameter::n=3
integer::lb,ib,kind,ptype,super,split,scenario,output,isdry,j,cases,ps,pe
integer,allocatable::imelt(:)
real(r8),allocatable::t(:),before(:),wi(:),wl(:),ice0(:),mass0(:),fact(:),brr(:),sabg(:),thaw(:),freeze(:)
real(r8)::por(n),psi(n),bsw(n),theta(n),alpha(n),nv(n),lv(n),sc(n),fc(n),dz(n)
real(r8)::scv,depth,sm,xmf,hs,signs,expected(3),melt_p(5),freeze_p(5),raw_p(5),weight(5)
#ifndef LEGACY
type(mt_thermal_diagnostics)::diag
type(mt_canopy_diagnostics)::canopy
type(mt_canopy_diagnostics),allocatable::pft(:)
! Explicit absent forwarding and resize/reset preserve canopy across layer reset.
call mt_thermal_reset(lb=-4,nl_soil=n)
call mt_phase_reset(lb=0,nl_soil=n)
call mt_phase_record(layer=1,ice_before=3._r8,ice_after=4._r8)
call mt_thermal_reset(diag,-4,n)
diag%canopy%melt=9._r8;diag%canopy%freeze=7._r8;diag%canopy%raw=5._r8
call mt_phase_reset(diag,1,n)
if(lbound(diag%soil_thaw,1)/=1.or.ubound(diag%soil_thaw,1)/=n) error stop 'phase resize'
if(diag%canopy%melt/=9._r8.or.diag%canopy%freeze/=7._r8.or.diag%canopy%raw/=5._r8) error stop 'canopy overwritten'
call mt_thermal_reset(diag,0,n)
if(any(diag%soil_thaw/=0._r8).or.any(diag%soil_freeze/=0._r8)) error stop 'phase reset'
if(diag%canopy%melt/=0._r8.or.diag%canopy%freeze/=0._r8.or.diag%canopy%raw/=0._r8) error stop 'thermal reset'
#endif
open(20,file='thermal.bin',access='stream',form='unformatted',status='replace')
! Scalar and PFT observations include values that expose addition-order changes.
do ps=1,3
 pe=ps+4
 melt_p=[1.e16_r8,1._r8,-1.e16_r8,3._r8,0.25_r8]
 freeze_p=[0.01_r8,1.e-20_r8,7._r8,2._r8,0.03_r8]
 raw_p=[11._r8,700._r8,0._r8,91._r8,3._r8]
 weight=[0.25_r8,0.125_r8,0.25_r8,0.25_r8,0.125_r8]
 expected=[sum(melt_p*weight),sum(freeze_p*weight),sum(raw_p*weight)]
#ifndef LEGACY
 allocate(pft(ps:pe))
 if(any(pft%melt/=0._r8).or.any(pft%freeze/=0._r8).or.any(pft%raw/=0._r8)) error stop 'PFT defaults'
 pft%melt=melt_p;pft%freeze=freeze_p;pft%raw=raw_p
 call mt_canopy_average(ps,pe,pft,weight,canopy)
 if(any(transfer([canopy%melt,canopy%freeze,canopy%raw],[0_8,0_8,0_8])/=transfer(expected,[0_8,0_8,0_8]))) &
  error stop 'PFT SUM bits'
 write(20) canopy%melt,canopy%freeze,canopy%raw
 deallocate(pft)
#else
 write(20) expected
#endif
enddo
por=0.45_r8;psi=-100._r8;bsw=4._r8;theta=0.05_r8
alpha=0.002_r8;nv=1.5_r8;lv=0.5_r8;dz=0.1_r8
#ifdef vanGenuchten_Mualem_SOIL_MODEL
do j=1,n
 call get_derived_parameters_vGM(psi(j),alpha(j),nv(j),sc(j),fc(j))
enddo
#endif
cases=0
do ib=1,3
 select case(ib)
 case(1);lb=-4
 case(2);lb=0
 case(3);lb=1
 end select
 allocate(t(lb:n),before(lb:n),wi(lb:n),wl(lb:n),ice0(lb:n),mass0(lb:n), &
  fact(lb:n),brr(lb:n),sabg(lb:1),thaw(lb:n),freeze(lb:n),imelt(lb:n))
 do kind=1,2
 do ptype=0,4,2
 do isdry=0,1
 do super=0,1
 do split=0,1
 do scenario=1,6
 do output=0,1
  DEF_USE_SUPERCOOL_WATER=super==1;DEF_SPLIT_SOILSNOW=split==1
  scv=0._r8;depth=0._r8
  if(lb==1.and.scenario==1) then
   scv=0.1_r8;depth=0.001_r8
  endif
  hs=0._r8;fact=0.01_r8;sabg=7._r8
  do j=lb,n
   signs=1._r8
   if(scenario==2.or.scenario==6) signs=-1._r8
   if(scenario==3.and.mod(j-lb,2)==0) signs=-1._r8
   if(scenario==4) signs=0._r8
   t(j)=tfrz+signs*2._r8;before(j)=t(j);brr(j)=signs*150._r8
   wi(j)=7._r8+j-lb;wl(j)=40._r8+j-lb
   if(scenario==5) wi(j)=0._r8
   if(scenario==6) wl(j)=0._r8
  enddo
  ice0=wi;mass0=wi+wl;thaw=0._r8;freeze=0._r8
  if(kind==1) then
   if(output==1) then
    call meltf(ptype,isdry==1,lb,n,1800._r8,fact,brr,hs,hs,hs,0.4_r8,-2._r8, &
     before,t,wl,wi,imelt,scv,depth,sm,xmf,por,psi, &
#ifdef Campbell_SOIL_MODEL
     bsw, &
#else
     theta,alpha,nv,lv,sc,fc, &
#endif
     dz, &
#ifdef LEGACY
     qphs_thaw_lay=thaw,qphs_frzc_lay=freeze)
#else
     thermal_diag=diag)
#endif
   else
    call meltf(ptype,isdry==1,lb,n,1800._r8,fact,brr,hs,hs,hs,0.4_r8,-2._r8, &
     before,t,wl,wi,imelt,scv,depth,sm,xmf,por,psi, &
#ifdef Campbell_SOIL_MODEL
     bsw, &
#else
     theta,alpha,nv,lv,sc,fc, &
#endif
     dz)
   endif
  else
   if(output==1) then
    call meltf_snicar(ptype,isdry==1,lb,n,1800._r8,fact,brr,hs,hs,hs,0.4_r8,sabg,-2._r8, &
     before,t,wl,wi,imelt,scv,depth,sm,xmf,por,psi, &
#ifdef Campbell_SOIL_MODEL
     bsw, &
#else
     theta,alpha,nv,lv,sc,fc, &
#endif
     dz, &
#ifdef LEGACY
     qphs_thaw_lay=thaw,qphs_frzc_lay=freeze)
#else
     thermal_diag=diag)
#endif
   else
    call meltf_snicar(ptype,isdry==1,lb,n,1800._r8,fact,brr,hs,hs,hs,0.4_r8,sabg,-2._r8, &
     before,t,wl,wi,imelt,scv,depth,sm,xmf,por,psi, &
#ifdef Campbell_SOIL_MODEL
     bsw, &
#else
     theta,alpha,nv,lv,sc,fc, &
#endif
     dz)
   endif
  endif
#ifndef LEGACY
  if(output==1) then
   thaw=diag%soil_thaw;freeze=diag%soil_freeze
  endif
#endif
  if(.not.all(ieee_is_finite([t,wi,wl,scv,depth,sm,xmf,thaw,freeze]))) error stop 'finite thermal outputs'
  if(any(abs(wi+wl-mass0)>1.e-12_r8)) error stop 'phase water conservation'
  if(output==1) then
   if(any(thaw/=max(ice0-wi,0._r8)).or.any(freeze/=max(wi-ice0,0._r8))) error stop 'phase observations'
  endif
  write(20) lb,kind,ptype,isdry,super,split,scenario,output,imelt,t,wi,wl,scv,depth,sm,xmf,thaw,freeze
  cases=cases+1
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 deallocate(t,before,wi,wl,ice0,mass0,fact,brr,sabg,thaw,freeze,imelt)
enddo
close(20)
print '(A,I0,A)', 'PASS ',cases,' actual thermal phase calls, canopy SUM/reset/bounds'
end program
'''


def run(path, baseline, model, optimization):
    path.mkdir()
    (path / 'define.h').write_text(f'#define MaterialTransport\n#define {model}\n')
    (path / 'stubs.F90').write_text('''module MOD_SPMD_Task
contains
subroutine CoLM_stop()
error stop 'host phase water balance'
end subroutine
end module
module MOD_Namelist
logical :: DEF_USE_SUPERCOOL_WATER, DEF_SPLIT_SOILSNOW
end module
''')
    (path / 'fixture.F90').write_text(FIXTURE)
    sources = [ROOT / 'share/MOD_Precision.F90', path / 'stubs.F90',
               ROOT / 'main/MOD_Const_Physical.F90', ROOT / 'main/HYDRO/MOD_Hydro_SoilFunction.F90']
    flags = [optimization, '-g', '-cpp', '-ffree-line-length-none', '-fdefault-real-8',
             '-fdefault-double-8', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-I' + str(path)]
    if baseline:
        flags += ['-DLEGACY']
        sources += [baseline / 'main/MOD_PhaseChange.F90']
    else:
        sources += [ROOT / 'main/MATERIAL/MOD_MT_ThermalDiagnostics.F90', ROOT / 'main/MOD_PhaseChange.F90']
    sources += [path / 'fixture.F90']
    subprocess.run(['gfortran', *flags, *map(str, sources), '-o', 'check'], cwd=path, check=True)
    subprocess.run([str(path / 'check')], cwd=path, check=True, timeout=30)
    return (path / 'thermal.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, help='frozen pre-phase-4 source tree')
    args = parser.parse_args()
    if not shutil.which('gfortran'):
        raise SystemExit('gfortran is required')
    baseline = args.baseline.resolve() if args.baseline else None
    with tempfile.TemporaryDirectory(prefix='colm-batch-thermal-') as tmp:
        for index, model in enumerate(['Campbell_SOIL_MODEL', 'vanGenuchten_Mualem_SOIL_MODEL']):
            for optimization in ['-O0', '-O2']:
                label = f'{index}-{optimization[1:]}'
                current = run(Path(tmp) / ('current-' + label), None, model, optimization)
                if baseline:
                    old = run(Path(tmp) / ('baseline-' + label), baseline, model, optimization)
                    assert current == old, f'thermal host/diagnostic bits differ: {model} {optimization}'
                    print(f'PASS thermal baseline bits: {model} {optimization}, {len(current)} bytes')


if __name__ == '__main__':
    main()
