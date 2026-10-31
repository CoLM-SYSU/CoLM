#!/usr/bin/env python3
"""Run actual standard/SNICAR remapping and optionally compare frozen host bits."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HOST = ROOT / 'main/MOD_SnowLayersCombineDivide.F90'


def fixture(lb=-4, ns=3):
    calls = []
    for routine in range(1, 5):
        combine = routine in (1, 3)
        name = ('snowlayerscombine', 'snowlayersdivide',
                'SnowLayersCombine_snicar', 'SnowLayersDivide_snicar')[routine - 1]
        args = 'lb,snl,z,dz,zi,wl,wi,t'
        if combine:
            args += ',scv,depth'
        if routine >= 3:
            args += ',' + ','.join(f'aerosol(:,{i})' for i in range(1, 9))
        args += ',mt_wliq=ml,mt_wice=mi,mt_solid=ms'
        if combine:
            args += ',mt_scv=mc'
        calls.append(f'case({routine})\n call {name}({args})')
    code = '''program check_snow_remap
use MOD_Precision
use MOD_SnowLayersCombineDivide
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer,parameter::lb=-4,ns=3
integer::routine,n,scenario,mask,maxmask,snl,j,k,unit,cases,repetition,i
real(r8)::z(lb:1),dz(lb:1),zi(lb-1:1),wl(lb:1),wi(lb:1),t(lb:1)
real(r8)::aerosol(lb:0,8),scv,depth,minima(5),maxima(5),before(ns),after(ns)
real(r8),allocatable::ml(:,:),mi(:,:),ms(:,:),mc(:),values(:)
minima=[0.010_r8,0.015_r8,0.025_r8,0.055_r8,0.115_r8]
maxima=[0.03_r8,0.07_r8,0.18_r8,0.41_r8,0.23_r8]
open(newunit=unit,file='results.bin',access='stream',form='unformatted',status='replace')
cases=0
do routine=1,4
 maxmask=7
 if(mod(routine,2)==1) maxmask=15
 do n=1,5
 do scenario=1,13
 do mask=0,maxmask
  snl=-n
  z=-9._r8; dz=0._r8; zi=-8._r8; wl=0._r8; wi=0._r8; t=260._r8
  aerosol=0._r8; scv=0._r8; depth=0._r8
  if(btest(mask,0)) allocate(ml(ns,lb:1))
  if(btest(mask,1)) allocate(mi(ns,lb:1))
  if(btest(mask,2)) allocate(ms(ns,lb:1))
  if(btest(mask,3)) allocate(mc(ns))
  if(allocated(ml)) ml=17._r8
  if(allocated(mi)) mi=19._r8
  if(allocated(ms)) ms=23._r8
  if(allocated(mc)) mc=0._r8
  do k=1,n
   j=k-n
   dz(j)=2._r8*minima(k)
   select case(scenario)
   case(2)
    dz(j)=0.009_r8/n
   case(3)
    dz(j)=0.6_r8
   case(8)
    dz(j)=nearest(minima(k),-1._r8)
   case(9)
    dz(j)=minima(k)
   case(10)
    dz(j)=nearest(minima(k),1._r8)
   case(11)
    dz(j)=maxima(k)
   case(12)
    dz(j)=nearest(maxima(k),1._r8)
   case(13)
    dz(j)=nearest(maxima(k),-1._r8)
   end select
   wi(j)=10._r8+k; wl(j)=0.1_r8*k; t(j)=258._r8+2._r8*k
   if((scenario==4.and.k==1).or.(scenario==5.and.k==n).or. &
      (scenario==6.and.k==(n+1)/2).or.scenario==7) wi(j)=0.05_r8
   if(allocated(ml)) ml(:,j)=[(0.1_r8*i*i*k,i=1,ns)]
   if(allocated(mi)) mi(:,j)=[(2.1_r8*i*k*k,i=1,ns)]
   if(allocated(ms)) ms(:,j)=[(0.07_r8*i*i*(n-k+1),i=1,ns)]
   aerosol(j,:)=[(0.001_r8*k*j*j, j=1,8)]
  enddo
  wi(1)=42._r8; wl(1)=17._r8; dz(1)=0.1_r8; t(1)=270._r8
  scv=sum(wi(lb:0)+wl(lb:0)); depth=sum(dz(lb:0))
  do repetition=1,2
  before=0._r8
  if(allocated(ml).and.allocated(mi).and.allocated(ms)) &
   before=sum(ml(:,snl+1:1)+mi(:,snl+1:1)+ms(:,snl+1:1),dim=2)
  if(allocated(mc)) before=before+mc
  select case(routine)
''' + '\n'.join(calls) + '''
  end select
  values=[z,dz,zi,wl,wi,t,scv,depth,reshape(aerosol,[size(aerosol)])]
  if(allocated(ml)) values=[values,reshape(ml,[size(ml)])]
  if(allocated(mi)) values=[values,reshape(mi,[size(mi)])]
  if(allocated(ms)) values=[values,reshape(ms,[size(ms)])]
  if(allocated(mc)) values=[values,mc]
  if(.not.all(ieee_is_finite(values))) error stop 'nonfinite remap output'
  ! Check complete material budgets where the thin-snow destination is supplied.
  if(mask==15.or.(mod(routine,2)==0.and.mask==7)) then
   after=sum(ml(:,snl+1:1)+mi(:,snl+1:1)+ms(:,snl+1:1),dim=2)
   if(allocated(mc)) after=after+mc
   if(any(abs(after-before)>1.e-12_r8*max(1._r8,abs(before)))) then
    print *,routine,n,scenario,mask,before,after
    error stop 'snow material conservation'
   endif
  endif
  write(unit) routine,n,scenario,mask,repetition,snl,size(values),values
  deallocate(values)
  cases=cases+1
  enddo
  if(allocated(ml)) deallocate(ml)
  if(allocated(mi)) deallocate(mi)
  if(allocated(ms)) deallocate(ms)
  if(allocated(mc)) deallocate(mc)
 enddo
 enddo
 enddo
enddo
close(unit)
print '(A,I0,A)', 'PASS ',cases,' snow remap cases: 4 routines, 1-5 layers, 13 boundaries, all optional combinations'
end program
'''
    return code.replace('lb=-4,ns=3', f'lb={lb},ns={ns}')


def run(source, target, compiler, lb=-4, ns=3):
    target.mkdir()
    (target / 'define.h').write_text('#define MaterialTransport\n')
    (target / 'host.F90').write_text(source.read_text())
    (target / 'fixture.F90').write_text(fixture(lb, ns))
    sources = [ROOT / 'share/MOD_Precision.F90', ROOT / 'main/MOD_Const_Physical.F90']
    if 'MOD_MT_SnowRemap'.lower() in source.read_text().lower():
        sources.append(ROOT / 'main/MATERIAL/MOD_MT_SnowRemap.F90')
    subprocess.run([compiler, '-O0', '-g', '-cpp', '-ffree-line-length-none',
                    '-fdefault-real-8', '-fdefault-double-8', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-I' + str(target),
                    *map(str, sources), 'host.F90', 'fixture.F90', '-o', 'check'],
                   cwd=target, check=True)
    subprocess.run([str(target / 'check')], cwd=target, check=True, timeout=30)
    return (target / 'results.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, help='frozen pre-refactor host source')
    args = parser.parse_args()
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='colm-snow-remap-') as tmp:
        path = Path(tmp)
        for lb, ns in ((-4, 3), (-5, 3), (-4, 0)):
            suffix = f'{lb}_{ns}'
            current = run(HOST, path / ('current' + suffix), compiler, lb, ns)
            if args.baseline:
                baseline = run(args.baseline.resolve(), path / ('baseline' + suffix), compiler, lb, ns)
                if current != baseline:
                    offset = next((i for i, pair in enumerate(zip(current, baseline))
                                   if pair[0] != pair[1]), min(len(current), len(baseline)))
                    raise SystemExit(f'FAIL snow remap differs at byte {offset}; '
                                     f'lb={lb} species={ns}, current={len(current)} baseline={len(baseline)} bytes')
                print(f'PASS frozen baseline bits: lb={lb}, species={ns} ({len(current)} bytes)')



if __name__ == '__main__':
    main()
