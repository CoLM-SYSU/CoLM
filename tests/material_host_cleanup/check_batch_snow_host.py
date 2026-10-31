#!/usr/bin/env python3
"""Check real snow adapters and optional stage consolidation against frozen host calls."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
STUBS = '''module MOD_MT_Defs
integer :: nmts=0
end module
module MOD_MT_Vars
use MOD_Precision
real(r8),allocatable :: mt_wliq_soisno(:,:,:),mt_wice_soisno(:,:,:),mt_solid_soisno(:,:,:),mt_scv(:,:)
end module
'''
FIXTURE = '''program snow_host_regression
use MOD_Precision
use MOD_MT_Defs
use MOD_MT_Vars
use MOD_MT_SnowHost
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer::lb,ib,ns,n,snl,kind,snicar,scenario,j,k,cases
real(r8),allocatable::z(:),dz(:),zi(:),wl(:),wi(:),t(:),aerosol(:,:),values(:)
real(r8)::scv,depth
open(20,file='snow.bin',access='stream',form='unformatted',status='replace')
cases=0
do ib=1,2
 lb=-3-ib
 allocate(z(lb:1),dz(lb:1),zi(lb-1:1),wl(lb:1),wi(lb:1),t(lb:1),aerosol(lb:0,8))
 do ns=0,3
  if(ns==2) cycle
  nmts=ns
  if(ns>0) allocate(mt_wliq_soisno(ns,lb:1,2),mt_wice_soisno(ns,lb:1,2), &
   mt_solid_soisno(ns,lb:1,2),mt_scv(ns,2))
  do n=1,5
  do kind=1,2
  do snicar=0,1
  do scenario=1,8
   snl=-n;z=-9._r8;dz=0._r8;zi=-7._r8;wl=0._r8;wi=0._r8;t=260._r8
   aerosol=0._r8;scv=0._r8;depth=0._r8
   if(ns>0) then
    mt_wliq_soisno=0._r8;mt_wice_soisno=0._r8;mt_solid_soisno=0._r8;mt_scv=0._r8
    mt_wliq_soisno(:,:,1)=101._r8;mt_wice_soisno(:,:,1)=103._r8
    mt_solid_soisno(:,:,1)=107._r8;mt_scv(:,1)=109._r8
   endif
   do k=1,n
    j=k-n;dz(j)=0.03_r8*k;wl(j)=0.2_r8*k;wi(j)=11._r8+k;t(j)=258._r8+k
    if(scenario==2) dz(j)=0.009_r8/n
    if(scenario==3) dz(j)=0.6_r8
    if(scenario==4.and.k==1) wi(j)=0.05_r8
    if(scenario==5.and.k==n) wi(j)=0.05_r8
    if(scenario==6) wi(j)=0.05_r8
    if(scenario==7) dz(j)=0.01_r8
    if(scenario==8) dz(j)=nearest(0.03_r8,1._r8)
    aerosol(j,:)=[(0.001_r8*k*j,j=1,8)]
    if(ns>0) then
     mt_wliq_soisno(:,k-n,2)=[(0.13_r8*j*k,j=1,ns)]
     mt_wice_soisno(:,k-n,2)=[(0.17_r8*j*k*k,j=1,ns)]
     mt_solid_soisno(:,k-n,2)=[(0.19_r8*j*(n-k+1),j=1,ns)]
    endif
   enddo
   wl(1)=20._r8;wi(1)=40._r8;dz(1)=0.1_r8;t(1)=270._r8
   scv=sum(wl(lb:0)+wi(lb:0));depth=sum(dz(lb:0))
   if(kind==1) then
    call mt_land_snow_combine(2,snicar==1,lb,snl,z,dz,zi,wl,wi,t,scv,depth, &
     aerosol(:,1),aerosol(:,2),aerosol(:,3),aerosol(:,4), &
     aerosol(:,5),aerosol(:,6),aerosol(:,7),aerosol(:,8))
   else
    call mt_land_snow_divide(2,snicar==1,lb,snl,z(lb:0),dz(lb:0),zi(lb-1:0), &
     wl(lb:0),wi(lb:0),t(lb:0),aerosol(:,1),aerosol(:,2),aerosol(:,3),aerosol(:,4), &
     aerosol(:,5),aerosol(:,6),aerosol(:,7),aerosol(:,8))
   endif
   values=[z,dz,zi,wl,wi,t,scv,depth,reshape(aerosol,[size(aerosol)])]
   if(ns>0) then
    values=[values,reshape(mt_wliq_soisno,[size(mt_wliq_soisno)]), &
     reshape(mt_wice_soisno,[size(mt_wice_soisno)]), &
     reshape(mt_solid_soisno,[size(mt_solid_soisno)]),reshape(mt_scv,[size(mt_scv)])]
    if(any(mt_wliq_soisno(:,:,1)/=101._r8).or.any(mt_wice_soisno(:,:,1)/=103._r8).or. &
       any(mt_solid_soisno(:,:,1)/=107._r8).or.any(mt_scv(:,1)/=109._r8)) error stop 'patch isolation'
   else
    if(allocated(mt_wliq_soisno).or.allocated(mt_wice_soisno).or.allocated(mt_solid_soisno).or. &
       allocated(mt_scv)) error stop 'zero material arrays must stay unallocated'
   endif
   if(.not.all(ieee_is_finite(values))) error stop 'finite snow adapter outputs'
   write(20) lb,ns,n,kind,snicar,scenario,snl,values
   deallocate(values);cases=cases+1
  enddo
  enddo
  enddo
  enddo
  if(ns>0) deallocate(mt_wliq_soisno,mt_wice_soisno,mt_solid_soisno,mt_scv)
 enddo
 deallocate(z,dz,zi,wl,wi,t,aerosol)
enddo
close(20)
print '(A,I0,A)','PASS ',cases,' snow adapter calls: standard/SNICAR, zero/one/three materials, two bounds, isolated patches'
end program
'''


def legacy_adapter(main):
    """Extract only the old direct calls; no current adapter code is reused."""
    calls = {}
    for match in re.finditer(r'(?im)^\s*CALL (snowlayers(?:combine|divide)(?:_snicar)?)\s*\(', main):
        start = match.start()
        pos = main.index('(', start) + 1
        depth = 1
        while depth:
            depth += (main[pos] == '(') - (main[pos] == ')')
            pos += 1
        calls.setdefault(match.group(1).lower(), []).append(main[start:pos].strip())
    assert len(calls) == 4 and all(len(value) == 2 for value in calls.values())
    text = '''module MOD_MT_SnowHost
use MOD_Precision
use MOD_MT_Defs
use MOD_MT_Vars
use MOD_SnowLayersCombineDivide
implicit none
contains
'''
    aerosol = ['mss_bcpho', 'mss_bcphi', 'mss_ocpho', 'mss_ocphi',
               'mss_dst1', 'mss_dst2', 'mss_dst3', 'mss_dst4']
    for kind, upper in [('combine', 1), ('divide', 0)]:
        text += f'''subroutine mt_land_snow_{kind}(ipatch,use_snicar,lb,snl,z_soisno,dz_soisno,zi_soisno, &
 wliq_soisno,wice_soisno,t_soisno, &
'''
        if kind == 'combine':
            text += ' scv,snowdp, &\n'
        text += ','.join(aerosol) + ')\n'
        text += 'integer :: ipatch,lb,snl\nlogical :: use_snicar\n'
        text += f'real(r8) :: z_soisno(lb:{upper}),dz_soisno(lb:{upper}),zi_soisno(lb-1:{upper}), &\n'
        text += f' wliq_soisno(lb:{upper}),wice_soisno(lb:{upper}),t_soisno(lb:{upper})\n'
        if kind == 'combine':
            text += 'real(r8)::scv,snowdp\n'
        text += 'real(r8)::' + ','.join(name + '(lb:0)' for name in aerosol) + '\n'
        text += 'if(use_snicar) then\n'
        for variant in ['_snicar', '']:
            if not variant:
                text += 'else\n'
            old = calls['snowlayers' + kind + variant]
            text += 'if(nmts>0)then\n' + old[0] + '\nelse\n' + old[1] + '\nendif\n'
        text += 'endif\nend subroutine\n'
    return text + 'end module\n'


def stage_fixture(baseline):
    """Compare the new entry with the exact frozen MAIN combine/divide sequence."""
    start = FIXTURE.index('   if(kind==1) then')
    end = FIXTURE.index('   values=', start)
    if baseline:
        main = (baseline / 'main/CoLMMAIN.F90').read_text()
        a = main.index('            CALL mt_land_snow_combine(')
        b = main.index('\n#else', a)
        sequence = main[a:b]
        for old, new in [('ipatch', '2'), ('DEF_USE_SNICAR', 'snicar==1'),
                         ('z_soisno', 'z'), ('dz_soisno', 'dz'), ('zi_soisno', 'zi'),
                         ('wliq_soisno', 'wl'), ('wice_soisno', 'wi'),
                         ('t_soisno', 't'), ('snowdp', 'depth')]:
            sequence = re.sub(r'\b' + old + r'\b', new, sequence)
        for i, name in enumerate(('mss_bcpho', 'mss_bcphi', 'mss_ocpho', 'mss_ocphi',
                                  'mss_dst1', 'mss_dst2', 'mss_dst3', 'mss_dst4'), 1):
            sequence = sequence.replace(name + '(lb:0)', f'aerosol(:,{i})')
    else:
        a = FIXTURE.index('    call mt_land_snow_combine(', start)
        b = FIXTURE.index('   else', a)
        sequence = FIXTURE[a:b].replace('mt_land_snow_combine', 'mt_land_snow_normalize')
    return (FIXTURE[:start] + sequence + '\n' + FIXTURE[end:]).replace(
        'do kind=1,2', 'do kind=3,3').replace('do n=1,5', 'do n=0,5')


def run(path, baseline, stage=False):
    path.mkdir()
    (path / 'define.h').write_text('#define MaterialTransport\n')
    (path / 'stubs.F90').write_text(STUBS)
    (path / 'fixture.F90').write_text(stage_fixture(baseline) if stage else FIXTURE)
    origin = baseline or ROOT
    sources = [ROOT / 'share/MOD_Precision.F90', ROOT / 'main/MOD_Const_Physical.F90',
               path / 'stubs.F90', origin / 'main/MATERIAL/MOD_MT_SnowRemap.F90',
               origin / 'main/MOD_SnowLayersCombineDivide.F90']
    if baseline and not stage:
        (path / 'legacy.F90').write_text(legacy_adapter((baseline / 'main/CoLMMAIN.F90').read_text()))
        sources.append(path / 'legacy.F90')
    else:
        sources.append(origin / 'main/MATERIAL/MOD_MT_SnowHost.F90')
    sources.append(path / 'fixture.F90')
    subprocess.run(['gfortran', '-O0', '-g', '-cpp', '-ffree-line-length-none',
                    '-fdefault-real-8', '-fdefault-double-8', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-I' + str(path),
                    *map(str, sources), '-o', 'check'], cwd=path, check=True)
    subprocess.run([str(path / 'check')], cwd=path, check=True, timeout=30)
    return (path / 'snow.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--stage-baseline', type=Path, help='frozen tree with sequential SnowHost calls in MAIN')
    args = parser.parse_args()
    if not shutil.which('gfortran'):
        raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='colm-snow-host-') as tmp:
        current = run(Path(tmp) / 'current', None)
        old = run(Path(tmp) / 'baseline', args.baseline.resolve())
        assert current == old, 'snow adapter changed host/material/aerosol output bits'
        print(f'PASS frozen MAIN snow calls: {len(current)} bytes, zero differences')
        if args.stage_baseline:
            current = run(Path(tmp) / 'stage-current', None, stage=True)
            old = run(Path(tmp) / 'stage-baseline', args.stage_baseline.resolve(), stage=True)
            assert current == old, 'snow normalization changed sequential host/material/aerosol output bits'
            print(f'PASS frozen sequential snow stage: {len(current)} bytes, zero differences')


if __name__ == '__main__':
    main()
