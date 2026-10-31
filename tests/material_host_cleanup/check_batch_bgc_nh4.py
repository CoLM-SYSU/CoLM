#!/usr/bin/env python3
"""Exercise actual NH4 helpers; --baseline compares frozen BGC statement blocks."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

DECL = '''
integer,parameter::nl_soil=4,numpatch=2
real(r8)::smin_nh4_vr(nl_soil,numpatch),wliq_soisno(nl_soil,numpatch),dz_soi(nl_soil)
real(r8)::f_n2o_nit(numpatch),smin_no3_leached(numpatch),smin_no3_runoff(numpatch)
'''

def old_block(path, marker):
    text = path.read_text()
    start = text.index(marker)
    return text[start:text.index('#endif', start)]


def model(baseline):
    if baseline:
        bgc = baseline / 'main/BGC'
        leach = old_block(bgc / 'MOD_BGC_Soil_BiogeochemNLeaching.F90',
                          '            ! Species-resolved NH4 leaching')
        update = old_block(bgc / 'MOD_BGC_CNNStateUpdate3.F90',
                           '            ! Apply the optional NH4 leaching loss')
        summary = old_block(bgc / 'MOD_BGC_CNSummary.F90', '         smin_nh4_leached(i)')
        balance = (bgc / 'MOD_BGC_CNBalanceCheck.F90').read_text()
        start = balance.index('         col_noutputs = col_noutputs + f_n2o_nit(i) + smin_nh4_leached(i)')
        balance = balance[start:balance.index('#else', start)]
        imports = '''use MOD_Precision
use MOD_MT_Namelist
implicit none
real(r8),allocatable::smin_nh4_leached_vr(:,:),smin_nh4_leached(:)
'''
    else:
        imports = 'use MOD_Precision\nuse MOD_MT_Namelist\nuse MOD_MT_BGC_NH4\nimplicit none\n'
        leach = 'call mt_nh4_leach_layer(j,i,sf,smin_nh4_vr(j,i),dz_soi(j),wliq_soisno(j,i),tot_water,drain_tot,deltim)'
        update = 'call mt_nh4_update_layer(j,i,deltim,smin_nh4_vr(j,i))'
        summary = 'call mt_nh4_sum_layer(j,i,dz_soi(j))'
        balance = 'call mt_nh4_balance_outputs(i,col_noutputs,f_n2o_nit(i),smin_no3_leached(i),smin_no3_runoff(i))'
    return 'module model\n' + imports + DECL + '''contains
subroutine leach(j,i,sf,tot_water,drain_tot,deltim)
integer::j,i
real(r8)::sf,tot_water,drain_tot,deltim,disn_conc
''' + leach + '''
end subroutine
subroutine update(j,i,deltim)
integer::j,i
real(r8)::deltim
''' + update + '''
end subroutine
subroutine summary(j,i)
integer::j,i
''' + summary + '''
end subroutine
subroutine balance(i,col_noutputs)
integer::i
real(r8)::col_noutputs
''' + balance + '\nend subroutine\nend module\n'

FIXTURE = '''program check
use model
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
integer::enabled,wet,total,drain,pool,j,i,n,unit
real(r8)::sf,tot_water,drain_tot,deltim,col_noutputs,old_pool,bound
real(r8),parameter::drains(3)=[-2._r8,0.03_r8,100._r8],pools(3)=[-3._r8,0._r8,1.7_r8]
allocate(smin_nh4_leached_vr(nl_soil,numpatch),smin_nh4_leached(numpatch))
open(newunit=unit,file='values.bin',access='stream',form='unformatted',status='replace')
n=0
sf=0.1_r8;deltim=1800._r8
 do enabled=0,1
 do wet=0,2
 do total=0,1
 do drain=1,3
 do pool=1,3
  DEF_USE_NH4_LEACHING=enabled==1
  smin_nh4_leached_vr=-9._r8;smin_nh4_leached=0._r8
  do i=1,numpatch
   do j=1,nl_soil
    dz_soi(j)=0.03_r8*j
    wliq_soisno(j,i)=wet*0.7_r8*(j+i)
    smin_nh4_vr(j,i)=pools(pool)*(j+i)
   enddo
  enddo
  tot_water=total*sum(wliq_soisno(:,1));drain_tot=drains(drain)
  do i=1,numpatch
   do j=1,nl_soil
    old_pool=smin_nh4_vr(j,i)
    call leach(j,i,sf,tot_water,drain_tot,deltim)
    if(smin_nh4_leached_vr(j,i)<0._r8) error stop 'negative NH4 loss'
    bound=sf*max(old_pool,0._r8)/deltim
    if(smin_nh4_leached_vr(j,i)>bound) error stop 'NH4 soluble cap'
    if(.not.DEF_USE_NH4_LEACHING.and.smin_nh4_leached_vr(j,i)/=0._r8) error stop 'disabled NH4 loss'
    call update(j,i,deltim)
    if(.not.DEF_USE_NH4_LEACHING.and.smin_nh4_vr(j,i)/=old_pool) error stop 'disabled NH4 state'
    call summary(j,i)
   enddo
   f_n2o_nit(i)=-1.e16_r8
   smin_no3_leached(i)=0.73_r8;smin_no3_runoff(i)=-0.31_r8
   col_noutputs=1.e16_r8
   call balance(i,col_noutputs)
   if(.not.ieee_is_finite(col_noutputs)) error stop 'nonfinite balance'
   write(unit)smin_nh4_leached_vr(:,i),smin_nh4_vr(:,i),smin_nh4_leached(i),col_noutputs
  enddo
  n=n+1
 enddo
 enddo
 enddo
 enddo
 enddo
close(unit)
print *, 'PASS NH4 scenarios:',n
end program
'''


def run(path, compiler, baseline=None):
    path.mkdir()
    (path/'define.h').write_text('#define MaterialTransport\n#define BGC\n')
    (path/'settings.F90').write_text('module MOD_MT_Namelist\nlogical::DEF_USE_NH4_LEACHING\nend module\n')
    (path/'model.F90').write_text(model(baseline))
    (path/'fixture.F90').write_text(FIXTURE)
    sources = [ROOT/'share/MOD_Precision.F90', path/'settings.F90']
    if not baseline:
        sources.append(ROOT/'main/MATERIAL/FAMILY_SOLUTE/MOD_MT_BGC_NH4.F90')
    sources += [path/'model.F90', path/'fixture.F90']
    subprocess.run([compiler,'-O0','-g','-cpp','-ffree-line-length-none','-fdefault-real-8',
                    '-fdefault-double-8','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                    '-I'+str(path),*map(str,sources),'-o','check'],cwd=path,check=True)
    subprocess.run([str(path/'check')],cwd=path,check=True,timeout=20)
    return (path/'values.bin').read_bytes()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,help='frozen full pre-phase-4 source tree')
    args=parser.parse_args()
    compiler=shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='colm-batch-nh4-') as tmp:
        current=run(Path(tmp)/'current',compiler)
        if args.baseline:
            previous=run(Path(tmp)/'baseline',compiler,args.baseline.resolve())
            assert current==previous,'NH4 flux/state/sum/balance arithmetic changed'
            print('PASS frozen NH4 statement bits:',len(current),'bytes')


if __name__=='__main__':
    main()
