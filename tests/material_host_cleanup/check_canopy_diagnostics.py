#!/usr/bin/env python3
"""Compare real canopy/PFT interception and drainage against a frozen host source."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HOST = ROOT / 'main/MOD_LeafInterception.F90'

STUBS = '''module MOD_Namelist
logical :: DEF_VEG_SNOW
integer :: DEF_Interception_scheme
end module
module MOD_LandPFT
use MOD_Precision
integer :: patch_pft_s(1)=[1],patch_pft_e(1)=[4],pftclass(4)=[1,4,9,2]
real(r8) :: pftfrac(4)=[0.137_r8,0.503_r8,0.223_r8,0.137_r8]
end module
module MOD_Vars_PFTimeInvariants
use MOD_Precision
real(r8) :: ncd_p(4),ncw_p(4),bcw_p(4),htop_p(4)
end module
module MOD_Vars_PFTimeVariables
use MOD_Precision
real(r8) :: sigf_p(4),lai_p(4),sai_p(4),tleaf_p(4),z0m_p(4)
real(r8) :: ldew_p(4),ldew_rain_p(4),ldew_snow_p(4)
end module
module MOD_Vars_1DPFTFluxes
use MOD_Precision
real(r8) :: qintr_p(4),qintr_rain_p(4),qintr_snow_p(4)
end module
module MOD_Const_PFT
use MOD_Precision
real(r8) :: chil_p(17)=0.1_r8
end module
'''

FIXTURE = '''program check
use MOD_Precision
use MOD_LeafInterception
use MOD_Namelist
use MOD_LandPFT
use MOD_Vars_PFTimeInvariants
use MOD_Vars_PFTimeVariables
use MOD_Vars_1DPFTFluxes
use drainage_adapter
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer :: route,scheme,snow,dtcase,scenario,step,unit,cases,k
real(r8) :: dt,dewmx,lai,sai,leaf,store,rain_store,snow_store
real(r8) :: pg_r,pg_s,qi,qr,qs,prc_r,prc_s,prl_r,prl_s,sprinkler,old_store
real(r8) :: values(40),times(3)=[1._r8,60._r8,1800._r8],drain(2),txr,txs
open(newunit=unit,file='canopy.bin',access='stream',form='unformatted',status='replace')
cases=0
do route=1,2
do scheme=1,8,7
do snow=0,1
do dtcase=1,3
do scenario=1,4
 DEF_Interception_scheme=scheme
 DEF_VEG_SNOW=snow==1
 dt=times(dtcase);dewmx=0.1_r8;leaf=275._r8;lai=2._r8;sai=0.2_r8
 prc_r=0.001_r8;prl_r=0.0007_r8;prc_s=0.0003_r8;prl_s=0.0002_r8;sprinkler=0.00013_r8
 rain_store=0.04_r8*scenario;snow_store=0.3_r8*(scenario-1)
 if(scenario==2) then
  prc_r=0._r8;prl_r=0._r8;prc_s=0._r8;prl_s=0._r8;sprinkler=0._r8
 endif
 if(scenario==3) then
  rain_store=1.2_r8;snow_store=12._r8;leaf=268._r8
 endif
 if(scenario==4) then
  rain_store=0._r8;snow_store=0._r8
 endif
 store=rain_store+snow_store
 do k=1,4
  sigf_p(k)=0.8_r8;lai_p(k)=0.5_r8+k;sai_p(k)=0.1_r8*k
  tleaf_p(k)=leaf+k-2;z0m_p(k)=0.2_r8+0.1_r8*k
  ncd_p(k)=3._r8+k;ncw_p(k)=3._r8+0.2_r8*k;bcw_p(k)=2._r8+k;htop_p(k)=10._r8+k
  ldew_rain_p(k)=rain_store*(0.3_r8+k);ldew_snow_p(k)=snow_store*(0.2_r8+k)
 enddo
 ldew_p=ldew_rain_p+ldew_snow_p
 qintr_p=0._r8;qintr_rain_p=0._r8;qintr_snow_p=0._r8
 do step=1,3
  if(step==3) then
   lai=0._r8;sai=0._r8;lai_p=0._r8;sai_p=0._r8
  endif
  old_store=store
  if(route==2) old_store=sum(ldew_p*pftfrac)
  call invoke(route,dt)
  drain=read_drainage()
  if(step==3.and.any(drain/=0._r8)) error stop 'old-canopy release must clear drainage'
  if(abs(store+(pg_r+pg_s)*dt-old_store-(prc_r+prl_r+prc_s+prl_s+sprinkler)*dt)>1.e-10_r8) &
   error stop 'canopy water budget'
  call record()
 enddo
enddo
enddo
enddo
enddo
enddo
! Zero time is defined for a bare canopy with no old store; vegetated physics
! itself divides by dt, so only the moved capture guard is tested at zero below.
do route=1,2
 DEF_Interception_scheme=1;DEF_VEG_SNOW=.true.;dt=0._r8
 lai=0._r8;sai=0._r8;lai_p=0._r8;sai_p=0._r8
 store=0._r8;rain_store=0._r8;snow_store=0._r8
 ldew_p=0._r8;ldew_rain_p=0._r8;ldew_snow_p=0._r8
 call seed_drainage(7._r8,9._r8)
 call invoke(route,dt)
 drain=read_drainage()
 if(any(drain/=0._r8)) error stop 'bare zero-time drainage'
 call record()
enddo
do dtcase=1,3
 dt=real(dtcase-2,r8)*60._r8
 do scenario=1,3
  txr=real(scenario-2,r8);txs=-txr
  call seed_drainage(7._r8,9._r8)
  call capture_guard(dt,txr,txs)
  drain=read_drainage()
  if(dt<=0._r8.and.any(drain/=0._r8)) error stop 'nonpositive time guard'
  if(any(drain<0._r8)) error stop 'negative drainage'
  write(unit) drain
 enddo
enddo
close(unit)
print '(A,I0,A)', 'PASS ',cases,' actual canopy/PFT calls and 9 capture guards'
contains
subroutine invoke(which,timestep)
 integer,intent(in)::which
 real(r8),intent(in)::timestep
 if(which==1) then
  call LEAF_interception_wrap(timestep,dewmx,3._r8,4._r8,0.1_r8,0.8_r8,lai,sai,270._r8,leaf, &
   prc_r,prc_s,prl_r,prl_s,100._r8,1,5._r8,4._r8,3._r8,15._r8, &
   store,rain_store,snow_store,0.3_r8,10._r8,pg_r,pg_s,qi,qr,qs,sprinkler)
 else
  call LEAF_interception_pftwrap(1,timestep,dewmx,3._r8,4._r8,270._r8, &
   prc_r,prc_s,prl_r,prl_s,100._r8, &
   store,rain_store,snow_store,0.3_r8,10._r8,pg_r,pg_s,qi,qr,qs,sprinkler)
 endif
end subroutine
subroutine record()
 values=[store,rain_store,snow_store,leaf,pg_r,pg_s,qi,qr,qs,drain, &
         ldew_p,ldew_rain_p,ldew_snow_p,tleaf_p,qintr_p,qintr_rain_p,qintr_snow_p,0._r8]
 if(.not.all(ieee_is_finite(values))) error stop 'nonfinite canopy output'
 write(unit) values
 cases=cases+1
end subroutine
end program
'''


def adapter(source, legacy):
    if legacy:
        fields = ('mt_leaf_qdrp_rain', 'mt_leaf_qdrp_snow')
        imports = 'use MOD_LeafInterception, only: ' + ','.join(fields)
        capture = re.search(r'         IF \(deltim > 0\._r8\) THEN\n.*?         ENDIF',
                            source, re.S).group()
    else:
        fields = ('mt_canopy_drain%rain', 'mt_canopy_drain%snow')
        imports = 'use MOD_MT_CanopyDiagnostics, only: mt_canopy_drain'
        capture = 'call mt_canopy_drain%capture(deltim,tex_rain,tex_snow)'
    return f'''module drainage_adapter
use MOD_Precision
{imports}
implicit none
contains
function read_drainage() result(values)
real(r8)::values(2)
values=[{fields[0]},{fields[1]}]
end function
subroutine seed_drainage(rain,snow)
real(r8),intent(in)::rain,snow
{fields[0]}=rain
{fields[1]}=snow
end subroutine
subroutine capture_guard(deltim,tex_rain,tex_snow)
real(r8),intent(in)::deltim,tex_rain,tex_snow
{capture}
end subroutine
end module
'''


def run(source, directory, compiler):
    directory.mkdir()
    text = source.read_text()
    legacy = 'USE MOD_MT_CanopyDiagnostics' not in text
    (directory / 'define.h').write_text('#define MaterialTransport\n#define LULC_IGBP_PFT\n#define CROP\n')
    (directory / 'stubs.F90').write_text(STUBS)
    (directory / 'host.F90').write_text(text)
    (directory / 'adapter.F90').write_text(adapter(text, legacy))
    (directory / 'check.F90').write_text(FIXTURE)
    sources = [ROOT / 'share/MOD_Precision.F90', ROOT / 'main/MOD_Const_Physical.F90',
               directory / 'stubs.F90']
    if not legacy:
        sources.append(ROOT / 'main/MATERIAL/MOD_MT_CanopyDiagnostics.F90')
    sources += [directory / f for f in ['host.F90', 'adapter.F90', 'check.F90']]
    subprocess.run([compiler, '-O0', '-g', '-cpp', '-ffree-line-length-none',
                    '-fdefault-real-8', '-fdefault-double-8', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-I'+str(directory),
                    *map(str, sources), '-o', 'check'], cwd=directory, check=True)
    subprocess.run([str(directory / 'check')], cwd=directory, check=True, timeout=30)
    return (directory / 'canopy.bin').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='colm-canopy-diagnostics-') as tmp:
        path = Path(tmp)
        current = run(HOST, path / 'current', compiler)
        if args.baseline:
            previous = run(args.baseline.resolve(), path / 'baseline', compiler)
            assert current == previous, 'canopy water or drainage output changed'
            print(f'PASS frozen canopy baseline: {len(current)} bytes exactly equal')


if __name__ == '__main__':
    main()
