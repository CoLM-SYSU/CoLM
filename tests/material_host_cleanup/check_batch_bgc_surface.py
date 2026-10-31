#!/usr/bin/env python3
"""Exercise methane clipping callbacks and verify moved aggregation source bodies."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]


def aggregate_body(text,name):
    return re.search(r'(?im)^SUBROUTINE '+name+r'\b.*?^END SUBROUTINE '+name+r'\s*$',
                     text,re.S|re.M).group().rstrip()


def legacy_module(baseline):
    text=(baseline/'mksrfdata/MOD_RegionClip.F90').read_text()
    start=text.index('               ! Methane surface data (lake-bottom soil carbon, soil pH)')
    block=text[start:text.index('#endif',start)]
    current=(ROOT/'main/MATERIAL/MOD_MT_SurfaceData.F90').read_text()
    header=current[:current.index('               ! Methane surface data')]
    return header+block+'   END SUBROUTINE mt_clip_methane_surface_data\nEND MODULE MOD_MT_SurfaceData\n#endif\n'


STUB='''module MOD_Block
implicit none
contains
subroutine get_filename_block(filename,iblk,jblk,fileblock)
character(len=*),intent(in)::filename
integer,intent(in)::iblk,jblk
character(len=*),intent(out)::fileblock
character(30)::suffix
write(suffix,'(A,I0,A,I0)')'.',iblk,'.',jblk
fileblock=trim(filename)//trim(suffix)
end subroutine
end module
'''
FIXTURE='''module recorder
implicit none
integer::unit,errors,clips
contains
subroutine clip_vector(file_in,file_out,iblk,jblk,varname,vecmask)
character(len=*),intent(in)::file_in,file_out,varname
integer,intent(in)::iblk,jblk
logical,intent(in)::vecmask(:)
clips=clips+1
write(unit)'clip',file_in,file_out,iblk,jblk,varname,vecmask
end subroutine
subroutine clip_error(message)
character(len=*),intent(in)::message
errors=errors+1
write(unit)'error',message
end subroutine
end module
program check
use MOD_Block
use MOD_MT_SurfaceData
use recorder
implicit none
integer::present_mask,required_mask,layout,f,iblk,jblk,k,expected_errors
character(len=256)::file_in,fileblock
character(len=4)::year_in,year_out
logical::mask(4),exists
open(newunit=unit,file='calls.bin',access='stream',form='unformatted',status='replace')
mask=[.true.,.false.,.true.,.false.]
do layout=1,2
year_in='2001';year_out='2010'
if(layout==2)then
 year_in='1999';year_out='1999'
endif
iblk=layout;jblk=layout+2
call system('mkdir -p input/soil/'//year_in)
do present_mask=0,3
do required_mask=0,3
 errors=0;clips=0;expected_errors=0
 do k=0,1
  if(k==0)then
   file_in='input/soil/'//year_in//'/lake_soilc_patches.nc'
  else
   file_in='input/soil/'//year_in//'/methane_ph_patches.nc'
  endif
  call get_filename_block(file_in,iblk,jblk,fileblock)
  inquire(file=trim(fileblock),exist=exists)
  if(exists)then
   open(newunit=f,file=trim(fileblock));close(f,status='delete')
  endif
  if(btest(present_mask,k))then
   open(newunit=f,file=trim(fileblock),status='new');close(f)
  else if(btest(required_mask,k))then
   expected_errors=expected_errors+1
  endif
 enddo
 write(unit)layout,present_mask,required_mask
 call mt_clip_methane_surface_data('input','output',year_in,year_out, &
  iblk,jblk,mask,btest(required_mask,0),btest(required_mask,1),clip_vector,clip_error)
 if(clips/=popcnt(present_mask)) error stop 'wrong existing-file clipping count'
 if(errors/=expected_errors) error stop 'wrong required-file error count'
enddo
enddo
enddo
close(unit)
print *, 'PASS 32 methane clipping presence/required/path scenarios'
end program
'''


def run(path,compiler,baseline=None):
    path.mkdir()
    (path/'define.h').write_text('#define MaterialTransport\n#define BGC\n')
    (path/'block.F90').write_text(STUB)
    (path/'surface.F90').write_text(legacy_module(baseline) if baseline else
                                  (ROOT/'main/MATERIAL/MOD_MT_SurfaceData.F90').read_text())
    (path/'fixture.F90').write_text(FIXTURE)
    subprocess.run([compiler,'-O0','-g','-cpp','-ffree-line-length-none','-fcheck=all',
                    '-I'+str(path),'block.F90','surface.F90','fixture.F90','-o','check'],
                   cwd=path,check=True)
    subprocess.run([str(path/'check')],cwd=path,check=True,timeout=20)
    return (path/'calls.bin').read_bytes()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,help='frozen full pre-phase-4 source tree')
    args=parser.parse_args()
    compiler=shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required')
    if args.baseline:
        args.baseline=args.baseline.resolve()
        current=(ROOT/'main/MATERIAL/FAMILY_GAS/MOD_MT_Methane_Aggregation.F90').read_text()
        for name in ['Aggregation_LakeSoilC','Aggregation_MethanePH']:
            previous=(args.baseline/'mksrfdata'/f'{name}.F90').read_text()
            assert aggregate_body(current,name)==aggregate_body(previous,name),name+' body changed'
        print('PASS both complete aggregation procedure bodies unchanged')
    with tempfile.TemporaryDirectory(prefix='colm-batch-surface-') as tmp:
        current=run(Path(tmp)/'current',compiler)
        if args.baseline:
            previous=run(Path(tmp)/'baseline',compiler,args.baseline)
            assert current==previous,'methane clip callbacks, paths, masks or errors changed'
            print('PASS frozen methane clipping callback bits:',len(current),'bytes')


if __name__=='__main__':
    main()
