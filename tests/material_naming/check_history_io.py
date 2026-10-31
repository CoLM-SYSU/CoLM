#!/usr/bin/env python3
"""Extract production history routines; exercise real NetCDF with MPI/empty ranks.
Adapters replace only CoLM block distribution, not field selection or validation.
"""
import pathlib, re, shlex, subprocess, tempfile
ROOT = pathlib.Path(__file__).resolve().parents[2]
def extract(text, name):
    return re.search(r'^\s*(?:logical\s+FUNCTION|SUBROUTINE|FUNCTION) '+name+r'\s*\(.*?^\s*END (?:SUBROUTINE|FUNCTION) '+name+r'\b', text, re.M|re.S|re.I).group()
src = (ROOT/'main/MOD_HistAccRestart.F90').read_text()
vector = (ROOT/'share/MOD_NetCDFVector.F90').read_text()
pre = '''
module MOD_Precision
 integer,parameter::r8=kind(1d0)
end module
module MOD_SPMD_Task
 use mpi
 integer::p_comm_glb=MPI_COMM_WORLD,p_err,myrank,nproc
 logical::p_is_io
end module
module MOD_Namelist
 integer::DEF_REST_CompressLevel=0
end module
module MOD_Pixelset
 type gs
 integer::vlen(2,1)=1
 end type
 type pixelset_type
 integer::nblkgrp=0
 integer::xblkgrp(2)=0,yblkgrp(2)=1
 type(gs)::vecgs
 end type
end module
module MOD_Block
contains
 subroutine get_filename_block(base,i,j,file)
 character(*)::base,file
 integer::i,j
 write(file,'(A,I0,A)') trim(base),i,'.nc'
 end subroutine
end module
module MOD_NetCDFSerial
 use netcdf
contains
 logical function ncio_var_exist(file,name,readflag)
 character(*)::file,name
 logical::readflag
 integer::nc,id,s
 ncio_var_exist=.false.
 s=nf90_open(file,NF90_NOWRITE,nc)
 if(s/=0)return
 ncio_var_exist=nf90_inq_varid(nc,name,id)==0
 s=nf90_close(nc)
 end function
end module
module MOD_NetCDFVector
 use netcdf
 use MOD_Precision
 use MOD_Pixelset
 use MOD_Block
 interface ncio_write_vector
 module procedure write2,write3
 end interface
 interface ncio_read_vector
 module procedure read2,read3
 end interface
contains
 subroutine ck(s)
 integer::s
 if(s/=0)then
 print *,nf90_strerror(s)
 error stop 2
 endif
 end subroutine
 subroutine ncio_define_dimension_vector(file,ps,name,n)
 character(*)::file,name
 type(pixelset_type)::ps
 integer::n,b,nc,id
 character(256)::path
 logical::exists
 do b=1,ps%nblkgrp
 call get_filename_block(file,ps%xblkgrp(b),1,path)
 inquire(file=trim(path),exist=exists)
 if(exists)then
 call ck(nf90_open(trim(path),NF90_WRITE,nc))
 call ck(nf90_redef(nc))
 else
 call ck(nf90_create(trim(path),NF90_CLOBBER,nc))
 call ck(nf90_def_dim(nc,'patch',1,id))
 endif
 call ck(nf90_def_dim(nc,name,n,id))
 call ck(nf90_close(nc))
 enddo
 end subroutine
'''
for rank in (2,3):
    dimargs = 'd1,n1,' + ('d2,n2,' if rank==3 else '')
    dims = 'n1,' + ('n2,' if rank==3 else '')
    shape = ':,'* (rank-1)+':'
    sl = ':,'*(rank-1)+'b:b'
    pre += f'''subroutine write{rank}(file,name,{dimargs}vec,ps,a,level)
 character(*)::file,name,d1,vec{',d2' if rank==3 else ''}
 integer::n1,level{',n2' if rank==3 else ''}
 type(pixelset_type)::ps
 real(r8)::a({shape})
 integer::b,nc,id,ids({rank})
 character(256)::path
 do b=1,ps%nblkgrp
 call get_filename_block(file,ps%xblkgrp(b),1,path)
 call ck(nf90_open(trim(path),NF90_WRITE,nc))
 call ck(nf90_redef(nc))
 call ck(nf90_inq_dimid(nc,d1,ids(1)))
 ''' + ('call ck(nf90_inq_dimid(nc,d2,ids(2)))\n' if rank==3 else '') + f'''
 call ck(nf90_inq_dimid(nc,vec,ids({rank})))
 call ck(nf90_def_var(nc,name,NF90_DOUBLE,ids,id))
 call ck(nf90_enddef(nc))
 call ck(nf90_put_var(nc,id,a({sl})))
 call ck(nf90_close(nc))
 enddo
 end subroutine
 subroutine read{rank}(file,name,{dims}ps,a)
 character(*)::file,name
 integer::n1{',n2' if rank==3 else ''}
 type(pixelset_type)::ps
 real(r8),allocatable::a({shape})
 integer::b,nc,id
 character(256)::path
 do b=1,ps%nblkgrp
 call get_filename_block(file,ps%xblkgrp(b),1,path)
 call ck(nf90_open(trim(path),NF90_NOWRITE,nc))
 call ck(nf90_inq_varid(nc,name,id))
 call ck(nf90_get_var(nc,id,a({sl})))
 call ck(nf90_close(nc))
 enddo
 end subroutine
'''
pre += extract(vector,'ncio_vector_exist')+'\nend module\nmodule tested\nuse MOD_Precision\nuse MOD_SPMD_Task\ninteger,parameter::HACC_WRITE=1,HACC_READ=2\ncontains\n'
for name in ('hacc_field2d','hacc_field3d','require_var','global_max'):
    pre += extract(src,name)+'\n'
history = (ROOT/'main/MATERIAL/MOD_MT_HistoryRestart.F90').read_text()
pre += extract(history,'mt_history_read_name').replace('mt_history_read_name','hacc_read_name')+'\n'
pre += '''subroutine CoLM_stop()
 call MPI_ABORT(p_comm_glb,77,p_err)
 end subroutine
end module
program check
 use tested
 use MOD_Pixelset
 use MOD_NetCDFVector,only:ck
 use netcdf
 implicit none
 type(pixelset_type)::ps,selected
 real(r8),allocatable::a(:,:),b(:,:,:)
 character(32)::mode,name,first
 character(256)::path
 integer::i,nc,id
 call MPI_INIT(p_err)
 call MPI_COMM_RANK(p_comm_glb,myrank,p_err)
 call MPI_COMM_SIZE(p_comm_glb,nproc,p_err)
 call get_command_argument(1,mode)
 call get_command_argument(2,first)
 p_is_io=nproc==1.or.myrank==1.or.myrank==2
 if(nproc==1)then
 ps%nblkgrp=2
 ps%xblkgrp=[1,2]
 elseif(p_is_io)then
 ps%nblkgrp=1
 ps%xblkgrp(1)=myrank
 endif
 allocate(a(2,ps%nblkgrp),b(2,3,ps%nblkgrp))
 a=5d0;b=7d0
 name='a_mt_test'
 if(mode/='new')name='a_trc_test'
 if(mode/='missing')then
 call hacc_field2d('two',HACC_WRITE,trim(name),a,ps,'patch')
 call hacc_field3d('three',HACC_WRITE,trim(name),b,ps,'patch')
 endif
 if(mode=='prefer'.or.mode=='mixed'.or.mode=='bad_dim'.or.mode=='bad_dim2')then
 selected=ps
 if(mode=='mixed')then
 selected%nblkgrp=0
 do i=1,ps%nblkgrp
 if(ps%xblkgrp(i)==1)then
 selected%nblkgrp=1;selected%xblkgrp(1)=1
 endif
 enddo
 endif
 a=11d0;b=13d0
 call hacc_field2d('two',HACC_WRITE,'a_mt_test',a,selected,'patch')
 call hacc_field3d('three',HACC_WRITE,'a_mt_test',b,selected,'patch')
 if(mode=='bad_dim'.or.mode=='bad_dim2')then
 do i=1,ps%nblkgrp
 if(mode=='bad_dim2')then
 write(path,'(A,I0,A)')'two',ps%xblkgrp(i),'.nc'
 else
 write(path,'(A,I0,A)')'three',ps%xblkgrp(i),'.nc'
 endif
 call ck(nf90_open(trim(path),NF90_WRITE,nc))
 call ck(nf90_redef(nc))
 if(mode=='bad_dim2')then
 call ck(nf90_inq_dimid(nc,'d1_a_mt_test',id))
 else
 call ck(nf90_inq_dimid(nc,'d2_a_mt_test',id))
 endif
 call ck(nf90_rename_dim(nc,id,'wrong_dimension'))
 call ck(nf90_close(nc))
 enddo
 endif
 endif
 call MPI_BARRIER(p_comm_glb,p_err)
 a=0;b=0
 if(first=='3')then
 call hacc_field3d('three',HACC_READ,'a_mt_test',b,ps,'patch')
 call hacc_field2d('two',HACC_READ,'a_mt_test',a,ps,'patch')
 else
 call hacc_field2d('two',HACC_READ,'a_mt_test',a,ps,'patch')
 call hacc_field3d('three',HACC_READ,'a_mt_test',b,ps,'patch')
 endif
 if(mode=='prefer')then
 if(any(a/=11d0).or.any(b/=13d0))call CoLM_stop()
 else
 if(any(a/=5d0).or.any(b/=7d0))call CoLM_stop()
 endif
 if(myrank==0)print *,'PASS ',trim(mode)
 call MPI_FINALIZE(p_err)
end program
'''
with tempfile.TemporaryDirectory(prefix='material-history-') as temp:
    root=pathlib.Path(temp)
    (root/'test.F90').write_text(pre)
    flags=shlex.split(subprocess.check_output(['nf-config','--fflags'],text=True))
    libs=shlex.split(subprocess.check_output(['nf-config','--flibs'],text=True))
    libs += shlex.split(subprocess.check_output(['nc-config','--libs'],text=True))
    subprocess.run(['mpif90','-cpp','-DUSEMPI','-ffree-line-length-none','-fcheck=all',*flags,'test.F90',*libs,'-o','check'],cwd=root,check=True)
    for ranks in (1,4):
        for mode,first in [('new','2'),('legacy','2'),('prefer','2'),('missing','2'),('mixed','2'),('bad_dim2','2'),('missing','3'),('mixed','3'),('bad_dim','3')]:
            case=root/f'{ranks}-{mode}-{first}';case.mkdir()
            p=subprocess.run(['mpiexec','-n',str(ranks),str(root/'check'),mode,first],cwd=case,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=25)
            good=mode in ('new','legacy','prefer')
            assert (p.returncode==0)==good, (ranks,mode,p.returncode,p.stdout)
            print(f'PASS ranks={ranks} case={mode} first={first}d '+('roundtrip' if good else 'rejected'))
