ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"/../../
if [ $# -eq 2 ];then
   if [ "${1:0:1}" == '/' ];then
      CASENAME1=`echo "${1##*/}"`
      CASEPATH1=`echo "${1%/*}"`
      echo $CASEPATH1
      echo $CASENAME1
   else
      CASEPATHNAME1=$PWD/$1
      CASENAME1=`echo "${CASEPATHNAME1##*/}"`
      CASEPATH1=`echo "${CASEPATHNAME1%/*}"`
      echo $CASEPATH1
      echo $CASENAME1
   fi
   if [ "${2:0:1}" == '/' ];then
      CASENAME2=`echo "${2##*/}"`
      CASEPATH2=`echo "${2%/*}"`
      echo $CASEPATH2
      echo $CASENAME2
   else
      CASEPATHNAME2=$PWD/$2
      CASENAME2=`echo "${CASEPATHNAME2##*/}"`
      CASEPATH2=`echo "${CASEPATHNAME2%/*}"`
      echo $CASEPATH2
      echo $CASENAME2
   fi

   echo diff between case $1 and $2
   cd $CASEPATH1/$CASENAME1/bld/main/
   for files in *F90
   do
#      echo main/$files
      tmp=`diff $files /$CASEPATH2/$CASENAME2/bld/main/`
      if [ ! -z "$tmp" ];then
	 echo "file differs !!!!!!!!!!" main/$files
	 diff $files $CASEPATH2/$CASENAME2/bld/main/
      fi
   done

   cd BGC
   for files in *F90
   do
#      echo main/BGC/$files
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/main/BGC`
      if [ ! -z "$tmp" ];then
	 echo "file differs !!!!!!!!!!" main/$files
         diff $files $CASEPATH2/$CASENAME2/bld/main/BGC
      fi
   done

   cd ../HYDRO
   for files in *F90
   do
#      echo main/HYDRO/$files
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/main/HYDRO`
      if [ ! -z "$tmp" ];then
	 echo "file differs !!!!!!!!!!" main/$files
         diff $files $CASEPATH2/$CASENAME2/bld/main/HYDRO
      fi
   done

   cd ../URBAN
   for files in *F90
   do
#      echo main/URBAN/$files
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/main/URBAN`
      if [ ! -z "$tmp" ];then
	 echo "file differs !!!!!!!!!!" main/$files
         diff $files $CASEPATH2/$CASENAME2/bld/main/URBAN
      fi
   done

   cd ../../mksrfdata/
   for files in *F90
   do
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/mksrfdata`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!!" mksrfdata/$files
         diff $files $CASEPATH2/$CASENAME2/bld/mksrfdata 
      fi
   done

   cd ../mkinidata/
   for files in *F90
   do
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/mkinidata`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" mkinidata/$files
         diff $files $CASEPATH2/$CASENAME2/bld/mkinidata
      fi
   done

   cd ../share/
   for files in *F90
   do
      tmp=`diff $files $CASEPATH2/$CASENAME2/bld/share`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" share/$files
         diff $files $CASEPATH2/$CASENAME2/bld/share
      fi
   done
fi

if [ $# -eq 1 ];then
   if [ "${1:0:1}" == '/' ];then
      CASENAME1=`echo "${1##*/}"`
      CASEPATH1=`echo "${1%/*}"`
      echo $CASEPATH1
      echo $CASENAME1
   else
      CASEPATHNAME1=$PWD/$1
      CASENAME1=`echo "${CASEPATHNAME1##*/}"`
      CASEPATH1=`echo "${CASEPATHNAME1%/*}"`
      echo $CASEPATH1
      echo $CASENAME1
   fi
   echo compare CASE1 $CASEPATH1/$CASENAME1 with ROOT: $ROOT
	
   echo diff between case $1 and root
   cd $CASEPATH1/$CASENAME1/bld/main
   for files in *F90
   do
      tmp=`diff $files $ROOT/main/`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" main/$files
         diff $files $ROOT/main/
      fi
   done

   cd BGC
   for files in *F90
   do
      tmp=`diff $files $ROOT/main/BGC`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" main/BGC/$files
         diff $files $ROOT/main/BGC
      fi
   done

   cd ../HYDRO
   for files in *F90
   do
      tmp=`diff $files $ROOT/main/HYDRO`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" main/HYDRO/$files
         diff $files $ROOT/main/HYDRO
      fi
   done

   cd ../URBAN
   for files in *F90
   do
      tmp=`diff $files $ROOT/main/URBAN`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" main/URBAN/$files
         diff $files $ROOT/main/URBAN
      fi
   done

   cd ../../mksrfdata/
   for files in *F90
   do
      tmp=`diff  $files $ROOT/mksrfdata`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" mksrfdata/$files
         diff $files $ROOT/mksrfdata
      fi
   done

   cd ../mkinidata/
   for files in *F90
   do 
      tmp=`diff $files $ROOT/mkinidata`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" mkinidata/$files
         diff $files $ROOT/mkinidata
      fi
   done

   cd ../share/
   for files in *F90
   do
      tmp=`diff $files $ROOT/share`
      if [ ! -z "$tmp" ];then
         echo "file differs !!!!!!!!" share/$files
         diff $files $ROOT/share
      fi
   done
fi
