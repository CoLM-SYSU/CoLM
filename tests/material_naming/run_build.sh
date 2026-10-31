#!/usr/bin/env bash
set -euo pipefail
export PATH=/home/zhwei/miniconda3/bin:$PATH
grid=GRID
case ${1:?pc, pft, off, lulcc, single, or vector} in
  pc)    lc=LULC_IGBP_PC;  bgc=BGCON;  mt=MTON;  lulcc=LULCCOFF ;;
  pft)   lc=LULC_IGBP_PFT; bgc=BGCON;  mt=MTON;  lulcc=LULCCOFF ;;
  off)   lc=LULC_IGBP;     bgc=BGCOFF; mt=MTOFF; lulcc=LULCCOFF ;;
  lulcc) lc=LULC_IGBP;     bgc=BGCOFF; mt=MTON;  lulcc=LULCCON ;;
  single) grid=SinglePoint; lc=LULC_IGBP; bgc=BGCOFF; mt=MTON; lulcc=LULCCOFF ;;
  vector) grid=UNSTRUCTURED; lc=LULC_IGBP_PC; bgc=BGCON; mt=MTON; lulcc=LULCCOFF ;;
  *) exit 2 ;;
esac
out=${OBREF_JOB_RESULTS:?runner output directory}
bash .github/workflows/create_defineh.bash "$grid" "$lc" URBANOFF vanGenu CaMaOFF "$bgc" CROPOFF "$mt" "$lulcc" > "$out/config.log"
cp include/define.h "$out/define.h"
make colm.x mkinidata.x mksrfdata.x FF=mpif90 NETCDF_INC="$(nf-config --includedir)" NETCDF_LIB="$(nc-config --libdir)" \
  FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O0 -g -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium -fcheck=bounds -ffpe-trap=invalid,zero,overflow' \
  > "$out/build.log" 2>&1
test -x run/colm.x
test -x run/mkinidata.x
test -x run/mksrfdata.x
sha256sum run/colm.x run/mkinidata.x run/mksrfdata.x > "$out/binaries.sha256"
printf 'PASS complete %s %s %s %s %s compile and link\n' "$grid" "$lc" "$bgc" "$mt" "$lulcc" | tee "$out/result.txt"
