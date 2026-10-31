#!/usr/bin/env python3
"""Run real soil physics; optionally compare with a pre-refactor source directory.

Only model-wide configuration and MPI diagnostics are stubbed. Hydraulic
functions, root/aquifer exchange, Richards solver and location helper are real.
"""
from pathlib import Path
import argparse
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=Path, help='directory containing pre-refactor MOD_Hydro_SoilWater.F90')
args = parser.parse_args()
compiler = shutil.which('gfortran')
if not compiler:
    raise SystemExit('gfortran is required')


def run(directory, legacy=False, baseline=None):
    directory.mkdir()
    (directory / 'define.h').write_text('#define MaterialTransport\n#define vanGenuchten_Mualem_SOIL_MODEL\n')
    (directory / 'stubs.f90').write_text('''module MOD_Namelist
logical :: DEF_USE_PLANTHYDRAULICS=.false.
end module
module MOD_SPMD_Task
end module
''')
    sources = [ROOT/'share/MOD_Precision.F90', directory/'stubs.f90',
               ROOT/'main/MOD_Const_Physical.F90', ROOT/'share/MOD_UserDefFun.F90',
               ROOT/'main/HYDRO/MOD_Hydro_SoilFunction.F90']
    if not legacy:
        sources += [(baseline/'MOD_MT_HydroDiagnostics.F90') if baseline
                    else ROOT/'main/MATERIAL/MOD_MT_HydroDiagnostics.F90']
    sources += [(baseline/'MOD_Hydro_SoilWater.F90') if baseline else ROOT/'main/HYDRO/MOD_Hydro_SoilWater.F90',
                HERE/'soil_column.F90']
    flags = ['-O0', '-g', '-cpp', '-fdefault-real-8', '-fdefault-double-8',
             '-ffree-line-length-none', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-I'+str(directory)]
    if legacy:
        flags += ['-DLEGACY_INTERFACE']
    subprocess.run([compiler, *flags, *map(str, sources), '-o', 'check'], cwd=directory, check=True)
    subprocess.run([str(directory/'check')], cwd=directory, check=True, timeout=30)
    return (directory/'column.bin').read_bytes()


with tempfile.TemporaryDirectory(prefix='colm-soil-column-') as tmp:
    current = run(Path(tmp)/'current')
    if args.baseline:
        baseline = args.baseline.resolve()
        legacy = 'mt_water_diagnostics' not in (baseline/'MOD_Hydro_SoilWater.F90').read_text().lower()
        old = run(Path(tmp)/'baseline', legacy=legacy, baseline=baseline)
        assert current == old, 'host water outputs or path diagnostics changed from baseline'
        print(f'PASS bitwise baseline comparison: {len(current)//8} real outputs')
