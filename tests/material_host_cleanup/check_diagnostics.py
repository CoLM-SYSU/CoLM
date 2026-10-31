#!/usr/bin/env python3
"""Check diagnostic object lifetime and bounds using the production module."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
compiler = shutil.which('gfortran')
if not compiler:
    raise SystemExit('gfortran is required')
with tempfile.TemporaryDirectory(prefix='colm-diagnostics-') as tmp:
    path = Path(tmp)
    (path/'define.h').write_text('#define MaterialTransport\n')
    subprocess.run([compiler, '-O0', '-g', '-cpp', '-ffree-line-length-none', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-I'+tmp,
                    str(ROOT/'share/MOD_Precision.F90'),
                    str(ROOT/'main/MATERIAL/MOD_MT_HydroDiagnostics.F90'),
                    str(Path(__file__).with_name('diagnostics.F90')), '-o', str(path/'check')],
                   cwd=path, check=True)
    subprocess.run([str(path/'check')], cwd=path, check=True, timeout=10)
