#!/usr/bin/env python3
"""Keep the upstream BGC exception local to its file and two identifiers."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BGC = 'main/BGC/MOD_BGC_Soil_BiogeochemLittVertTransp.F90'

with tempfile.TemporaryDirectory(prefix='colm-naming-') as tmp:
    root = Path(tmp)
    script = root / 'tests/material_naming/check_inventory.py'
    script.parent.mkdir(parents=True)
    shutil.copyfile(Path(__file__).with_name('check_inventory.py'), script)
    (root / 'main/BGC').mkdir(parents=True)
    shutil.copyfile(ROOT / 'main/MOD_HistAccRestart.F90', root / 'main/MOD_HistAccRestart.F90')
    (root / 'main/MATERIAL').mkdir(parents=True)
    shutil.copyfile(ROOT / 'main/MATERIAL/MOD_MT_HistoryRestart.F90', root / 'main/MATERIAL/MOD_MT_HistoryRestart.F90')
    bgc = root / BGC
    bgc.write_text('real :: conc_trcr_c, CONC_TRCR_N\n')

    def check(accepted, token=''):
        result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
        assert (result.returncode == 0) == accepted, result.stdout + result.stderr
        if not accepted:
            assert token in result.stderr, result.stderr

    check(True)
    bgc.write_text('real :: conc_trcr_c, conc_trcr_extra\n')
    check(False, 'conc_trcr_extra')
    bgc.write_text('real :: conc_trcr_c, conc_trcr_n\n')
    other = root / 'main/other.F90'
    other.write_text('real :: conc_trcr_c\n')
    check(False, 'conc_trcr_c')
    other.write_text('real :: a_trc_precip\n')
    check(False, 'a_trc_precip')
    other.write_text("! a_trc_comment\ncharacter(*) :: label = 'a_trc_read_compat'\n")
    check(True)

print('PASS upstream BGC names allowed only in their original file; legacy names still rejected')
