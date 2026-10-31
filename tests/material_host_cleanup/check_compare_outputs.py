#!/usr/bin/env python3
"""Exercise strict comparison failure exits with small, independent fixtures."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset

from compare_outputs import compare


def write(path, *, dimension='time', unlimited=False, dtype='f8', values=(1., 2.),
          extra=False, group=False):
    with Dataset(path, 'w') as nc:
        nc.createDimension(dimension, None if unlimited else len(values))
        var = nc.createVariable('value', dtype, (dimension,), fill_value=-9999)
        var[:] = values
        if extra:
            nc.createVariable('extra', 'i4')[:] = 1
        if group:
            nc.createGroup('nested').createVariable('nested_value', 'f8')[:] = 3.


checks = 0
with tempfile.TemporaryDirectory(prefix='colm-compare-') as tmp:
    left, right = Path(tmp) / 'left', Path(tmp) / 'right'
    left.mkdir()
    right.mkdir()
    script = Path(__file__).with_name('compare_outputs.py')

    def check(success):
        global checks
        result = subprocess.run([sys.executable, str(script), str(left), str(right)],
                                capture_output=True, text=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        checks += 1

    check(False)  # empty directories cannot prove equivalence
    write(left / 'data.nc')
    check(False)
    write(right / 'data.nc')
    check(True)
    for change in [dict(extra=True), dict(dimension='level'), dict(unlimited=True),
                   dict(dtype='f4'), dict(values=(1.,)), dict(values=(1., 3.)),
                   dict(group=True)]:
        write(right / 'data.nc', **change)
        check(False)
    write(left / 'data.nc', group=True)
    write(right / 'data.nc', group=True)
    check(True)
    with Dataset(right / 'data.nc', 'a') as nc:
        nc.groups['nested'].variables['nested_value'][:] = 4.
    check(False)
    # Equal NaNs pass, while payload changes and signed zero remain detectable.
    nan1 = np.array([0x7ff8000000000001], dtype='u8').view('f8')[0]
    nan2 = np.array([0x7ff8000000000002], dtype='u8').view('f8')[0]
    write(left / 'data.nc', values=(nan1, 0.))
    write(right / 'data.nc', values=(nan1, 0.))
    check(True)
    write(right / 'data.nc', values=(nan2, 0.))
    check(False)
    write(right / 'data.nc', values=(nan1, -0.))
    check(False)
    # Different raw fill values are not hidden by automatic masking/scaling.
    for path, fill in [(left / 'data.nc', -9999), (right / 'data.nc', -8888)]:
        with Dataset(path, 'w') as nc:
            nc.createDimension('time', 1)
            var = nc.createVariable('value', 'i2', ('time',), fill_value=fill)
            var.scale_factor = .5
            var.set_auto_maskandscale(False)
            var[:] = fill
    check(False)
    shutil.copyfile(left / 'data.nc', right / 'data.nc')
    check(True)
    # Equal raw data still fail when masking, packing or time interpretation differs.
    for attr, value in [('missing_value', 1.), ('scale_factor', 2.), ('add_offset', 3.),
                        ('units', 'days since 2005-06-01'), ('calendar', 'noleap')]:
        write(left / 'data.nc')
        write(right / 'data.nc')
        with Dataset(right / 'data.nc', 'a') as nc:
            nc.variables['value'].setncattr(attr, value)
        check(False)
    shutil.copyfile(left / 'data.nc', right / 'data.nc')
    # Reporting truncation must not stop counting differences.
    for index in range(25):
        write(left / f'extra{index}.nc')
    differences, variables = compare(left, right)
    assert differences == 25 and variables == 1, (differences, variables)
    checks += 1

print(f'PASS {checks} strict NetCDF comparison checks')
