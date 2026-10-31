#!/usr/bin/env python3
"""Strictly compare NetCDF structure and raw variable values in two directories."""
from pathlib import Path
import argparse

import numpy as np
from netCDF4 import Dataset


def same_values(left, right):
    left, right = np.asarray(left), np.asarray(right)
    if left.shape != right.shape or left.dtype != right.dtype:
        return False
    if left.dtype.hasobject:
        return all(same_values(x, y) for x, y in zip(left.flat, right.flat))
    return left.tobytes() == right.tobytes()


def groups(dataset, prefix=''):
    result = {prefix: dataset}
    for name, group in dataset.groups.items():
        result.update(groups(group, prefix + '/' + name))
    return result


def compare(left, right):
    left, right = Path(left), Path(right)
    differences = variables = 0

    def mismatch(message):
        nonlocal differences
        differences += 1
        if differences <= 20:
            print('DIFF ' + message)

    if not left.is_dir() or not right.is_dir():
        mismatch('both inputs must be existing directories')
        return differences, variables
    files_left = {p.relative_to(left) for p in left.rglob('*.nc') if p.is_file()}
    files_right = {p.relative_to(right) for p in right.rglob('*.nc') if p.is_file()}
    if not files_left or not files_right:
        mismatch('both inputs must contain NetCDF files')
    for path in sorted(files_left ^ files_right):
        mismatch(f'file only in {"left" if path in files_left else "right"}: {path}')
    for path in sorted(files_left & files_right):
        with Dataset(left / path) as a, Dataset(right / path) as b:
            a.set_auto_maskandscale(False)
            b.set_auto_maskandscale(False)
            a.set_auto_chartostring(False)
            b.set_auto_chartostring(False)
            ga, gb = groups(a), groups(b)
            for name in sorted(ga.keys() ^ gb.keys()):
                mismatch(f'{path}: group {name}')
            for group in sorted(ga.keys() & gb.keys()):
                x, y = ga[group], gb[group]
                dx = {n: (len(d), d.isunlimited()) for n, d in x.dimensions.items()}
                dy = {n: (len(d), d.isunlimited()) for n, d in y.dimensions.items()}
                if dx != dy:
                    mismatch(f'{path}:{group}: dimensions')
                for name in sorted(x.variables.keys() ^ y.variables.keys()):
                    mismatch(f'{path}:{group}/{name}: variable missing')
                for name in sorted(x.variables.keys() & y.variables.keys()):
                    u, v = x[name], y[name]
                    variables += 1
                    for attr in ['_FillValue', 'missing_value', 'scale_factor', 'add_offset',
                                 '_Unsigned', 'units', 'calendar']:
                        present = attr in u.ncattrs(), attr in v.ncattrs()
                        if present[0] != present[1] or (all(present) and
                                not same_values(u.getncattr(attr), v.getncattr(attr))):
                            mismatch(f'{path}:{group}/{name}: {attr}')
                    if (u.dtype, u.dimensions, u.shape) != (v.dtype, v.dimensions, v.shape):
                        mismatch(f'{path}:{group}/{name}: type/dimensions/shape')
                    elif not same_values(u[...], v[...]):
                        mismatch(f'{path}:{group}/{name}: raw values')
    print(f'files left={len(files_left)} right={len(files_right)}, '
          f'compared variables={variables}, differences={differences}')
    return differences, variables


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    args = parser.parse_args()
    count, _ = compare(args.left, args.right)
    raise SystemExit(1 if count else 0)
