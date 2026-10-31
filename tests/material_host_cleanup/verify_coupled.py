#!/usr/bin/env python3
"""Verify the fixed two-day NEC fixture and compare both model revisions."""
from contextlib import redirect_stdout
import json
import os
from pathlib import Path
import shutil

import numpy as np
from netCDF4 import Dataset

from compare_outputs import compare, groups


def validate(case):
    history = case / 'case/history'
    expected = {'case_hist_2005-06.nc', 'case_hist_material_2005-06.nc',
                'case_hist_unitcat_2005-06.nc'}
    assert {p.name for p in history.glob('*.nc')} == expected, history
    for name in expected:
        with Dataset(history / name) as nc:
            time = nc.variables['time']
            assert time.units == 'minutes since 1900-1-1 0:0:0', (name, time.units)
            assert np.array_equal(time[:], [55443600, 55445040]), (name, time[:])
            if 'history_window_seconds' in nc.variables:
                assert np.array_equal(nc['history_window_seconds'][:], [86400, 86400]), name
    with Dataset(history / 'case_hist_unitcat_2005-06.nc') as nc:
        for name in ['f_mt_flux_H2_18O', 'f_mt_flux_HDO', 'f_mt_flux_CL', 'f_sedout_1']:
            for day in range(2):
                values = np.ma.asarray(nc.variables[name][day]).compressed()
                assert values.size and np.isfinite(values).all(), (case, name, day)
    with Dataset(history / 'case_hist_material_2005-06.nc') as nc:
        for species in ['H2_18O', 'HDO', 'CL']:
            assert any(name.endswith('_' + species) for name in nc.variables), (case, species)
    final = case / 'case/restart/2005-154-00000'
    expected_restart = {
        'case_restart_2005-154-00000_lc2005_e125_n45.nc',
        'case_restart_bgc_2005-154-00000_lc2005_e125_n45.nc',
        'case_restart_gridriver_2005-154-00000_lc2005.nc',
        'case_restart_hist_2005-154-00000_e125_n45.nc',
        'case_restart_pft_2005-154-00000_lc2005_e125_n45.nc',
    }
    assert {p.name for p in final.glob('*.nc')} == expected_restart, final
    for directory in [history, case / 'case/restart']:
        for path in directory.rglob('*.nc'):
            with Dataset(path) as nc:
                for group, node in groups(nc).items():
                    for name, variable in node.variables.items():
                        if np.dtype(variable.dtype).kind in 'fc':
                            values = np.ma.asarray(variable[...]).compressed()
                            assert np.isfinite(values).all(), (path, group, name, 'nonfinite valid data')


def main():
    out = Path(os.environ['OBREF_JOB_RESULTS'])
    variants = {name: Path(os.environ[f'OBREF_INPUT_{name.upper()}'])
                for name in ['baseline', 'current']}
    for name, directory in variants.items():
        runner_log = Path(os.environ[f'OBREF_INPUT_{name.upper()}_LOG'])
        text = runner_log.read_text(errors='replace')
        assert 'CONTINUOUS_PASS 96 steps' in text, (name, 'continuous exit not confirmed')
        assert 'SPLIT_PASS 48+48 steps' in text, (name, 'split exits not confirmed')
        for mode in ['continuous', 'split']:
            case = directory / 'run' / mode
            logs_to_check = [('colm.log', 96)] if mode == 'continuous' else [
                ('colm.first.log', 48), ('colm.second.log', 48)]
            for filename, steps in logs_to_check:
                text = (case / filename).read_text(errors='replace')
                assert text.count('TIMESTEP') == steps, (case, filename, steps)
                assert 'CoLM Execution Completed' in text, (case, filename)
            assert 'CoLM Initialization Execution Completed' in (case / 'init.log').read_text(), case
            validate(case)
            logs = out / 'model-logs' / name / mode
            logs.mkdir(parents=True)
            for path in case.glob('*.log'):
                shutil.copyfile(path, logs / path.name)
            for filename in ['case.nml', 'forcing.nml']:
                shutil.copyfile(case / filename, logs / filename)
        shutil.copyfile(runner_log, out / 'model-logs' / name / 'runner.log')
        for filename in ['build.log', 'define.log', 'define.h', 'isotope.nml',
                         'binaries.sha256', 'build-source-changes.json', 'result.json', 'compare-history.log']:
            path = directory / filename
            if path.is_file():
                shutil.copyfile(path, out / 'model-logs' / name / filename)

    checks = []

    def check(label, left, right):
        with (out / (label + '.log')).open('w') as stream, redirect_stdout(stream):
            differences, variables = compare(left, right)
        checks.append({'name': label, 'differences': differences, 'variables': variables})
        print(label, differences, 'differences;', variables, 'variables', flush=True)

    for mode in ['continuous', 'split']:
        for category in ['history', 'restart']:
            check(f'revisions-{mode}-{category}',
                  variants['baseline'] / 'run' / mode / 'case' / category,
                  variants['current'] / 'run' / mode / 'case' / category)
    for variant, directory in variants.items():
        for category in ['history', 'restart/2005-154-00000']:
            check(f'{variant}-split-' + category.split('/')[0],
                  directory / 'run/continuous/case' / category,
                  directory / 'run/split/case' / category)
    (out / 'comparison.json').write_text(json.dumps(checks, indent=2) + '\n')
    assert all(item['differences'] == 0 for item in checks), checks
    print('PASS both revisions: full history coverage, species output and restart equivalence')


if __name__ == '__main__':
    main()
