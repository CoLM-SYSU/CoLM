#!/usr/bin/env python3
"""Run the declared NEC IsoGSM fixture in a private run-7920 snapshot."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def setting(text, key, value):
    text, count = re.subn(r'^\s*' + re.escape(key) + r'\s*=.*$',
                          lambda _: f'{key} = {value}', text, flags=re.M | re.I)
    assert count == 1, (key, count)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        assert setting('x = 1\nx_extra = 2\n', 'x', '3') == 'x = 3\nx_extra = 2\n'
        assert setting('DEF%day = 1\n', 'DEF%day', '2') == 'DEF%day = 2\n'
        for text in ['', 'x=1\nx=2\n']:
            try:
                setting(text, 'x', '3')
            except AssertionError:
                continue
            raise AssertionError('missing or duplicate setting accepted')
        print('PASS exact namelist replacement and missing/duplicate rejection')
        return

    root = Path.cwd()
    out = Path(os.environ['OBREF_JOB_RESULTS'])
    inputs = {key: Path(os.environ[f'OBREF_INPUT_{key}']) for key in ['CASE', 'RAW', 'RT', 'ISO']}
    template = (inputs['CASE'] / 'case.nml').read_text()
    expected = {'DEF_CASE_NAME': "'case'", 'DEF_simulation_time%timestep': '1800.',
                'DEF_simulation_time%spinup_repeat': '0'}
    for phase in ['start', 'end', 'spinup']:
        for field, value in [('year', '2005'), ('month', '6'), ('sec', '0')]:
            expected[f'DEF_simulation_time%{phase}_{field}'] = value
    for key, value in expected.items():
        values = re.findall(r'^\s*' + re.escape(key) + r'\s*=\s*(.*?)\s*$',
                            template, flags=re.M | re.I)
        assert values == [value], (key, values, value)
    if args.baseline:
        manifest = json.loads((root / 'validation-baseline/manifest.json').read_text())
        for name, digest in manifest.items():
            source = root / 'validation-baseline' / name
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest, name
            shutil.copy2(source, root / name)
    os.environ['PATH'] = '/home/zhwei/miniconda3/bin:' + os.environ['PATH']
    os.environ['OMP_NUM_THREADS'] = '1'

    def command(argv, cwd, log, steps=None):
        with log.open('w') as stream:
            subprocess.run(argv, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=1500)
        if steps is not None:
            text = log.read_text(errors='replace')
            assert text.count('TIMESTEP') == steps, (log, text.count('TIMESTEP'), steps)
            assert 'CoLM Execution Completed' in text, log

    inc = subprocess.check_output(['nf-config', '--includedir'], text=True).strip()
    lib = subprocess.check_output(['nc-config', '--libdir'], text=True).strip()
    # NetCDF/HDF5 may probe floating-point values internally; keep bounds checks.
    flags = ('-fdefault-real-8 -fdefault-double-8 -ffree-form -O0 -g -cpp '
             '-ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium -fcheck=bounds')
    make = ['make', 'colm.x', 'mkinidata.x', 'mksrfdata.x', 'FF=mpif90', f'NETCDF_INC={inc}',
            f'NETCDF_LIB={lib}', 'FOPTS=' + flags]
    if 'OBREF_INPUT_BUILD' in os.environ:
        cached = Path(os.environ['OBREF_INPUT_BUILD'])
        shutil.copy2(cached / 'include/define.h', root / 'include/define.h')
        changed = []
        for directory in ['main', 'share', 'include', 'mkinidata', 'mksrfdata']:
            for source in (root / directory).rglob('*'):
                if source.is_file() and source.suffix.lower() in {'.f90', '.f', '.h'}:
                    previous = cached / source.relative_to(root)
                    if previous.is_file() and source.read_bytes() == previous.read_bytes():
                        stamp = previous.stat().st_mtime_ns
                        os.utime(source, ns=(stamp, stamp))
                    else:
                        changed.append(str(source.relative_to(root)))
                        os.utime(source, None)
        shutil.copytree(cached / '.bld', root / '.bld')
        # Entry points must honor these flags even if the cache enabled FPE traps.
        for source in ['main/CoLM.F90', 'mkinidata/CoLMINI.F90']:
            os.utime(root / source, None)
        (out / 'build-source-changes.json').write_text(json.dumps(changed, indent=2) + '\n')
        command(make, root, out / 'build.log')
    else:
        command(['bash', '.github/workflows/create_defineh.bash', 'GRID', 'LULC_IGBP_PFT',
                 'URBANOFF', 'vanGenu', 'CaMaOFF', 'BGCON', 'CROPON', 'MTON', 'LULCCOFF'],
                root, out / 'define.log')
        command(make, root, out / 'build.log')
    shutil.copyfile(root / 'include/define.h', out / 'define.h')
    print('BUILD_PASS', flush=True)
    (out / 'binaries.sha256').write_text('\n'.join(
        f'{hashlib.sha256((root / "run" / name).read_bytes()).hexdigest()}  {name}'
        for name in ['colm.x', 'mkinidata.x', 'mksrfdata.x']) + '\n')

    isotope = (root / 'run/standard_isotope_parameter.nml').read_text()
    isotope = setting(isotope, 'DEF_ISO_USE_FRACTIONATION', '.true.')
    (out / 'isotope.nml').write_text(isotope)

    def prepare(name, end_day):
        case = out / 'run' / name
        (case / 'case').mkdir(parents=True)
        shutil.copytree(inputs['CASE'] / 'case/landdata', case / 'case/landdata')
        forcing = setting((root / 'run/forcing/IsoGSM.nml').read_text(),
                          'DEF_dir_forcing', f"'{inputs['ISO']}/'")
        (case / 'forcing.nml').write_text(forcing)
        fields = {
            'DEF_dir_rawdata': f"'{inputs['RAW']}/'", 'DEF_dir_runtime': f"'{inputs['RT']}/'",
            'DEF_dir_output': f"'{case}/'", 'DEF_rawdata_namelist': f"'{root}/run/rawdata/colm2024.nml'",
            'DEF_forcing_namelist': f"'{case}/forcing.nml'",
            'DEF_UnitCatchment_file': f"'{inputs['RT']}/unitcatchment/grid_routing_data_15min.nc'",
            'DEF_ReservoirPara_file': f"'{inputs['RT']}/unitcatchment/grid_routing_data_15min.nc'",
            'DEF_ISOTOPE_PARAM_FILE': f"'{out}/isotope.nml'",
            'DEF_USE_ISOTOPE': '.true.', 'DEF_USE_SEDIMENT': '.true.',
            'DEF_SOLUTE_NAMES': "'CL'",
            'DEF_SOLUTE_PARAM_FILES': f"'CL:{root}/run/standard_chloride_parameter.nml'",
            'DEF_simulation_time%start_day': '1', 'DEF_simulation_time%end_day': str(end_day),
        }
        text = template
        for key, value in fields.items():
            text = setting(text, key, value)
        (case / 'case.nml').write_text(text)
        command(['mpiexec', '-n', '16', str(root / 'run/mkinidata.x'), 'case.nml'],
                case, case / 'init.log')
        return case

    continuous = prepare('continuous', 3)
    command(['mpiexec', '-n', '16', str(root / 'run/colm.x'), 'case.nml'],
            continuous, continuous / 'colm.log', steps=96)
    print('CONTINUOUS_PASS 96 steps', flush=True)
    split = prepare('split', 2)
    command(['mpiexec', '-n', '16', str(root / 'run/colm.x'), 'case.nml'],
            split, split / 'colm.first.log', steps=48)
    nml = split / 'case.nml'
    nml.write_text(setting(setting(nml.read_text(), 'DEF_simulation_time%start_day', '2'),
                           'DEF_simulation_time%end_day', '3'))
    command(['mpiexec', '-n', '16', str(root / 'run/colm.x'), 'case.nml'],
            split, split / 'colm.second.log', steps=48)
    print('SPLIT_PASS 48+48 steps', flush=True)
    (out / 'result.json').write_text(json.dumps({
        'variant': 'baseline' if args.baseline else 'current', 'continuous_steps': 96,
        'split_steps': [48, 48],
    }, indent=2) + '\n')
    print('PASS all model segments completed; run verify_coupled.py for numerical comparison', flush=True)


if __name__ == '__main__':
    main()
