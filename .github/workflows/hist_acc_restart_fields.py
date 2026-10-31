#!/usr/bin/env python3
"""Keep the field list of main/MOD_HistAccRestart.F90 in step with
main/MOD_Vars_1DAccFluxes.F90.

MOD_HistAccRestart saves every history accumulator allocated in
allocate_acc_fluxes so that a restart in the middle of a history window
continues the same window. Its field list is written out by hand (no include
files), between the two marker lines below. This script

  --check   (default) fails if the list differs from what allocate_acc_fluxes
            allocates, e.g. after a new a_* accumulator was added upstream;
  --write   rewrites the list in place.

Every allocated array is listed except the per-step work arrays in SKIP.
The preprocessor conditions around each allocate statement are copied, so a
field is compiled exactly when it is allocated.
"""
import re
import sys

ACC = 'main/MOD_Vars_1DAccFluxes.F90'
MOD = 'main/MOD_HistAccRestart.F90'
BEGIN = '      ! >>> fields from allocate_acc_fluxes'
END = '      ! <<< fields from allocate_acc_fluxes'

# filter_dt: per-step daylight mask; decomp_vr_tmp: per-step work array.
SKIP = {'filter_dt', 'decomp_vr_tmp'}

ALLOC = re.compile(r'^\s*allocate\s*\(\s*(\w+)\s*\((.*)\)\s*\)\s*(!.*)?$', re.I)
PP = re.compile(r'^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b', re.I)


def fields():
    src = open(ACC).read().splitlines()
    i0 = next(i for i, l in enumerate(src) if re.match(r'\s*SUBROUTINE allocate_acc_fluxes', l, re.I))
    i1 = next(i for i, l in enumerate(src) if re.match(r'\s*END SUBROUTINE allocate_acc_fluxes', l, re.I))
    out = []
    for ln in src[i0 + 1:i1]:
        if PP.match(ln):
            out.append(('pp', ln.strip()))
            continue
        if re.match(r'\s*CALL\s+mt_allocate_frcsat_count\s*\(', ln, re.I):
            out.append(('f', 'nac_frcsat', 1, 'landpatch', 'patch'))
            continue
        m = ALLOC.match(ln)
        if not m or m.group(1) in SKIP:
            continue
        dims = [d.strip() for d in m.group(2).split(',')]
        last = dims[-1].lower()
        if last == 'numpatch':
            pset, dname = 'landpatch', 'patch'
        elif last == 'numurban':
            pset, dname = 'landurban', 'urban'
        else:
            sys.exit(f'{ACC}: {m.group(1)} has an unknown vector dimension {dims[-1]}')
        if len(dims) > 3:
            sys.exit(f'{ACC}: {m.group(1)} has rank {len(dims)}; only ranks 1-3 are supported')
        out.append(('f', m.group(1), len(dims), pset, dname))
    return prune(out)


def prune(items):
    """Drop #if...#endif groups that contain no field, so the list stays readable."""
    changed = True
    while changed:
        changed = False
        for i in range(len(items)):
            if items[i][0] != 'pp' or not re.match(r'#\s*if', items[i][1]):
                continue
            depth, j = 0, i
            while j < len(items):
                if items[j][0] == 'pp':
                    if re.match(r'#\s*if', items[j][1]):
                        depth += 1
                    elif re.match(r'#\s*endif', items[j][1]):
                        depth -= 1
                        if depth == 0:
                            break
                j += 1
            if not any(t[0] == 'f' for t in items[i:j + 1]):
                del items[i:j + 1]
                changed = True
                break
    return items


def render(items):
    lines = []
    for t in items:
        if t[0] == 'pp':
            lines.append(t[1])
        else:
            _, name, rank, pset, dname = t
            if name == 'nac_frcsat':
                # Older partial windows have no trustworthy valid-sample count.
                lines.append('      CALL mt_history_frcsat_count(file, mode, HACC_WRITE, nac, a_frcsat, &')
                lines.append('         hacc_field1d, write_marker, read_marker, require_var, global_max)')
            else:
                lines.append(f"      CALL hacc_field{rank}d (file, mode, '{name}', {name}, {pset}, '{dname}')")
    return lines


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else '--check'
    text = open(MOD).read().splitlines()
    try:
        b = text.index(BEGIN)
        e = text.index(END)
    except ValueError:
        sys.exit(f'{MOD}: marker lines not found')
    new = render(fields())
    if mode == '--write':
        open(MOD, 'w').write('\n'.join(text[:b + 1] + new + text[e:]) + '\n')
        print(f'{MOD}: {sum(1 for l in new if "CALL" in l)} fields written')
    elif text[b + 1:e] != new:
        old = set(l for l in text[b + 1:e] if 'CALL' in l)
        cur = set(l for l in new if 'CALL' in l)
        for l in sorted(cur - old):
            print('missing in MOD_HistAccRestart:', l.strip())
        for l in sorted(old - cur):
            print('not allocated in allocate_acc_fluxes:', l.strip())
        if old == cur:
            print('same fields, but the order or the preprocessor conditions differ')
        print(f'Run: python3 {sys.argv[0]} --write')
        sys.exit(1)
    else:
        print(f'{MOD}: {len(cur := [l for l in new if "CALL" in l])} fields match {ACC}')


if __name__ == '__main__':
    main()
