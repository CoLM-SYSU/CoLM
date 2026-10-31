#!/usr/bin/env python3
"""Reject old material identifiers and legacy names in new history writes."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
allowed = {'etrc', 'etrc_p', 'strcyc'}
# Upstream BGC concentration terms are unrelated to the material naming migration.
upstream_allowed = {
    'main/BGC/MOD_BGC_Soil_BiogeochemLittVertTransp.F90': {'conc_trcr_c', 'conc_trcr_n'},
}
bad = []
for directory in ['main', 'share', 'mksrfdata', 'mkinidata', 'include']:
    for p in (root / directory).rglob('*'):
        if not p.is_file() or p.suffix.lower() not in ['.f90', '.h']:
            continue
        file_allowed = allowed | upstream_allowed.get(p.relative_to(root).as_posix(), set())
        for n, line in enumerate(p.read_text().splitlines(), 1):
            # The old prefix remains only as a read-only compatibility string.
            code = re.sub(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"", '', line).split('!', 1)[0]
            tokens = re.findall(r'\b[A-Za-z][A-Za-z0-9_]*\b', code)
            bad += [(str(p.relative_to(root)), n, t) for t in tokens
                    if ('trc' in t.lower() or t.lower() == 'ntr') and t.lower() not in file_allowed]
assert not bad, bad
s = (root / 'main/MATERIAL/MOD_MT_HistoryRestart.F90').read_text()
fields = re.findall(r"CALL field([23])d\s*\(file, mode, '(a_mt_[^']+)'", s)
assert len(fields) == 27 and sum(rank == '3' for rank, _ in fields) == 2
assert not re.search(r"CALL (?:hacc_)?field[123]d\s*\(file, mode, 'a_trc_", s)
print('PASS no legacy material identifiers; 27 new history fields including 2 three-dimensional fields')
