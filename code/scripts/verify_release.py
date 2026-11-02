"""Check v1.1.0 release packaging and SHA-256 manifest without altering scientific outputs."""
from pathlib import Path
import csv,hashlib,json
root=Path(__file__).resolve().parents[2]
meta=json.loads((root/'.zenodo.json').read_text(encoding='utf-8'))
assert meta['version']=='1.1.0'
assert 'version: 1.1.0' in (root/'CITATION.cff').read_text(encoding='utf-8')
assert (root/'LICENSE').is_file()
manifest=root/'MANIFEST.csv'
with manifest.open(encoding='utf-8-sig',newline='') as fh:
    rows=list(csv.DictReader(fh))
assert rows, 'Missing manifest entries'
registered={}
for row in rows:
    path=row['relative_path']
    assert path not in registered, f'Duplicate manifest path: {path}'
    p=root/path
    assert p.is_file(),f'Missing released file: {path}'
    actual=hashlib.sha256(p.read_bytes()).hexdigest()
    assert actual==row['sha256'],f'SHA256 mismatch: {path}'
    assert p.stat().st_size==int(row['size_bytes']),f'Size mismatch: {path}'
    registered[path]=actual
ignore_names={'__pycache__','.venv','.git','.ipynb_checkpoints','.pytest_cache'}
actual_files={p.relative_to(root).as_posix() for p in root.rglob('*')
    if p.is_file() and p!=manifest and not (set(p.relative_to(root).parts)&ignore_names)}
assert actual_files==set(registered),f'File list differs from manifest: extra={actual_files-set(registered)}, missing={set(registered)-actual_files}'
for disallowed in ('submission_drafts','review_notes','archive','outputs/step06','outputs/narrative_blue'):
    assert not (root/disallowed).exists(),f'Internal staging directory included: {disallowed}'
print(f'PASS: release metadata v1.1.0, {len(rows)} files checked against complete SHA-256 manifest; no Word manuscript/SI included.')
