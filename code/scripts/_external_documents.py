"""Optional checks against separately supplied Word files (not archived with the code)."""
from pathlib import Path
import argparse

def get_external_documents():
    p=argparse.ArgumentParser(description='Numerical tests; optional separate manuscript/SI cross-checks')
    p.add_argument('--manuscript', type=Path, help='Full path to current external manuscript .docx')
    p.add_argument('--supplementary', type=Path, help='Full path to current external Supplementary Data .docx')
    a=p.parse_args()
    if bool(a.manuscript) != bool(a.supplementary):
        p.error('Supply BOTH --manuscript and --supplementary, or neither')
    if not a.manuscript:
        print('INFO: External DOCX cross-check NOT RUN (manuscript/SI intentionally excluded from archive).')
        return None, None
    for path in (a.manuscript,a.supplementary):
        if not path.is_file() or path.suffix.lower() != '.docx':
            p.error(f'External DOCX not found: {path}')
    try:
        from docx import Document
    except ImportError as e:
        p.error('Optional document checks require pip install python-docx==1.2.0')
    return Document(a.manuscript), Document(a.supplementary)
