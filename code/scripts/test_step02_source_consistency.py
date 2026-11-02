"""Source-derived input/locked-output checks; optional current manuscript/SI audit."""
from pathlib import Path
import numpy as np
import pandas as pd
from _external_documents import get_external_documents
ROOT=Path(__file__).resolve().parents[2]
canon=pd.read_csv(ROOT/'data/provenance/STEP02_CANONICAL_PORCINE_PROFILE.csv')
audit=pd.read_csv(ROOT/'data/provenance/STEP02_18_DISCREPANCY_AUDIT.csv')
assert len(audit)==18, len(audit)
assert (audit.SI_STEP01_ug_cm2 != audit.ZIP_STEP01_ug_cm2).all()
assert audit.groupby('compartment').size().to_dict()=={'EP':13,'R':2,'SC':3}
assert len(canon)==42, canon.groupby('compartment').size()
assert not canon.duplicated(['formulation','compartment','time_h']).any()
for a in audit.itertuples():
    current=canon[(canon.formulation==a.formulation)&(canon.compartment==a.compartment)&(canon.time_h==a.time_h)]
    assert len(current)==1,(a.formulation,a.compartment,a.time_h)
    assert abs(float(current.iloc[0].canonical_ug_cm2)-a.reconciled_STEP02_ug_cm2)<1e-8
full=pd.read_csv(ROOT/'data/processed/full_profile_data_clean.csv')
for x in canon.itertuples():
    rec=full[(full.formulation==x.formulation)&(full.compartment==x.compartment)&(np.isclose(full.time_h,x.time_h))]
    assert len(rec)==1,('missing full',x)
    assert abs(float(rec.iloc[0].observed_ug_cm2)-x.canonical_ug_cm2)<1e-8,('value different',x)
fit=pd.read_csv(ROOT/'data/processed/fem_fit_data_master.csv')
for f,g in fit.groupby('formulation'):
    assert len(g)==7
    assert set(g[g.compartment=='R'].time_h)=={4.,6.,8.,12.,24.}
    assert len(g[g.compartment.isin(['SC','EP'])])==2
    assert bool(g.use_for_fitting.all())
assert canon[(canon.formulation=='Free KC')&(canon.compartment=='R')&(canon.time_h==24)].iloc[0].canonical_ug_cm2==26.81
assert canon[(canon.formulation=='KLZ-NPs')&(canon.compartment=='R')&(canon.time_h==24)].iloc[0].canonical_ug_cm2==33.19
p=pd.read_csv(ROOT/'outputs/step02/source_reconciled_all_cases_parameters.csv')
assert len(p)==4 and p.success.all(), p[['formulation','model_case','success']]
profile=pd.read_csv(ROOT/'outputs/models/selected_final_model_dense_simulation.csv')
assert len(profile[profile.time_h==24])==2
for z in profile.itertuples():
    assert abs(z.Q_receptor_cumulative_ug_cm2-z.Q4_receptor_ug_cm2-z.Q_removed_ug_cm2)<1e-8
    amount=float(p[(p.formulation==z.formulation)&(p.model_case==z.model_case)].iloc[0].Q0_eff_ug_cm2)
    assert abs(z.Q_total_with_samples_ug_cm2-amount)<1e-4
print('PASS: 18-source-discrepancy audit, canonical/processed inputs, seven targets per formulation, four stored fits and sampling mass ledger.')

m,s=get_external_documents()
if m is not None:
    assert len(s.tables)==20, len(s.tables)
    assert len(s.tables[4].columns)==5, 'S4 should have five columns'
    for compartment,ti in [('R',2),('SC',3),('EP',4)]:
        t=s.tables[ti]
        if compartment=='R':
            assert t.cell(2,1).text=='Not resolved' and t.cell(2,2).text=='Not resolved'
        for row in t.rows[1:]:
            try:time=float(row.cells[0].text)
            except ValueError:continue
            if compartment=='R' and time<4:continue
            for col,form in [(1,'Free KC'),(2,'KLZ-NPs')]:
                actual=float(row.cells[col].text.split(' +/-')[0].split(' ±')[0].strip())
                x=canon[(canon.formulation==form)&(canon.compartment==compartment)&(np.isclose(canon.time_h,time))]
                assert len(x)==1,(compartment,form,time)
                assert abs(actual-float(x.iloc[0].canonical_ug_cm2))<1e-7,(compartment,form,time,actual)
    fulltext=' '.join(p.text for p in m.paragraphs)
    assert 'not statistically significant' in fulltext
    assert '18 transcription discrepancies' in fulltext
    assert 'approximate' in fulltext
    print('PASS: external current manuscript source-precision statements and all Supplementary S2-S4 numeric entries match locked inputs.')
