"""Numerical equation/objective checks; optional editable equation/SI S19 Word cross-check."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from _external_documents import get_external_documents
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/scripts'))
import sampling_aware_fem as forward
import fit_sampling_aware as fit
assert forward.D1==forward.D4==10 and forward.K12==forward.K23==100
assert abs(forward.A-2.2)<1e-12 and abs(forward.VR-6.5)<1e-12 and abs(forward.VS-5)<1e-12
assert fit.names==['Q0_eff_ug_cm2','D2_SC_cm2_h','D3_EP_cm2_h','P12_donor_SC','P23_SC_EP','P34_EP_receptor','K34_cm_h']
obs=pd.read_csv(ROOT/'data/processed/fem_fit_data_master.csv')
pars=pd.read_csv(ROOT/'outputs/step02/source_reconciled_all_cases_parameters.csv')
for row in pars.itertuples(index=False):
    p=np.array([float(getattr(row,n)) for n in fit.names],dtype=float)
    df,_=fit.sim(p)
    lookup=df.set_index('time_h')
    current=obs[(obs.formulation==row.formulation)&(obs.use_for_fitting==True)]
    nR=current[current.compartment=='R'].sort_values('time_h')
    nSC=current[current.compartment=='SC']
    nEP=current[current.compartment=='EP']
    assert len(nR)==5 and len(nSC)==1 and len(nEP)==1
    rr=np.array([(lookup.loc[float(t),'Q_receptor_cumulative_ug_cm2']-float(y))/float(nR.observed_ug_cm2.max())
         for t,y in zip(nR.time_h,nR.observed_ug_cm2)])
    dsc=(lookup.loc[24.0,'Q2_SC_ug_cm2']-float(nSC.observed_ug_cm2.iloc[0]))/float(nSC.observed_ug_cm2.iloc[0])
    dep=(lookup.loc[24.0,'Q3_EP_ug_cm2']-float(nEP.observed_ug_cm2.iloc[0]))/float(nEP.observed_ug_cm2.iloc[0])
    objective=float(np.dot(rr,rr)/5+dsc*dsc+dep*dep)
    assert abs(objective-row.objective)<2e-9,(row.formulation,row.model_case,objective,row.objective)
    print('PASS OBJECTIVE',row.formulation,row.model_case,'value:',objective)
print('PASS: all four objectives match the five-receptor/two-endpoint weighted objective; numerical fixed coefficients unchanged.')

man,si=get_external_documents()
if man is not None:
    def math_text(k):
        elems=man.tables[k].rows[0].cells[0]._tc.xpath('.//m:oMath')
        assert elems,('missing editable equation',k)
        return ' '.join(''.join(n.text or '' for n in el.iter() if n.tag.endswith('}t')) for el in elems)
    assert len(man.tables)>=16 and len(si.tables)==20
    for k in range(1,13):assert math_text(k),k
    e3=math_text(3)
    assert e3.count('∂x')>=2 and '∂t' not in e3,e3
    e8,e9=math_text(10),math_text(11)
    assert all(s in e8 for s in ('Q0,eff','D2','D3','P12','P34','K34')) and 'P23' not in e8
    assert all(s in e9 for s in ('Q0,eff','D2','D3','P12','P23','P34','K34'))
    e10=math_text(12)
    assert 'Φ' in e10 and 'ΔQR' in e10 and 'ΔQSC' in e10 and 'ΔQEP' in e10 and '24' in e10
    assert si.tables[19].rows[0].cells[0].text=='Equation'
    crosswalk=' '.join(c.text for row in si.tables[19].rows for c in row.cells)
    for target in ['(1)','(2)','(3)','(4)-(5)','(6)-(7)','(8)','(9)','(10)','1/sqrt(5)','100 cm/h']:
        assert target in crosswalk,target
    print('PASS: editable equations in external manuscript and crosswalk Table S19 agree with solver assumptions.')
