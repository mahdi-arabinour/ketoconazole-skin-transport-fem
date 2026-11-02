"""Locked local-identifiability diagnostics; optional current manuscript/SI audit."""
from pathlib import Path
import numpy as np
import pandas as pd
from _external_documents import get_external_documents
root=Path(__file__).resolve().parents[2]
x=root/'outputs/step03'
svd=pd.read_csv(x/'local_svd_summary.csv')
assert len(svd)==4 and set(svd.model_case)=={'Case 1','Case 2'}
assert (svd.n_obs==7).all()
assert set(svd[svd.model_case=='Case 1'].n_parameters)=={6}
assert set(svd[svd.model_case=='Case 2'].n_parameters)=={7}
assert (svd.condition_number>1e8).all()
assert (svd.rank_rel_1e_minus_3<svd.n_parameters).all()
assert np.allclose(svd.condition_number.to_numpy(),svd.s1.to_numpy()/svd.s_min.to_numpy(),rtol=1e-8)
pert=pd.read_csv(x/'fixed_other_parameter_perturbations.csv')
assert all((pert.factor>0)&(pert.new_objective>0))
assert np.max(abs(pert.new_objective-pert.baseline_objective-pert.delta_objective))<1e-10
assert len(pd.read_csv(x/'parameter_local_sensitivities.csv'))==26
p=pd.read_csv(x/'P23_profiled_objective.csv')
assert not p.optimizer_converged.any(), 'Unexpected profile convergence status: review before publication'
assert len(pd.read_csv(root/'outputs/step02/source_reconciled_all_cases_parameters.csv'))==4
assert len(pd.read_csv(root/'data/provenance/STEP02_18_DISCREPANCY_AUDIT.csv'))==18
print('PASS: four SVD diagnostics, 26 local-parameter sensitivities, receptor-transfer perturbations, and exploratory nonconvergence flags.')

main,si=get_external_documents()
if main is not None:
    text=' '.join(x.text for x in main.paragraphs)
    for phrase in ['Practical identifiability and case comparison','severely ill-conditioned','AICc','not uniquely']:
        assert phrase in text,phrase
    assert len(si.tables)==20,len(si.tables)
    for phrase in ['Table S17.','Table S18.','Figure S8.']:
        assert any(x.text.startswith(phrase) for x in si.paragraphs),phrase
    assert len(si.inline_shapes)>=8
    print('PASS: external manuscript discusses ill-conditioning and external SI includes sensitivity diagnostics S17/S18/S8.')
