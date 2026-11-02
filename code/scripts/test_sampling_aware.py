"""Self-contained sampling correction acceptance tests, relative to the repository root."""
import numpy as np
import pandas as pd
from sampling_aware_fem import ROOT,run_model,VS,VR,A,layer_w
from fit_sampling_aware import sim,names

par=pd.read_csv(ROOT/'outputs/step02/source_reconciled_all_cases_parameters.csv')
for f,c in [('Free KC','Case 1'),('KLZ-NPs','Case 2')]:
 row=par[(par.formulation==f)&(par.model_case==c)].iloc[0]
 p=np.array([row[k] for k in names],float)
 df,states=sim(p,return_states=True)
 assert len(df)==7
 assert (df.Q_receptor_cumulative_ug_cm2>=df.Q4_receptor_ug_cm2).all()
 assert np.max(np.abs(df.Q_receptor_cumulative_ug_cm2-(df.Q4_receptor_ug_cm2+df.Q_removed_ug_cm2)))<1e-11
 assert np.max(abs(df.Q_total_with_samples_ug_cm2-p[0]))/p[0]<5e-8
 assert np.all(np.diff(df.Q_removed_ug_cm2)>=-1e-10)
 assert np.max(abs(df.C_receptor_pre_sampling_ug_mL-(df.Q4_receptor_ug_cm2*A/VR)))<1e-10
 assert df.iloc[-1].Q_removed_ug_cm2>0
 # The alternate direct sparse time integrator must return the same sampled outputs.
 direct,_=run_model(p,[1,2,4,6,8,12,24],sampled=True,rtol=1e-7,atol=1e-9)
 err=np.max(np.abs(direct.Q_receptor_cumulative_ug_cm2.to_numpy()-df.Q_receptor_cumulative_ug_cm2.to_numpy()))
 assert err<3e-4,(f,err)
 # Sample removal at 24 h is not executed until after the terminal measurement.
 # The total at t=24 must include all earlier collected masses and current receptor only once.
 s=df.iloc[-1];total=s.Q1_donor_ug_cm2+s.Q2_SC_ug_cm2+s.Q3_EP_ug_cm2+s.Q4_receptor_ug_cm2+s.Q_removed_ug_cm2
 assert abs(total-p[0])/p[0]<5e-8
 print('PASS',f,c,'max difference vs independent BDF solver',err,'total rel error',abs(total-p[0])/p[0])
print('Sampling correction tests all passed.')
