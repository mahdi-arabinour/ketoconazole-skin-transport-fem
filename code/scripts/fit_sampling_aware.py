import sys,time,json,argparse
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import eig,solve
from scipy.optimize import least_squares
import sampling_aware_fem as m

pdir=m.ROOT/'outputs/step02';pdir.mkdir(parents=True,exist_ok=True);root=m.ROOT
obs=pd.read_csv(root/'data/processed/fem_fit_data_master.csv');obs=obs[obs.use_for_fitting==True].copy()
_base=pd.read_csv(root/'outputs/step02/source_reconciled_all_cases_parameters.csv')
prior={1:_base[_base.model_case=='Case 1'].copy(),2:_base[_base.model_case=='Case 2'].copy()}
names=['Q0_eff_ug_cm2','D2_SC_cm2_h','D3_EP_cm2_h','P12_donor_SC','P23_SC_EP','P34_EP_receptor','K34_cm_h']
lo=np.array([1,1e-10,1e-8,1e-3,1e-3,1e-3,1e-5]); hi=np.array([5000,1e-2,1,1e4,1e4,1e4,100])
receptor_times=[4,6,8,12,24]
all_times=[1,2,4,6,8,12,24]

def sim(p,times=all_times,sampling=True,return_states=False):
 q0,d2,d3,p12,p23,p34,k34=p
 op=m.make_linear_operator(*p)
 vals,vecs=eig(op,overwrite_a=True,check_finite=False)
 if np.max(np.abs(np.imag(vals)))>1e-6 or np.max(vals.real)>1e-5:
  raise FloatingPointError('Unstable eigen spectrum')
 inv=solve(vecs,np.eye(m.n),check_finite=False,overwrite_a=False)
 y=np.zeros(m.n);y[m.reg[1]]=q0/m.length[1]
 removed=0.;t0=0.;rows=[];states={}
 events=sorted(set(map(float,times))|({z for z in m.SAMPLE_TIMES if z<=max(times)} if sampling else set())|{0.})
 for t in events:
  if t!=0.:
   y=np.real_if_close(vecs@(np.exp(vals*(t-t0))*(inv@y)),tol=1e4)
   if np.iscomplexobj(y): raise FloatingPointError('complex sim')
  if t in times:
   parts=[float(w@y) for w in (m.layer_w[1],m.layer_w[2],m.layer_w[3],m.layer_w[4])]
   rows.append({'time_h':t,'Q1_donor_ug_cm2':parts[0],'Q2_SC_ug_cm2':parts[1],
    'Q3_EP_ug_cm2':parts[2],'Q4_receptor_ug_cm2':parts[3],
    'Q_removed_ug_cm2':removed,'Q_receptor_cumulative_ug_cm2':parts[3]+removed,
    'Q_total_with_samples_ug_cm2':sum(parts)+removed,'C_receptor_pre_sampling_ug_mL':parts[3]*m.A/m.VR})
   if return_states: states[t]=y.copy()
  if sampling and t in m.SAMPLE_TIMES and t<max(times):
   removed+=float(m.layer_w[4]@y)*(m.VS/m.VR)
   y[m.reg[4]]*=1-m.VS/m.VR
  t0=t
 return pd.DataFrame(rows),states

def as7(row,case):
 p=[float(row[k]) if k!='P23_SC_EP' else (1 if case==1 else float(row[k])) for k in names]
 return np.array(p)

def fit(form,case, max_nfev=150):
 d=obs[obs.formulation==form]
 scale={'R':float(d[d.compartment=='R'].observed_ug_cm2.max()),'SC':float(d[d.compartment=='SC'].observed_ug_cm2.max()),'EP':float(d[d.compartment=='EP'].observed_ug_cm2.max())}
 data=d.sort_values(['time_h','compartment'])
 years=data.time_h.to_numpy(float)
 targets=data.observed_ug_cm2.to_numpy(float)
 compartments=data.compartment.to_numpy()
 ws=np.array([1/np.sqrt(5) if c=='R' else 1 for c in compartments]);ss=np.array([scale[c] for c in compartments]);
 mask=np.array([0,1,2,3,5,6]) if case==1 else np.arange(7)
 lower=np.log(lo[mask]);upper=np.log(hi[mask])
 nfev=0
 def resid(z):
  nonlocal nfev
  p=np.ones(7);p[mask]=np.exp(z)
  df,_=sim(p)
  indexed=df.set_index('time_h')
  pred=np.array([indexed.loc[t,('Q_receptor_cumulative_ug_cm2' if c=='R' else 'Q2_SC_ug_cm2' if c=='SC' else 'Q3_EP_ug_cm2')] for t,c in zip(years,compartments)])
  nfev+=1
  return ws*(pred-targets)/ss
 starts=[]
 orig=prior[case].loc[prior[case].formulation==form].iloc[0]
 starts.append(as7(orig,case))
 alternate=prior[3-case].loc[prior[3-case].formulation==form].iloc[0]
 a=as7(alternate,3-case)
 if case==1:a[4]=1
 starts.append(a)
 if case==2:
  for val in (0.5,2.):
   s=as7(orig,case);s[4]=val;starts.append(s)
 else:
  s=as7(orig,case);s[0]*=.80;s[2]*=1.2;starts.append(s)
 # Avoid overexploration: 3/4 starts and bound iteration caps
 fits=[]
 for ix,x0 in enumerate(starts[:2]):
  x0=np.clip(x0,lo*1.01,hi/1.01)
  st=time.perf_counter()
  try:
   res=least_squares(resid,np.log(x0[mask]),bounds=(lower,upper),max_nfev=max_nfev,method='trf',ftol=1e-7,gtol=1e-7,xtol=1e-7,x_scale='jac')
   f=float(np.sum(res.fun**2))
   print('FIT',form,case,'start',ix,'obj',f,'cost',res.nfev,'seconds',time.perf_counter()-st,'message',res.message,flush=True)
   fits.append((f,res,ix))
  except Exception as e:
   print('FAILED',form,case,ix,repr(e),flush=True)
 best=min(fits,key=lambda x:x[0]);out=np.ones(7);out[mask]=np.exp(best[1].x)
 pred,_=sim(out)
 print('BEST',form,case,best[0],np.array2string(out,precision=8),flush=True)
 return {'formulation':form,'model_case':f'Case {case}','objective':best[0],'success':bool(best[1].success),'message':best[1].message,
         'start_index':best[2],'nfev':best[1].nfev,**dict(zip(names,out))}, pred

if __name__=='__main__':
 parser=argparse.ArgumentParser(description='Sampling-aware FEM: locked STEP02 verification or explicit refit.')
 parser.add_argument('--refit',action='store_true',help='Run expensive exploratory refit and overwrite Step02 locked parameters. Review resulting manuscript values manually.')
 opt=parser.parse_args()
 if not opt.refit:
  old=pd.read_csv(pdir/'source_reconciled_all_cases_parameters.csv')
  for _,r in old.iterrows():
   result,_=sim([float(r[n]) for n in names])
   assert np.isfinite(result.Q_receptor_cumulative_ug_cm2).all()
  print('PASS: all four STEP02 locked parameter sets produce finite sampling-aware outputs; no refit performed. Use --refit to change coefficients.')
  sys.exit(0)
 st=time.perf_counter()
 rows=[];profiles=[]
 for case in (1,2):
  for form in ('Free KC','KLZ-NPs'):
   r,profile=fit(form,case)
   rows.append(r)
   profile.insert(0,'formulation',form);profile.insert(1,'model_case',f'Case {case}');profiles.append(profile)
 result=pd.DataFrame(rows);result.to_csv(pdir/'source_reconciled_all_cases_parameters.csv',index=False)
 pd.concat(profiles).to_csv(pdir/'source_reconciled_all_cases_profiles.csv',index=False)
 print('FINAL SEC',time.perf_counter()-st,flush=True)
 print(result.to_string(index=False),flush=True)
