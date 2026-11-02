"""Step 03 diagnostic for sampling-aware, source-reconciled ketoconazole FEM.

Analyzes identifiability without assigning unmeasured experimental variances.
This script NEVER changes the Step 02 fitted parameters or primary data.
Usage from repository root (BLAS threads recommended =1):
    python code/scripts/step03_identifiability_audit.py
"""
from pathlib import Path
import sys, os, json, csv, warnings
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from numpy.linalg import svd
import fit_sampling_aware as fit

ROOT = fit.root
OUT = ROOT/'outputs/step03'
OUT.mkdir(parents=True,exist_ok=True)
PARNAMES=list(fit.names)
BASE=pd.read_csv(ROOT/'outputs/step02/source_reconciled_all_cases_parameters.csv')
FORMS=('Free KC','KLZ-NPs')
REF_TIMES=(4.,6.,8.,12.,24.)
BASELINES={}

class Audit:
 def __init__(self, form, case):
  self.form=form;self.case=case
  row=BASE[(BASE.formulation==form)&(BASE.model_case==f'Case {case}')].iloc[0]
  self.p=np.array([float(row[x]) for x in PARNAMES]);
  self.mask=np.array([0,1,2,3,5,6]) if case==1 else np.arange(7)
  data=fit.obs[fit.obs.formulation==form].sort_values(['time_h','compartment'])
  self.times=data.time_h.to_numpy(float);self.comps=data.compartment.to_numpy();self.y=data.observed_ug_cm2.to_numpy(float)
  scales={c:float(data.loc[data.compartment==c,'observed_ug_cm2'].max()) for c in ('R','SC','EP')}
  self.weights=np.array([1/np.sqrt(5) if c=='R' else 1 for c in self.comps])/np.array([scales[c] for c in self.comps])
  self.cache={}
 def residual(self,z,mask=None,fixed=None):
  if mask is None:mask=self.mask
  p=self.p.copy()
  p[mask]=np.exp(z)
  if fixed is not None:p[4]=float(fixed)
  tup=tuple(np.round(np.log(p),12))
  if tup in self.cache: pred=self.cache[tup]
  else:
   d,_=fit.sim(p)
   idx=d.set_index('time_h')
   pred=np.array([idx.loc[t,'Q_receptor_cumulative_ug_cm2' if c=='R' else 'Q2_SC_ug_cm2' if c=='SC' else 'Q3_EP_ug_cm2'] for t,c in zip(self.times,self.comps)])
   if len(self.cache)<25000:self.cache[tup]=pred
  return self.weights*(pred-self.y)
 def objective(self,z,mask=None,fixed=None):
  r=self.residual(z,mask,fixed)
  return float(np.dot(r,r))
 def fit_from(self,start,mask=None,fixed=None,max_nfev=70):
  if mask is None:mask=self.mask
  low=np.log(fit.lo[mask]);high=np.log(fit.hi[mask]);
  p=np.clip(start.copy(),fit.lo*1.001,fit.hi/1.001)
  if fixed is not None:p[4]=fixed
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',RuntimeWarning)
   res=least_squares(lambda z:self.residual(z,mask,fixed),np.log(p[mask]),bounds=(low,high),max_nfev=max_nfev,method='trf',ftol=5e-9,xtol=5e-9,gtol=5e-9,x_scale='jac')
  new=self.p.copy();new[mask]=np.exp(res.x)
  if fixed is not None:new[4]=fixed
  return new,float(np.dot(res.fun,res.fun)),int(res.nfev),bool(res.success)

# SVD rank, conditioning and column strength in log parameter space
srows=[];crows=[];weakrows=[];objrows=[]
for form in FORMS:
 for case in (1,2):
  m=Audit(form,case)
  x=np.log(m.p[m.mask]);H=1e-3
  jj=[]
  for j in range(len(x)):
   dz=np.zeros_like(x);dz[j]=H
   jj.append((m.residual(x+dz)-m.residual(x-dz))/(2*H))
  J=np.stack(jj,axis=1)
  u,ss,vt=svd(J,full_matrices=False)
  ref=float(ss[0]);cond=float(ss[0]/ss[-1]);obj=m.objective(x)
  effective_rank=int(np.sum(ss/ss[0]>=1e-3))
  srows.append(dict(formulation=form,model_case=f'Case {case}',n_obs=7,n_parameters=len(x),objective=obj,s1=ss[0],s_min=ss[-1],condition_number=cond,rank_rel_1e_minus_3=effective_rank,rank_rel_1e_minus_6=int(np.sum(ss/ss[0]>=1e-6)),boundary_values=';'.join(PARNAMES[k] for k in m.mask if (m.p[k]/fit.lo[k]<1.01 or fit.hi[k]/m.p[k]<1.01))))
  for j,k in enumerate(m.mask):
   crows.append(dict(formulation=form,model_case=f'Case {case}',parameter=PARNAMES[k],value=m.p[k],relative_J_column_norm=float(np.linalg.norm(J[:,j])),weakest_direction_loading=float(vt[-1,j])))
  for rank,val in enumerate(ss,1):weakrows.append(dict(formulation=form,model_case=f'Case {case}',singular_number=rank,value=float(val),value_relative_to_largest=float(val/ref)))
  print('SVD',form,case,'rank@1e-3',effective_rank,'condition',f'{cond:.2e}','obj',obj,flush=True)

pd.DataFrame(srows).to_csv(OUT/'local_svd_summary.csv',index=False)
pd.DataFrame(crows).to_csv(OUT/'parameter_local_sensitivities.csv',index=False)
pd.DataFrame(weakrows).to_csv(OUT/'svd_spectrum.csv',index=False)

# Deterministic perturbation (diagnostic only) of weakly-sensitive receptor-interface coefficients
pert=[]
for form in FORMS:
 for case in (1,2):
  m=Audit(form,case);base=m.objective(np.log(m.p[m.mask]))
  for k in [5,6]:
   for factor in [0.1,0.25,0.5,2.,4.,10.]:
    p=m.p.copy();p[k]*=factor
    if p[k]>=fit.hi[k] or p[k]<=fit.lo[k]:continue
    z=np.log(p[m.mask]);o=m.objective(z)
    pert.append(dict(formulation=form,model_case=f'Case {case}',parameter=PARNAMES[k],factor=factor,baseline_objective=base,new_objective=o,delta_objective=o-base))
pd.DataFrame(pert).to_csv(OUT/'fixed_other_parameter_perturbations.csv',index=False)

# A small controlled multistart check. Optimization is exploratory, not global.
mult=[]
for form in FORMS:
 for case in (1,2):
  m=Audit(form,case)
  starts=[('locked_step02',m.p.copy())]
  for v in [.1,10.]:
   p=m.p.copy();p[6]=np.clip(p[6]*v,fit.lo[6]*1.01,fit.hi[6]/1.01)
   starts.append((f'K34_x_{v:g}',p))
  for label,p in starts:
   try:
    new,o,nfev,suc=m.fit_from(p,max_nfev=15)
    mult.append(dict(formulation=form,model_case=f'Case {case}',start=label,objective=o,nfev=nfev,solver_success=suc,**dict(zip(PARNAMES,new))))
    print('MULTI',form,case,label,f'{o:.8g}',nfev,flush=True)
   except Exception as e:
    mult.append(dict(formulation=form,model_case=f'Case {case}',start=label,objective=np.nan,nfev=-1,solver_success=False,error=repr(e)))
pd.DataFrame(mult).to_csv(OUT/'multistart_diagnostics.csv',index=False)

# Profiled residual objective. Unlike a likelihood, it does not support confidence intervals.
prof=[]
for form in FORMS:
 m=Audit(form,2);mask=np.array([0,1,2,3,5,6]);
 grid=sorted(set([0.2,0.5,1.,2.,4.,8.,12.,float(m.p[4])]))
 seen=[]
 for val in sorted(grid,key=lambda z:abs(np.log(z/m.p[4]))):
  initial=m.p.copy()
  if seen:initial=min(seen,key=lambda v:abs(np.log(v[0]/val)))[1].copy()
  try:
   p,obj,nfev,suc=m.fit_from(initial,mask=mask,fixed=val,max_nfev=15)
   seen.append((val,p));prof.append(dict(formulation=form,P23_fixed=val,profiled_objective=obj,nfev=nfev,optimizer_converged=suc,**dict(zip(PARNAMES,p))))
   print('PROFILE',form,val,f'{obj:.8g}',nfev,suc,flush=True)
  except Exception as e:
   print('PROFILE FAILED',form,val,repr(e),flush=True)
profdf=pd.DataFrame(prof).sort_values(['formulation','P23_fixed']);profdf.to_csv(OUT/'P23_profiled_objective.csv',index=False)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(7.6,4.4))
for form in FORMS:
 d=profdf[profdf.formulation==form]
 ax.plot(d.P23_fixed,d.profiled_objective,'-o',markersize=4,label=form)
ax.set_xscale('log');ax.set_yscale('log');ax.axvline(1,color='gray',linestyle=':',linewidth=1)
ax.set_xlabel('Fixed P23 (dimensionless; nuisance parameters re-optimized)')
ax.set_ylabel('Weighted least-squares objective (not a likelihood)')
ax.grid(alpha=.25,which='both');ax.legend(fontsize=9);fig.tight_layout()
fig.savefig(OUT/'EXPLORATORY_P23_profile_objective.png',dpi=200);fig.savefig(OUT/'EXPLORATORY_P23_profile_objective.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.6,4.4))
for (form,case),g in pd.DataFrame(weakrows).groupby(['formulation','model_case']):
 ax.semilogy(g.singular_number,g.value_relative_to_largest,'o-',label=f'{form}, {case}')
ax.axhline(1e-3,color='gray',linestyle=':',linewidth=0.8)
ax.set_xlabel('Singular-value index');ax.set_ylabel('Singular value / largest singular value')
ax.set_ylim(1e-11,2);ax.set_xticks(range(1,8));ax.grid(alpha=.25,which='both');ax.legend(fontsize=8,ncol=2);fig.tight_layout()
fig.savefig(OUT/'Figure_S8_sensitivity_spectrum.png',dpi=200);fig.savefig(OUT/'Figure_S8_sensitivity_spectrum.svg');plt.close(fig)
summary={
 'status':'DIAGNOSTIC, NOT INDEPENDENT VALIDATION',
 'scope':'Identifiability checks on locked Step02 model, not additional experimental data',
 'observations_per_formulation':7,'receptor_timepoints':5,'skin_endpoints':2,
 'case1_free_parameters':6,'case2_free_parameters':7,
 'uncertainty':'Experimental error variance and graphical digitization distribution unavailable; diagnostics cannot establish confidence intervals, statistical model preference, posterior distributions, or unique coefficients.',
 'svd':srows,'profile_grid':'0.2,0.5,1,2,4,8,12 and each Step02 baseline value; six other parameters refitted with limited iterations.',
 'source_data_unchanged':True,'locked_step02_parameters_unchanged':True}
(OUT/'STEP03_DIAGNOSTICS.json').write_text(json.dumps(summary,indent=2))
print('ANALYSIS COMPLETE',flush=True)
