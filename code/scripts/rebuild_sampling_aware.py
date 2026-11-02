"""Rebuild sampling-aware four-compartment FEM outputs from locked STEP02 parameters.

This is a staging implementation. Source Figure 2 B/D/F transcriptions have been audited; raw replicate data remain unavailable.
Q_receptor_cumulative = Q_receptor_present_pre + Q_samples_withdrawn.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eig,solve
from sampling_aware_fem import ROOT, A, VR, VS, SAMPLE_TIMES, D1, D4, K12, K23, nodes, reg, layer_w
from fit_sampling_aware import sim,names

OUT=ROOT/'outputs'; TAB=OUT/'tables'; MOD=OUT/'models'; FIG=OUT/'figures'; ST=OUT/'step02'
for d in (TAB,MOD,FIG,ST):d.mkdir(parents=True,exist_ok=True)
params=pd.read_csv(ST/'source_reconciled_all_cases_parameters.csv')
obs=pd.read_csv(ROOT/'data/processed/fem_fit_data_master.csv')
selected={'Free KC':'Case 1','KLZ-NPs':'Case 2'}
forms=['Free KC','KLZ-NPs']
source_csv=ST/'source_reconciled_all_cases_profiles.csv'
raw_profiles=pd.read_csv(source_csv)
records=[]; performance=[];full_states={};profile_rows=[]
for f in forms:
 for case in ('Case 1','Case 2'):
  p=params[(params.formulation==f)&(params.model_case==case)].iloc[0]
  v=[float(p[k]) for k in names]
  pr,states=sim(v,return_states=True)
  expected=raw_profiles[(raw_profiles.formulation==f)&(raw_profiles.model_case==case)]
  diff=np.max(np.abs(pr['Q_receptor_cumulative_ug_cm2'].to_numpy()-expected['Q_receptor_cumulative_ug_cm2'].to_numpy()))
  if diff>1e-6:raise AssertionError(f'Cached result differs from regenerated sampling-aware FEM: {f},{case},{diff}')
  pr.insert(0,'formulation',f);pr.insert(1,'model_case',case); profile_rows.append(pr)
  full_states[(f,case)]=states
  dd=obs[(obs.formulation==f)&(obs.compartment=='R')]
  exp=dd.sort_values('time_h'); qq=pr.set_index('time_h').loc[exp.time_h,'Q_receptor_cumulative_ug_cm2'].to_numpy(); y=exp.observed_ug_cm2.to_numpy()
  rmse=np.sqrt(np.mean((qq-y)**2));mae=np.mean(abs(qq-y));r2=1-np.sum((qq-y)**2)/np.sum((y-y.mean())**2)
  perf={'formulation':f,'model_case':case,'RMSE_receptor_ug_cm2':rmse,'MAE_receptor_ug_cm2':mae,'R2_receptor':r2,'n_receptor':len(y),'objective':p.objective};performance.append(perf)
  pr.to_csv(MOD/f'{f.replace(" ","_")}_{case.replace(" ", "_")}_sampling_aware_profile.csv',index=False)
  for row in pr.to_dict(orient='records'):records.append({'formulation':f,'model_case':case,**row})
  p24=pr[pr.time_h==24].iloc[0]
  for c,x in [('EP','Q3_EP_ug_cm2'),('R','Q_receptor_cumulative_ug_cm2'),('SC','Q2_SC_ug_cm2')]:
   actual=float(obs[(obs.formulation==f)&(obs.compartment==c)].query('time_h==24').iloc[0].observed_ug_cm2)
   model=float(p24[x]); records.append({'formulation':f,'model_case':case,'time_h':24.,'compartment':c,'experimental_ug_cm2':actual,'predicted_ug_cm2':model,'signed_error_ug_cm2':model-actual,'kind':'endpoint'})
pr_all=pd.concat(profile_rows,ignore_index=True);pr_all.to_csv(MOD/'case1_case2_sampling_aware_profiles.csv',index=False)
performance=pd.DataFrame(performance);performance.to_csv(TAB/'case1_case2_performance_comparison.csv',index=False)
params.to_csv(MOD/'sampling_aware_case1_case2_parameters.csv',index=False)
for i,case in ((1,'Case 1'),(2,'Case 2')):params[params.model_case==case].to_csv(MOD/f'case{i}_fitted_parameters.csv',index=False)
params[params.apply(lambda r:selected[r.formulation]==r.model_case,axis=1)].to_csv(MOD/'selected_final_model_parameters.csv',index=False)
performance[performance.apply(lambda r:selected[r.formulation]==r.model_case,axis=1)].to_csv(MOD/'selected_final_model_metrics.csv',index=False)
selected_pr=pr_all[pr_all.apply(lambda r:selected[r.formulation]==r.model_case,axis=1)].copy(); selected_pr.to_csv(MOD/'selected_final_model_dense_simulation.csv',index=False)
endpoint=pd.DataFrame([r for r in records if r.get('kind')=='endpoint' and selected[r['formulation']]==r['model_case']]);endpoint.to_csv(TAB/'endpoint_experimental_vs_model_24h_measured_compartments.csv',index=False)

# Updated compartment balance, including mass removed by sampling. Amounts normalized to area.
comp_bal=[];mass_bal=[]
for f in forms:
 case=selected[f];p=params[(params.formulation==f)&(params.model_case==case)].iloc[0]
 ss=selected_pr[(selected_pr.formulation==f)&(selected_pr.time_h==24)].iloc[0]
 Q=float(p.Q0_eff_ug_cm2)
 for label,key,depth in [('Donor','Q1_donor_ug_cm2',1./A),('SC','Q2_SC_ug_cm2',.002),('EP/dermis','Q3_EP_ug_cm2',.088),('Receptor (present)','Q4_receptor_ug_cm2',VR/A),('Collected samples','Q_removed_ug_cm2',None)]:
  val=float(ss[key]);comp_bal.append({'formulation':f,'model_case':case,'compartment':label,'amount_ug_cm2':val,'mean_concentration_ug_mL_or_cm3':(val/depth if depth else np.nan),'percent_initial':100*val/Q})
 total=sum(float(ss[k]) for k in ('Q1_donor_ug_cm2','Q2_SC_ug_cm2','Q3_EP_ug_cm2','Q4_receptor_ug_cm2','Q_removed_ug_cm2'))
 mass_bal.append({'formulation':f,'model_case':case,'initial_ug_cm2':Q,'terminal_sum_including_samples_ug_cm2':total,'absolute_error_ug_cm2':total-Q,'relative_error':(total-Q)/Q})
pd.DataFrame(comp_bal).to_csv(TAB/'endpoint_compartment_balance_24h_long.csv',index=False)
pd.DataFrame(mass_bal).to_csv(TAB/'mass_conservation_including_samples.csv',index=False)

# Mesh check: explicitly assemble the same four uncoupled layer meshes at 80/160/320 elements.
lengths=[1/A,.002,.088,VR/A]
def mesh_predict(counts,v):
 x=[];regions={};edges=[];current=0.0
 for layer,(N,L) in enumerate(zip(counts,lengths),1):
  s=len(x)
  vals=current+np.linspace(0,L,N+1)
  x.extend(vals)
  regions[layer]=np.arange(s,s+N+1)
  for j in range(s,s+N):edges.append((layer,j,j+1,L/N))
  current+=L
 n=len(x);M=np.zeros((n,n));K=np.zeros((n,n));lw={k:np.zeros(n) for k in range(1,5)}
 q0,d2,d3,p12,p23,p34,k34=v
 D={1:D1,2:d2,3:d3,4:D4}
 for layer,i,j,h in edges:
  for u,ni in enumerate((i,j)):
   lw[layer][ni]+=h/2
   for w,nj in enumerate((i,j)):
    M[ni,nj]+=h*(2 if u==w else 1)/6
    K[ni,nj]+=D[layer]*(1 if u==w else -1)/h
 F=np.zeros_like(M)
 def itf(l,r,k,P):
  i=regions[l][-1];j=regions[r][0]
  F[i,i]-=k/P;F[i,j]+=k;F[j,i]+=k/P;F[j,j]-=k
 itf(1,2,K12,p12);itf(2,3,K23,p23);itf(3,4,k34,p34)
 O=solve(M,F-K,assume_a='pos',check_finite=False)
 w,z=eig(O,check_finite=False);inv=solve(z,np.eye(n),check_finite=False)
 y=np.zeros(n);y[regions[1]]=q0/lengths[0]
 removed=0;t0=0
 for t in SAMPLE_TIMES:
  y=np.real_if_close(z@(np.exp(w*(t-t0))*(inv@y)),tol=1e4)
  if np.iscomplexobj(y):raise ValueError('complex mesh result')
  if t < 24:
   removed+=(lw[4]@y)*VS/VR
   y[regions[4]]*=1-VS/VR
  t0=t
 return {'Q_cumulative_receptor_24h_ug_cm2':float(lw[4]@y+removed),'Q_SC_24h_ug_cm2':float(lw[2]@y),'Q_EP_24h_ug_cm2':float(lw[3]@y)}
mesh_vals=[]
for f in forms:
 case=selected[f];row=params[(params.formulation==f)&(params.model_case==case)].iloc[0];v=[float(row[k]) for k in names]
 for label,count in [('Coarse',[10,5,50,15]),('Base',[20,10,100,30]),('Fine',[40,20,200,60])]:
  print('mesh',f,label,flush=True)
  rec=mesh_predict(count,v)
  for name,val in rec.items():mesh_vals.append({'formulation':f,'mesh':label,'elements':sum(count),'output':name,'value_ug_cm2':val})
mesh_df=pd.DataFrame(mesh_vals)
for (f,name),sub in mesh_df.groupby(['formulation','output']):
 fine=float(sub[sub.mesh=='Fine'].value_ug_cm2.iloc[0]);idx=sub.index
 mesh_df.loc[idx,'fine_mesh_ug_cm2']=fine;mesh_df.loc[idx,'relative_difference_percent']=100*np.abs(mesh_df.loc[idx,'value_ug_cm2']-fine)/max(abs(fine),1e-30)
mesh_df.to_csv(TAB/'mesh_independence_summary.csv',index=False)
# Cross-check base mesh against canonical supplied discretization.
for f in forms:
 mm=mesh_df[(mesh_df.formulation==f)&(mesh_df.mesh=='Base')]
 s=selected_pr[(selected_pr.formulation==f)&(selected_pr.time_h==24)].iloc[0]
 for k,col in [('Q_cumulative_receptor_24h_ug_cm2','Q_receptor_cumulative_ug_cm2'),('Q_SC_24h_ug_cm2','Q2_SC_ug_cm2'),('Q_EP_24h_ug_cm2','Q3_EP_ug_cm2')]:
  err=abs(float(mm[mm.output==k].value_ug_cm2.iloc[0])-float(s[col]))
  if err>2e-4:raise AssertionError(f'Independent base mesh disagrees: {f},{k},{err}')

# Figure writing. Keep all output PNGs at 300 dpi and vector SVG/PDF.
def save(name):
 plt.tight_layout()
 for ext in ('png','svg','pdf'):plt.savefig(FIG/(name+'.'+ext),dpi=300,bbox_inches='tight')
 plt.close()
colors={'Free KC':'tab:blue','KLZ-NPs':'tab:orange'}
obsR=obs[obs.compartment=='R']

# Main Figure 3: receptor cumulative (not instantaneous) plus skin retention.
fig,axs=plt.subplots(2,2,figsize=(12,8))
for j,f in enumerate(forms):
 ax=axs[0,j];g=selected_pr[selected_pr.formulation==f]
 exp=obsR[obsR.formulation==f].sort_values('time_h')
 ax.plot(g.time_h,g.Q_receptor_cumulative_ug_cm2,'-',color=colors[f],label='Sampling-aware cumulative model')
 ax.scatter(exp.time_h,exp.observed_ug_cm2,color='black',marker='o',label='Published/digitized cumulative data')
 ax.set(xlabel='Time (h)',ylabel='Cumulative receptor recovery (ug/cm2)',title=('A. ' if j==0 else 'B. ')+f)
 ax.legend(fontsize=8)
for j,(cc,title,key) in enumerate([('SC','C. Stratum corneum','Q2_SC_ug_cm2'),('EP','D. Viable EP/dermis','Q3_EP_ug_cm2')]):
 ax=axs[1,j];ex=[float(obs[(obs.formulation==f)&(obs.compartment==cc)].iloc[0].observed_ug_cm2) for f in forms]; mo=[float(selected_pr[(selected_pr.formulation==f)&(selected_pr.time_h==24)].iloc[0][key]) for f in forms];xx=np.arange(2)
 ax.bar(xx-.18,ex,width=.35,label='Experimental');ax.bar(xx+.18,mo,width=.35,label='Model, 24 h');ax.set_xticks(xx,forms);ax.set_ylabel('Retention (ug/cm2)');ax.set_title(title);ax.legend(fontsize=8)
save('figure3_sampling_aware_selected_models')

# Figure 4: five time point concentration profiles in SC/EP from model states.
fig,axarr=plt.subplots(2,2,figsize=(12,8))
for j,f in enumerate(forms):
 states=full_states[(f,selected[f])]
 for kk,(layer,name) in enumerate([(2,'SC'),(3,'EP/dermis')]):
  ax=axarr[j,kk];ids=reg[layer]; xx=nodes.iloc[ids].x_global_cm.to_numpy();xx=(xx-(.4545454545 if layer==2 else .4565454545))*1e4
  for t in (1,4,8,12,24):ax.plot(xx,states[float(t)][ids],label=f'{t} h')
  ax.set(xlabel=f'Position within {name} (micrometers)',ylabel='Modeled concentration (ug/cm3)',title=f'{"A" if j==0 and kk==0 else "B" if j==0 else "C" if kk==0 else "D"}. {f}, {name}')
  ax.legend(fontsize=8)
save('figure4_sampling_aware_skin_profiles')

# S1: mesh sensitivity
fig,ax=plt.subplots(figsize=(8,4.5))
for f in forms:
 s=mesh_df[(mesh_df.formulation==f)&(mesh_df.mesh!='Fine')]
 for name,g in s.groupby('output'):
  ax.plot(g.elements,g.relative_difference_percent,marker='o',label=f+' '+name.replace('Q_','').replace('_24h_ug_cm2',''))
ax.set_yscale('log');ax.set(xlabel='Mesh elements',ylabel='Difference from 320-element mesh (%)',title='S1. Sampling-aware mesh sensitivity');ax.legend(fontsize=7,ncol=2)
save('S1_sampling_aware_mesh_sensitivity')

# S2 experimental vs modeled terminal endpoints
fig,axarr=plt.subplots(1,2,figsize=(11,4.5))
for j,f in enumerate(forms):
 ax=axarr[j];s=endpoint[endpoint.formulation==f].set_index('compartment').loc[['R','SC','EP']]
 x=np.arange(3);ax.bar(x-.18,s.experimental_ug_cm2,width=.35,label='Experimental');ax.bar(x+.18,s.predicted_ug_cm2,width=.35,label='Model cumulative R');ax.set_xticks(x,['Receptor (cumulative)','SC','EP/dermis']);ax.set(title=f, ylabel='Amount (ug/cm2)');ax.tick_params(axis='x',labelrotation=12);ax.legend(fontsize=8)
save('S2_sampling_aware_endpoint_comparison')

# S3 cumulative receptor residuals at receptor fitting times.
fig,axarr=plt.subplots(1,2,figsize=(10,4))
for j,f in enumerate(forms):
 ax=axarr[j];g=selected_pr[selected_pr.formulation==f].set_index('time_h');exp=obsR[obsR.formulation==f].sort_values('time_h');time=exp.time_h.to_numpy();res=g.loc[time,'Q_receptor_cumulative_ug_cm2'].to_numpy()-exp.observed_ug_cm2.to_numpy()
 ax.axhline(0,color='grey',lw=.9);ax.plot(time,res,'o-',color=colors[f]);ax.set(title=f,xlabel='Time (h)',ylabel='Model minus experimental (ug/cm2)')
save('S3_sampling_aware_receptor_residuals')

# S4 full domain profiles; note independent geometry for each layer.
fig,axarr=plt.subplots(1,2,figsize=(11,4.7))
for j,f in enumerate(forms):
 ax=axarr[j];states=full_states[(f,selected[f])]
 for t in (4,12,24):ax.plot(nodes.x_global_cm.to_numpy(),states[float(t)],label=f'{t} h')
 for x in (.4545454545,.4565454545,.5445454545):ax.axvline(x,color='0.5',lw=.6,linestyle=':')
 ax.set(title=f,xlabel='Full equivalent computational position (cm)',ylabel='Model concentration (ug/cm3)');ax.legend(fontsize=8)
save('S4_sampling_aware_full_domain_profiles')

# S5: model-only physical receptor concentration and mass, not falsely inferred experimental concentration.
fig,axarr=plt.subplots(1,2,figsize=(11,4.7))
for j,f in enumerate(forms):
 ax=axarr[j];g=selected_pr[selected_pr.formulation==f];ax.plot(g.time_h,g.C_receptor_pre_sampling_ug_mL,color=colors[f],marker='o',label='Pre-sample receptor concentration (model only)')
 ax.set(title=f,xlabel='Time (h)',ylabel='Concentration at sampling (ug/mL)');ax.legend(fontsize=8)
save('S5_sampling_aware_model_only_receptor_concentration')

# S6 terminal modeled drug amounts including collected samples.
fig,axarr=plt.subplots(1,2,figsize=(11,4.7))
for j,f in enumerate(forms):
 ax=axarr[j];g=pd.DataFrame(comp_bal);g=g[g.formulation==f];ax.bar(range(len(g)),g.amount_ug_cm2,color=colors[f]);ax.set_yscale('log');ax.set_xticks(range(len(g)),['Donor','SC','EP/dermis','Receptor','Removed'],rotation=15);ax.set(title=f,ylabel='Amount at 24 h (ug/cm2), log scale')
save('S6_sampling_aware_compartment_and_removed_amount')

# S7 model-predicted mean concentrations of remaining compartments, not removed sample amount.
fig,axarr=plt.subplots(1,2,figsize=(11,4.7))
for j,f in enumerate(forms):
 ax=axarr[j];g=pd.DataFrame(comp_bal);g=g[(g.formulation==f)&(g.compartment!='Collected samples')];ax.bar(range(len(g)),g.mean_concentration_ug_mL_or_cm3,color=colors[f]);ax.set_xticks(range(len(g)),['Donor','SC','EP/dermis','Receptor'],rotation=15);ax.set(title=f,ylabel='Model mean concentration (ug/cm3 or ug/mL)')
save('S7_sampling_aware_mean_concentrations')

# Figure 2: redraw from auditable source-reconciled experiment-derived values.
# These are approximate bar-height digitizations (except source-reported 24 h receptor means).
raw=pd.read_csv(ROOT/'data/provenance/STEP02_CANONICAL_PORCINE_PROFILE.csv')
fig,axs=plt.subplots(1,3,figsize=(14.6,4.75))
for j,(cc,title) in enumerate([('R','A. Receptor, cumulative permeation'),('SC','B. Stratum corneum'),('EP','C. Viable epidermis/dermis')]):
 ax=axs[j]
 for f in forms:
  g=raw[(raw.compartment==cc)&(raw.formulation==f)].sort_values('time_h')
  ax.plot(g.time_h,g.canonical_ug_cm2,'o-',color=colors[f],label=f)
 ax.set(xlabel='Time (h)',ylabel='Digitized amount (ug/cm2)',title=title)
 ax.set_xticks([4,6,8,12,24] if cc=='R' else [0,4,8,12,24])
 ax.legend(fontsize=8)
 ax.grid(axis='y',alpha=.18)
fig.suptitle('Porcine skin: digitized source Fig. 2B/D/F (approximate values)',fontsize=11)
save('figure2_experimental_input_profiles')

checks={'mesh_base_agreement':'passed','cumulative_recovery_includes_removed_samples':True,'source_digitization_mismatches_after_step02':0,'selected_cases':selected,'sampler_VR_mL':VR,'sampled_VS_mL':VS,'times_h':list(SAMPLE_TIMES),'total_mass_error':mass_bal,'performance':performance.to_dict(orient='records')}
import json
(ST/'SOURCE_RECONCILIATION_VALIDATION.json').write_text(json.dumps(checks,indent=2,default=float))
print('REBUILD COMPLETE:',ROOT)
for row in mass_bal:print('Mass',row)
