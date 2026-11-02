"""Sampling-aware, four-layer, linear one-dimensional finite-element transport solver.

Reimplements original FEM mass and stiffness weak form with the same mesh, fixed D1/D4
and finite interface transfer coefficients; only receptor withdrawal/replacement is added.
"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor, cho_solve
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
A=2.2
VR=6.5
VS=5.0
SAMPLE_TIMES=(1., 2., 4., 6., 8., 12., 24.)
K12=K23=100.
D1=D4=10.0
DATA=ROOT/'data/processed'

nodes=pd.read_csv(DATA/'four_layer_mesh_nodes.csv')
elems=pd.read_csv(DATA/'four_layer_mesh_elements.csv')
print('mesh cols',nodes.columns.tolist(),elems.columns.tolist()) if __name__ == '__main__' else None
n=len(nodes)
ids=nodes.layer_id.to_numpy(dtype=int)
reg={k:np.flatnonzero(ids==k) for k in [1,2,3,4]}
endpts={k:(v[0],v[-1]) for k,v in reg.items()}
length={1:1.0/A,2:0.002,3:0.088,4:VR/A}
M=np.zeros((n,n));Ks={i:np.zeros((n,n)) for i in range(1,5)}
layer_w={i:np.zeros(n) for i in range(1,5)}
for el in elems.itertuples(index=False):
 k=int(el.layer_id)
 i=int(el.left_node);j=int(el.right_node);h=float(el.element_length_cm)
 m=h/6*np.array([[2,1],[1,2]])
 kk=(1./h)*np.array([[1,-1],[-1,1]])
 for aa,ni in enumerate((i,j)):
  layer_w[k][ni]+=h/2
  for bb,nj in enumerate((i,j)):
   M[ni,nj]+=m[aa,bb]
   Ks[k][ni,nj]+=kk[aa,bb]
_mfac=cho_factor(M)
# M symmetric SPD; constant in all simulations

def params_from_dict(p):
 return (p['Q0_eff_ug_cm2'],p['D2_SC_cm2_h'],p['D3_EP_cm2_h'],p['P12_donor_SC'],p.get('P23_SC_EP',1.),p['P34_EP_receptor'],p['K34_cm_h'])

def make_linear_operator(q0,d2,d3,p12,p23,p34,k34):
 K=D1*Ks[1]+d2*Ks[2]+d3*Ks[3]+D4*Ks[4]
 F=np.zeros_like(M)
 def add_interface(kleft,kright,transfer,partition):
  i=endpts[kleft][1];j=endpts[kright][0]
  # positive flux J=transfer*(C_left/partition - C_right)
  F[i,i]-=transfer/partition
  F[i,j]+=transfer
  F[j,i]+=transfer/partition
  F[j,j]-=transfer
 add_interface(1,2,K12,p12)
 add_interface(2,3,K23,p23)
 add_interface(3,4,k34,p34)
 return cho_solve(_mfac,F-K,check_finite=False)

def run_model(params, times_h, sampled=True, rtol=3e-7, atol=5e-9, method='BDF'):
 q0,d2,d3,p12,p23,p34,k34= params if isinstance(params,(tuple,list,np.ndarray)) else params_from_dict(params)
 times=sorted(set(map(float,times_h)))
 if len(times)==0 or min(times)<0: raise ValueError('Invalid times')
 max_time=max(times)
 Aop=make_linear_operator(q0,d2,d3,p12,p23,p34,k34)
 y=np.zeros(n);y[reg[1]]=q0/length[1]
 removed=0.0
 results=[];receptor_nodes=reg[4]
 events=sorted(set(times) | ({x for x in SAMPLE_TIMES if x <= max_time} if sampled else set()) | {0.0})
 t0=0.0
 def record(t,pre):
  parts={k:float(layer_w[k]@pre) for k in [1,2,3,4]}
  rec={
   'time_h':float(t),'Q1_donor_ug_cm2':parts[1], 'Q2_SC_ug_cm2':parts[2],
   'Q3_EP_ug_cm2':parts[3], 'Q4_receptor_ug_cm2':parts[4],
   'Q_removed_ug_cm2':removed,
   'Q_receptor_cumulative_ug_cm2':parts[4]+removed,
   'Q_total_with_samples_ug_cm2':sum(parts.values())+removed,
   'C_receptor_pre_sampling_ug_mL':parts[4]/length[4]
  }
  results.append(rec)
 if 0 in times:record(0.,y)
 for t in events[1:]:
  sol=solve_ivp(lambda tx,cy:Aop@cy,(t0,t),y,method=method,jac=Aop,
             rtol=rtol,atol=atol,t_eval=[t])
  if not sol.success: raise RuntimeError(sol.message)
  y=sol.y[:,-1]
  if t in times:record(t,y)
  # Withdrawal and replacement happen just AFTER the measurement at t;
  # final 24-h sample need not be withdrawn to compute terminal values.
  if sampled and t in SAMPLE_TIMES and t<max_time:
   current=float(layer_w[4]@y)
   removed+=current*(VS/VR)
   y[receptor_nodes]*=(1.-VS/VR)
  t0=t
 return pd.DataFrame(results), y

if __name__=='__main__':
 import time
 p=(485.5,4.16e-6,8.11e-4,1.307,1,0.002985,.840)
 for sampled in (False,True):
  st=time.perf_counter();out,_=run_model(p,[0,1,2,4,6,8,12,24],sampled=sampled)
  print('mode',sampled,'duration',time.perf_counter()-st);print(out.to_string(index=False))
