"""Replot source-reconciled approximate source Figure 2 input profiles (no model refit)."""
from pathlib import Path
import pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2]
raw=pd.read_csv(R/'data/provenance/STEP02_CANONICAL_PORCINE_PROFILE.csv')
cols={'Free KC':'tab:blue','KLZ-NPs':'tab:orange'}
fig,axs=plt.subplots(1,3,figsize=(14.6,4.75))
for j,(comp,title) in enumerate([('R','A. Receptor, cumulative permeation'),('SC','B. Stratum corneum'),('EP','C. Viable epidermis/dermis')]):
 ax=axs[j]
 for f in ('Free KC','KLZ-NPs'):
  g=raw[(raw.compartment==comp)&(raw.formulation==f)].sort_values('time_h')
  ax.plot(g.time_h,g.canonical_ug_cm2,'o-',label=f,color=cols[f],lw=1.5,markersize=4.5)
 ax.set(xlabel='Time (h)',ylabel='Digitized amount (ug/cm2)',title=title)
 ax.set_xticks([4,6,8,12,24] if comp=='R' else [0,4,8,12,24]);ax.tick_params(axis='x',labelsize=9)
 ax.legend(fontsize=8);ax.grid(axis='y',alpha=.18)
fig.suptitle('Porcine skin: digitized source Fig. 2B/D/F (approximate values)',fontsize=11)
fig.tight_layout()
for ext in ('png','pdf','svg'):fig.savefig(R/'outputs/figures'/('figure2_experimental_input_profiles.'+ext),dpi=300,bbox_inches='tight')
plt.close(fig)
print('Figure 2 rebuilt without overlapping x-axis labels.')
