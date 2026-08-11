#!/usr/bin/env python3
"""Compact diagnostic plots for the two non-production Peak PD controls."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("runs/v9_pd_production")
OUT = ROOT / "figures" / "Peak_bridge_diagnostics"
OUT.mkdir(parents=True, exist_ok=True)
p4 = ROOT / "packages/Peak_blunt_4000MPa_N10B_generation_80f16e02"
p12 = ROOT / "packages/Peak_blunt_12000MPa_N3M_generation_0c0a2960"
d4 = json.loads((p4 / "terminal_exact_diagnostics.json").read_text())
h12 = pd.read_csv(p12 / "sn_stateful_pd_history.csv")
tail = pd.read_csv(ROOT / "diagnostics/Peak_4000MPa_K2_suppressed_tail/Peak_4000MPa_K2_suppressed_tail_matched_checkpoints.csv")

def save(fig, name):
    fig.tight_layout(); fig.savefig(OUT / name, dpi=190); plt.close(fig)

# 1 raw cleavage versus gated birth
fig, ax = plt.subplots(figsize=(6.2,4.0))
raw = [d4["raw_cleavage_rate_peak_s"], h12.pd_nucleation_rate_max_s.iloc[-1]]
gated = [d4["actual_site_birth_rate_peak_per_cycle"], h12.pd_birth_rate_max_per_cycle.iloc[-1]]
x=np.arange(2); w=.34
ax.bar(x-w/2,raw,w,label="raw cleavage [s$^{-1}$]");ax.bar(x+w/2,gated,w,label="gated birth [cycle$^{-1}$]")
ax.set_yscale('log');ax.set_xticks(x,["4 GPa\nK2-suppressed","12 GPa\noverstress"]);ax.legend();ax.set_ylabel("rate (mixed units; see legend)")
save(fig,"01_raw_cleavage_vs_gated_birth.png")

# 2 memory and completion
fig,axs=plt.subplots(1,2,figsize=(9,3.8))
axs[0].loglog(h12.cycles_total,h12.pd_delivery_memory_max,label="12 GPa")
axs[1].loglog(h12.cycles_total,np.maximum(h12.pd_completion_max,1e-320),label="12 GPa")
for label,g in tail.groupby('checkpoint'):
    n=float(g.cycles.iloc[0]); axs[0].scatter(n,g.delivery_memory.max(),c='C0',marker='s');axs[1].scatter(n,max(g.K2_completion.max(),1e-320),c='C0',marker='s')
axs[0].set_ylabel(r"delivery memory $\Lambda$");axs[1].set_ylabel("K=2 completion Q")
for a in axs:a.set_xlabel("cycles N");a.grid(True,which='both',alpha=.25)
save(fig,"02_delivery_memory_K2_completion.png")

# 3 cumulative action and no-embryo survival diagnostic
fig,axs=plt.subplots(1,2,figsize=(9,3.8))
mx=tail.groupby('cycles').cumulative_action.max(); n=mx.index.to_numpy(float);H=mx.to_numpy(float)
axs[0].loglog(n,np.maximum(H,1e-320),'o-',label='4 GPa max node action')
axs[0].loglog(h12.cycles_total,np.maximum(h12.pd_expected_births_cumulative,1e-320),label='12 GPa expected births')
axs[1].semilogx(n,np.exp(-H),'o-',label='4 GPa exp(-max node action)')
axs[0].set_ylabel("cumulative diagnostic");axs[1].set_ylabel("no-embryo survival diagnostic")
for a in axs:a.set_xlabel("cycles N");a.grid(True,which='both',alpha=.25);a.legend(fontsize=8)
save(fig,"03_cumulative_hazard_survival.png")

# 4 site-resolved remaining action and wait
g=tail[tail.checkpoint=='N1e14'].sort_values('remaining_threshold_action')
fig,axs=plt.subplots(1,2,figsize=(9,3.8)); ids=g.site_id.to_numpy()
axs[0].semilogy(ids,np.maximum(g.remaining_threshold_action,1e-320),'.',ms=2)
axs[1].semilogy(ids,np.maximum(g.instantaneous_wait_cycles,1e-320),'.',ms=2)
axs[0].set_ylabel("remaining threshold action");axs[1].set_ylabel("instantaneous wait [cycles]")
for a in axs:a.set_xlabel("persistent site ID");a.grid(True,alpha=.2)
save(fig,"04_site_remaining_action_wait.png")

# 5 topology evolution
fig,ax=plt.subplots(figsize=(7,4));ax.semilogx(h12.cycles_total,h12.pd_bond_damage_max,label='max bond damage')
ax.semilogx(h12.cycles_total,h12.pd_broken_bonds,label='broken bonds')
ax.semilogx(h12.cycles_total,h12.pd_active_front.astype(int),label='front captured flag')
ax.axvline(1528000,color='k',ls='--',lw=1,label='root connection')
ax.set_xlabel('cycles N');ax.set_ylabel('damage / count / flag');ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.25)
save(fig,"05_damage_topology.png")

# 6 event cycles versus stress, explicitly diagnostic
fig,ax=plt.subplots(figsize=(6.5,4.2))
events={'embryo':81050.8859089888,'stable site':81053.2753569623,'root connection':1528000.0}
for i,(name,N) in enumerate(events.items()):ax.scatter(N,12000,marker=['o','s','^'][i],s=55,label=name)
ax.scatter(1e10,4000,marker='>',facecolors='none',edgecolors='C0',s=70,label='4 GPa runout diagnostic')
ax.scatter(3e6,12000,marker='x',c='r',s=70,label='12 GPa front capture censored')
ax.set_xscale('log');ax.set_xlabel('cycles N');ax.set_ylabel('stress amplitude [MPa]');ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.25)
ax.set_title('Diagnostic controls — not a physical S–N bracket')
save(fig,"06_stress_event_diagnostics.png")

# 7 acceleration validation
c=json.loads((p4/'v9_pd_high_cycle_controller.json').read_text()); rows=[]
for m in c['mode_history']:
    if m['mode']=='exact_private_window' and m['accepted']:
        rows.append([m['cycles'],m['detail']['state_error'],m['detail']['log_hazard_error'],m['detail']['projected_cycles_per_exact_map']])
a=np.asarray(rows,float);fig,axs=plt.subplots(1,2,figsize=(9,3.8))
axs[0].semilogy(np.arange(len(a)),a[:,1],label='state error');axs[0].semilogy(np.arange(len(a)),np.maximum(a[:,2],1e-320),label='hazard error');axs[0].legend()
axs[1].plot(np.arange(len(a)),a[:,3]);axs[1].axhline(16,color='r',ls='--',label='admission floor');axs[1].legend()
axs[0].set_ylabel('validation error');axs[1].set_ylabel('accepted cycles / exact map')
for z in axs:z.set_xlabel('accepted exact window');z.grid(True,alpha=.25)
save(fig,"07_acceleration_validation.png")

print(OUT)
