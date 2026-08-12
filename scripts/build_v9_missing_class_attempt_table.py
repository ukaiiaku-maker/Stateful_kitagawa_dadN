#!/usr/bin/env python3
"""Collect exact completed H-action trajectories into a raw attempt S-N table."""
import argparse, csv, math
from pathlib import Path

TARGETS=(3.,60.,1.1e3,2.1e4,4e5,7.5e6,1.4e8,2.7e9,5e10,1e12)
Q=(('N10_attempt',-math.log(.9)),('N50_attempt',math.log(2)),('N90_attempt',-math.log(.1)))

def crossing(rows,target):
    n=[float(r['cycles_total']) for r in rows];h=[float(r['H_attempt']) for r in rows]
    for i,x in enumerate(h):
        if x>=target:
            if i==0:return n[0]*target/x
            return n[i-1]+(target-h[i-1])/(x-h[i-1])*(n[i]-n[i-1])
    raise RuntimeError(f'action target {target} not crossed')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--class-name',required=True);ap.add_argument('--option-id',required=True);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    out=[]
    for history in a.root.rglob('sn_stateful_pd_history.csv'):
        rows=list(csv.DictReader(history.open()))
        if not rows or float(rows[-1]['H_attempt'])<Q[-1][1]-1e-12:continue
        sigma=float(rows[-1]['sigma_a_MPa']); vals={name:crossing(rows,t) for name,t in Q}
        target=min(TARGETS,key=lambda x:abs(math.log10(vals['N50_attempt']/x)))
        out.append(dict(material_class=a.class_name,option_id=a.option_id,sigma_a_MPa=sigma,
          sigma_min_MPa=2*sigma*.1/.9,sigma_max_MPa=2*sigma/.9,**vals,target_N50=target,
          target_error_decades=math.log10(vals['N50_attempt']/target),source_history=str(history)))
    unique={r['sigma_a_MPa']:r for r in out};out=sorted(unique.values(),key=lambda r:r['N50_attempt'])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    fields=list(out[0]);
    with a.out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
    print({'rows':len(out),'out':str(a.out)})
if __name__=='__main__':main()
