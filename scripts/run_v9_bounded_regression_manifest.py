#!/usr/bin/env python3
"""Run each bounded test file separately and record a machine-readable manifest."""
from __future__ import annotations
import json,subprocess,sys,tempfile
from pathlib import Path
import xml.etree.ElementTree as ET

def main():
    root=Path(__file__).resolve().parents[1]; rows=[]
    for path in sorted((root/"tests").glob("test_*.py")):
        with tempfile.NamedTemporaryFile(suffix=".xml") as tmp:
            result=subprocess.run([sys.executable,"-m","pytest","-q",str(path),f"--junitxml={tmp.name}"],cwd=root,text=True,capture_output=True)
            suite=ET.parse(tmp.name).getroot()
            if suite.tag=="testsuites":suite=next(iter(suite))
            row={"test_group":str(path.relative_to(root)),"collected":int(suite.attrib.get("tests",0)),
                 "failed":int(suite.attrib.get("failures",0))+int(suite.attrib.get("errors",0)),
                 "skipped":int(suite.attrib.get("skipped",0)),"exit_code":result.returncode}
            row["passed"]=row["collected"]-row["failed"]-row["skipped"];rows.append(row)
            if result.returncode:break
    manifest={"schema":"V9_BOUNDED_REGRESSION_MANIFEST_1",
      "source_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True).strip(),
      "groups":rows,"totals":{k:sum(r[k] for r in rows) for k in ("collected","passed","failed","skipped")},
      "all_passed":all(r["exit_code"]==0 for r in rows)}
    out=root/"runs/sn_v9_shared_root_m1_peak/conditioned_premark_1500_N50_v1/regression_manifest.json"
    out.write_text(json.dumps(manifest,indent=2)+"\n")
    if not manifest["all_passed"]:raise SystemExit(1)

if __name__=="__main__":main()
