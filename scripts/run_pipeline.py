"""Execute the current fixed-partition protocol on two CUDA devices in dependency order."""
import importlib.metadata
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import os
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"


def parallel(script,phase,stem):
    jobs=[]
    for worker in [0,1]:
        log=(OUT/f"logs/{stem}_worker{worker}.log").open("w")
        args=[sys.executable,"-u",str(ROOT/"scripts"/script),"--phase",phase,
              "--worker",str(worker),"--workers","2","--device",f"cuda:{worker}"]
        process=subprocess.Popen(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        jobs.append((process,log,args))
    for process,log,args in jobs:
        status=process.wait()
        log.close()
        if status:
            raise subprocess.CalledProcessError(status,args)
    print(json.dumps({"completed_stage":stem,"time_utc":datetime.now(timezone.utc).isoformat()}),flush=True)


def serial(script,phase):
    subprocess.run([sys.executable,"-u",str(ROOT/"scripts"/script),"--phase",phase],cwd=ROOT,check=True)


def main():
    import torch
    (OUT/"logs").mkdir(parents=True,exist_ok=True)
    fingerprint_files = [ROOT/'configs/experiment.json', ROOT/'data/processed/model_arrays.npz',
                         ROOT/'data/evidence/integrated_features.npz', ROOT/'data/evidence/graphs.npz',
                         ROOT/'data/evidence/genes_annotation.csv']
    fingerprint_files += list((ROOT/'metagene').glob('*.py'))
    fingerprint_files += [ROOT/'scripts'/name for name in ['run_benchmarks.py','run_active.py','finalize_revision.py']]
    fingerprints = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in fingerprint_files}
    (OUT/'input_fingerprints.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
    environment={"started_utc":datetime.now(timezone.utc).isoformat(),"python":sys.version,
        "packages":{name:importlib.metadata.version(name) for name in
                    ["torch","numpy","pandas","scipy","scikit-learn","xgboost","matplotlib"]},
        "cuda":torch.version.cuda,"gpus":[torch.cuda.get_device_name(i) for i in [0,1]],
        "benchmark_fits":sum(len(v) * (5 if k in {"pair_holdout", "gene_holdout"} else 3) for k,v in json.loads((ROOT/"configs/experiment.json").read_text())["protocols"].items()),"al_sequences":20,"al_fits":260,"final_fits":5,
        "llm_api_calls":0}
    (OUT/"runtime.json").write_text(json.dumps(environment,indent=2)+"\n")
    parallel("run_benchmarks.py","fit","fixed_fit")
    serial("run_benchmarks.py","select")
    parallel("run_benchmarks.py","evaluate","fixed_evaluate")
    serial("run_benchmarks.py","summarize")
    serial("finalize_revision.py","audit")
    parallel("finalize_revision.py","external","external")
    serial("finalize_revision.py","external-summary")
    parallel("finalize_revision.py","fit-final","final_fit")
    serial("finalize_revision.py","predict")
    parallel("run_active.py","run","fixed_al")
    serial("run_active.py","summarize")
    subprocess.run([sys.executable,"-u",str(ROOT/"scripts/summarize_acquisition.py")],cwd=ROOT,check=True)
    environment["completed_utc"]=datetime.now(timezone.utc).isoformat()
    (OUT/"runtime.json").write_text(json.dumps(environment,indent=2)+"\n")
    print(json.dumps({"completed":"fixed-partition revision"}),flush=True)


if __name__=="__main__":
    if '--legacy-gi' in sys.argv:
        sys.argv.remove('--legacy-gi')
        main()
    else:
        import argparse
        sys.path.insert(0,str(ROOT))
        from metagene.reproducibility import resolve_devices
        parser=argparse.ArgumentParser(description='Run the current SL experiment using fixed inputs.')
        parser.add_argument('--devices',nargs='+',help='CUDA devices or cpu; defaults to METAGENE_DEVICES or automatic selection.')
        args=parser.parse_args()
        environment={**os.environ,'METAGENE_DEVICES':' '.join(resolve_devices(args.devices))}
        for script in ['prepare_sl_labels.py','run_sl_experiments.py','deliver_sl.py','analyze_sl.py','audit_sl.py','make_sl_report.py','finalize_sl_release.py']:
            subprocess.run([sys.executable,'-u',str(ROOT/'scripts'/script)],cwd=ROOT,env=environment,check=True)
