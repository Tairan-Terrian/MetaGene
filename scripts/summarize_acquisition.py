"""Describe realized acquisition scores without changing any frozen decision."""
import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/active_learning"


def main():
    records=[]
    for run_path in sorted((OUT/"runs").glob("*/run.json")):
        run=json.loads(run_path.read_text())
        for round_number in range(run["query_rounds"]):
            frame=pd.read_csv(run_path.parent/f"round_{round_number:02d}/query_decisions.csv")
            chosen=frame[frame.selected_rank.gt(0)]
            cutoff=chosen.acquisition_score.min()
            records.append({"strategy":run["strategy"],"seed":run["seed"],"round":round_number,
                "available_pairs":len(frame),"selected_pairs":len(chosen),
                "score_zero_fraction":float(frame.acquisition_score.eq(0).mean()),
                "selected_score_zero_fraction":float(chosen.acquisition_score.eq(0).mean()),
                "distinct_score_fraction":frame.acquisition_score.nunique()/len(frame),
                "cutoff_score":cutoff,
                "n_pairs_tied_at_cutoff":int(frame.acquisition_score.eq(cutoff).sum()),
                "n_selected_at_cutoff":int(chosen.acquisition_score.eq(cutoff).sum()),
                "acquisition_all_zero":bool(frame.acquisition_score.eq(0).all()),
                "method_scope":"Post hoc score diagnostic; no hidden outcome labels read and no decisions changed"})
    frame=pd.DataFrame(records)
    frame.to_csv(OUT/"acquisition_diagnostics.csv",index=False)
    summary=frame.groupby("strategy").agg(rounds=("round","size"),
        mean_score_zero_fraction=("score_zero_fraction","mean"),
        mean_selected_score_zero_fraction=("selected_score_zero_fraction","mean"),
        all_zero_rounds=("acquisition_all_zero","sum"),
        mean_cutoff_ties=("n_pairs_tied_at_cutoff","mean"))
    summary.to_csv(OUT/"acquisition_diagnostic_summary.csv")
    print(summary.to_string())


if __name__=="__main__":
    main()
