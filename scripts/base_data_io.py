"""Load measured-pair arrays plus independent gene covariates."""
import numpy as np
import pandas as pd
from metagene.training import PairData

def load_data(root):
    arrays = dict(np.load(root / "data/processed/model_arrays.npz"))
    enriched = root / "data/evidence/integrated_features.npz"
    arrays["gene_features"] = np.load(enriched)["features"]
    data = PairData(**arrays, gi_target_name="additive_gi_proxy")
    return data, pd.read_csv(root / "data/processed/model_pairs.csv")
