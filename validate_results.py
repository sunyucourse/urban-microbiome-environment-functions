"""Compare regenerated tables with preserved historical reference results."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
KEYS = {
    "spatial": ["index", "env"], "category_screen": ["category", "index"],
    "twostage": ["index", "env"], "category_twostage": ["category", "index"],
    "attribution": ["index", "block"], "composites": ["x", "y"],
    "replication_factor": ["index", "env"], "replication_category": ["category", "index"],
    "replication_concordance": ["index_discovery", "env_discovery"],
    "replication_enrichment": ["index_discovery", "env_discovery"], "S14": ["factor"],
}

def validate(output):
    out = Path(output).resolve()
    comparisons = []
    for reference in sorted((ROOT/"reference_results").glob("*.csv")):
        actual = pd.read_csv(out/reference.name)
        expected = pd.read_csv(reference)
        keys = next(v for k,v in KEYS.items() if reference.name.startswith(k))
        actual = actual.sort_values(keys).reset_index(drop=True)
        expected = expected.sort_values(keys).reset_index(drop=True)
        assert set(actual.columns) == set(expected.columns), reference.name
        actual = actual[expected.columns]
        # LMG permutation order has no scientific meaning.
        if "sig_categories" in actual:
            for df in [actual, expected]:
                df["sig_categories"] = df.sig_categories.fillna("").map(
                    lambda s: "|".join(sorted(str(s).split("|"))))
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                      check_exact=False, rtol=1e-7, atol=1e-10)
        comparisons.append({"file": reference.name, "rows": len(actual),
                            "columns": len(actual.columns), "status": "PASS"})
    fac = pd.read_csv(out/"twostage_v7.csv")
    cat = pd.read_csv(out/"category_twostage_v7.csv")
    conc = pd.read_csv(out/"replication_concordance_v8.csv")
    amr = conc[conc["index"] == "sw_amr_depth_pm"]
    agree = int((np.sign(amr.rho_discovery) == np.sign(amr.rho_sewage)).sum())
    summary = {
        "confirmed_factor_pairs": int(((fac.fdr_sample < .10)&fac.direction_match).sum()),
        "confirmed_category_pairs": int(((cat.fdr_sample < .10)&cat.direction_match).sum()),
        "matched_pairs": len(conc), "amr_depth_direction_agreement": [agree, len(amr)],
        "nominal_binomial_p": stats.binomtest(agree, len(amr), .5).pvalue,
    }
    assert summary["confirmed_factor_pairs"] == 47
    assert summary["confirmed_category_pairs"] == 14
    assert summary["matched_pairs"] == 185
    assert summary["amr_depth_direction_agreement"] == [29,37]
    report = {"scope": "Post-profiling only; not an end-to-end raw-read validation",
              "tolerance": {"rtol": 1e-7, "atol": 1e-10},
              "comparisons": comparisons, "summary": summary}
    (out/"validation_report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
    return report

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results")
    validate(parser.parse_args().output)
