import json
import gdown

import numpy as np
from scipy.io import loadmat

from task_and_baseline import baseline, build_task_helpers

# Download the dataset
url = "https://drive.google.com/file/d/1BBHVSI4KB-B8OX46eN1Nm4ARCeq6Rui4/view?usp=sharing"
downloaded_file = "challenge.mat"
gdown.download(url, downloaded_file, quiet=False, fuzzy=True)

data = loadmat("challenge.mat", simplify_cells=True)
tx = data["tx"].astype(np.complex128)
rx = data["rx"].astype(np.complex128)
Fs = float(data["Fs"])
N, _ = tx.shape

tx_n = tx / (np.sqrt(np.mean(np.abs(tx) ** 2, axis=0, keepdims=True)) + 1e-30)
helpers = build_task_helpers(tx_n, Fs, N)


def your_canceller(tx_n, rx, richardson_steps: int = 2, alpha: float = 0.70):
    """TX IMD3 LS fit + rank-1 external interferer subtraction with Richardson
    compensation for the scorer's double band-pass. See SOLUTION.md."""
    del tx_n
    fit_tx = helpers["fit_tx_prediction"]
    score_filter = helpers["score_filter"]
    n_ch = rx.shape[1]

    tx_pred = fit_tx(rx)
    resid = rx - tx_pred
    resid_band = np.column_stack(
        [score_filter(resid[:, c]) for c in range(n_ch)]
    )

    # Rank-1 projection; formula matches task_and_baseline.rank1_from_band_matrix.
    cov = resid_band.conj().T @ resid_band / resid_band.shape[0]
    _, vecs = np.linalg.eigh(cov)
    shared = resid_band @ vecs[:, -1]
    denom = np.vdot(shared, shared) + 1e-30
    target = np.column_stack(
        [
            (np.vdot(shared, resid_band[:, c]) / denom) * shared
            for c in range(n_ch)
        ]
    )

    # Richardson: solves H·e_hat ≈ target so the scorer's second band-pass pass
    # does not silently attenuate our subtraction.
    e_hat = target.copy()
    for _ in range(richardson_steps):
        h_e = np.column_stack(
            [score_filter(e_hat[:, c]) for c in range(n_ch)]
        )
        e_hat = e_hat + alpha * (target - h_e)

    return rx - tx_pred - e_hat


print("\n=== Baseline ===")
baseline_reds, baseline_avg = helpers["score"](
    rx, baseline(tx_n, rx, helpers["fit_tx_prediction"]), label="baseline"
)

print("=== Your Solution ===")
yours_reds, yours_avg = helpers["score"](rx, your_canceller(tx_n, rx), label="yours")

results = {
    "baseline": {
        "per_channel_db": baseline_reds,
        "average_db": baseline_avg,
    },
    "yours": {
        "per_channel_db": yours_reds,
        "average_db": yours_avg,
    },
}

with open("results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
