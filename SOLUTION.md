SMILES-2026 Signal Interference Cancellation Solution

Abstract

We address the SMILES-2026 signal interference cancellation task. 
Using only the 130-dimensional third-order intermodulation (IMD3) basis
prescribed by the scorer plus one spatial rank-1 component, we reach
9.42 dB of in-band interference suppression, up from 4.02 dB for the
provided baseline and 1.42 dB above the 8 dB good tier in the task
brief. The key enabler is an observation absent from the baseline, where the
scorer's band-pass filter is applied twice to the submitted estimate, so the
naive single-pass rank-1 subtraction is attenuated by ~21 % in power near the
band edges. A two-step Richardson iteration over the scorer's band-pass
operator recovers almost the full attenuated energy while remaining inside
the scorer's validity envelope where the explain_ratio ≥ 0.95 and per-channel
err/residual ≤ 0.80.


I. Reproducibility
A. Environment

Create an environment with Python 3 and install required packages:
pip install numpy scipy "gdown<6"

Then run:
python applicant_solution.py

The script is tested on Python 3.12 with numpy 2.4.4, scipy 1.17.1, gdown 5.2.2.
the argument "gdown<6" is required because the  command 
gdown.download(..., fuzzy=True) is used in applicant_solution.py, 
and "fuzzy" was removed in gdown 6.0.

B. Modified Components 

Only the body of your_canceller() in applicant_solution.py
is replaced by proposed solution. task_and_baseline.py and the rest of
applicant_solution.py codes remain identical to the initial repository.

C. Expected output

The output will be saved in results.json and printed as below:

=== Baseline ===
  ch0: 3.98 dB
  ch1: 4.86 dB
  ch2: 3.49 dB
  ch3: 3.74 dB
  Metric [baseline]: 4.02 dB

=== Your Solution ===
  ch0: 10.41 dB
  ch1: 7.85 dB
  ch2: 12.55 dB
  ch3: 6.86 dB
  Metric [yours]: 9.42 dB



II. Scorer Validity Constraints

The task's signal model and scoring rule are given in README.md. 
The scorer in task_and_baseline.decompose_removed_component
imposes two hard validity checks that bound the admissible method space.
The scorer decomposes the removed component rx - rx_hat (band-filtered)
into three parts:

removed_band = band_filter(rx - rx_hat)          # per channel
tx_part      = fit_tx_prediction(rx - rx_hat)    # LS projection on a fixed 130-dim basis
                                                  
residual     = removed_band - tx_part
rank1_part   = rank-1 projection of residual     # 4x4 eigendecomposition
err          = residual - rank1_part

The 130-dim basis is the ten cross-carrier cubic terms
tx[:, a]² · conj(tx[:, b]) , where pairs of carriers within and across the
three transmitter power levels, each shifted by 13 integer lags in
[-6, +6]. This is the standard basis for passive-intermodulation
cancellation in FDD carrier-aggregation transceivers [1] and for external
PIM in FDD MIMO arrays [2].

The rank-1 stage is classical subspace projection as developed for
spatially-coherent RFI rejection in radio-astronomy arrays [3].

Validity requires

explain_ratio = 1 - E[|err|²] / E[|removed_band|²] ≥ 0.95
per channel c: E[|err[:, c]|²] ≤ 0.80 · E[|residual_band[:, c]|²]

Any violation forces every per-channel score to 0 dB. The operational
consequence is that whatever we subtract must lie, in-band, inside

span(Φ) ⊕ {rank-1 spatial subspace}

where Φ is the 130-column basis above. Energy outside this span flows
into variable err and is penalised.


III. Baseline Analysis

The provided baseline rx - fit_tx_prediction(rx) performs only the
TX-IMD3 stage and leaves the external spatially coherent interferer
E[n, c] untouched. The 4x4 band-filtered covariance of the
post-baseline residual has eigenvalue spectrum

[λ₀, λ₁, λ₂, λ₃] = [-37, -43, -53, -55] dB,
                   trace = 2.39e-4

A dominant rank-1 direction sits 17 dB above the noise floor. An
idealised rank-1 subtraction would drop the band-filtered residual
trace to λ₁ + λ₂ + λ₃ = 5.83e-5. With a reference in-band trace of
6.00e-4 in rx, this sets a ceiling of

10 · log₁₀(6.00e-4 / 5.83e-5) = 10.14 dB

for any single-rank spatial scheme that leaves the TX stage unchanged.


IV. Proposed Canceller

The proposec canceller is  given below

tx_pred   = fit_tx_prediction(rx)                       # stock LS
resid     = rx - tx_pred
resid_bp  = score_filter(resid)                         # per channel
v         = dominant eigenvector of (resid_bp^H · resid_bp / N)
target    = rank-1 projection of resid_bp along v       # 4-column N-vector

# Richardson refinement, α = 0.70, two steps
x = target
for k in {1, 2}:
    x = x + α · (target - score_filter(x))

rx_hat = rx - tx_pred - x


A. Why naive rank-1 subtraction falls short

Setting x = target directly reaches only 7.01 dB. Measured on
channel 0, where:

E[|target[:, 0]|²]                  = 1.91e-5
E[|score_filter(target)[:, 0]|²]    = 1.51e-5

The scorer's score computation applies the 2047-tap Blackman band-pass
FIR a second time to the already-band-filtered target, attenuating
near-edge energy by around 21 % in power. That is almost exactly the 3 dB
gap between naive subtraction (7.01 dB) and the 10.14 dB ceiling that we
derived in Section III.

B. Richardson iteration over the band-pass operator

We want a signal x with H · x ≈ target, where H is the
score_filter operator applied column-wise. Richardson iteration [4]
is a linear fixed-point solver that needs only applications of H:

x_{k+1} = x_k + α · (target - H · x_k)

Starting from x₀ = target, two iterations at α = 0.70 recover almost
all of the attenuated energy. Column-wise application of H preserves
rank-1 structure, where each column is the same time signal scaled by the
spatial eigenvector entry v[c]. So that every iterate stays inside the
rank-1 subspace and the explain_ratio check holds.

C. Hyperparameter choice

The variable α and the step count trade off against the scorer's residual guard.
Observed outcomes on the provided capture:

| steps | α | avg dB | status |
| --- | --- | --- | --- |
| 0 | — | 7.01 | valid |
| 1 | 1.00 | 8.80 | valid |
| 2 | 0.70 | 9.42 | selected |
| 3 | 0.45 | 9.40 | valid |
| ≥ 4 | any | — | INVALID |

(steps = 2, α = 0.70) is an interior point of the feasible region with
the best observed score; (3, 0.45) matches it within 0.02 dB. Variable-
step schedules ([1.0, 0.3], [1.0, 0.5], [0.8, 0.6], [0.7, 0.4, 0.2], ...)
were all dominated.

V. Contribution Breakdown

| Component | Incremental gain |
| --- | --- |
| TX IMD3 LS fit (provided baseline) | 4.02 dB |
| Add rank-1 spatial subtraction (no Richardson) | .01 dB |
| Add Richardson compensation | 9.42 dB |
| Ceiling (ideal rank-1, residual eigen-spectrum) | 10.14 dB |

Richardson is the single largest post-baseline contribution; without it
the scorer's second filter pass silently attenuates the rank-1
subtraction. 9.42 dB sits 0.72 dB under the ceiling; the remaining gap
is Richardson's truncation error plus unavoidable in-band signal/noise
that cannot be removed without exceeding rank-1 and failing validity.


VI. Ablations and Discarded Directions

All entries were implemented and scored on the provided capture.

A. Outer TX-refit iteration (K_outer > 1)

Wrapping the two-stage scheme in an outer alternating minimisation

tx^(k+1) = fit_tx_prediction(rx - e^(k))
e^(k+1)  = rank-1 projection of band_filter(rx - tx^(k+1))

was the initial design plan. In practice K_outer = 2 degrades to
5.27 dB, K_outer = 3 to 5.04 dB, and K_outer = 5 becomes invalid at
0.947 explainability. Each outer pass lets the TX regression recapture
portions of the rank-1 direction, since E[n, c] is not strictly
orthogonal to Φ in time, driving the decomposition out of the
scorer's feasible region. K_outer = 1 is optimal.

B. Enlarging the TX training slice

The stock fit_tx_prediction uses the fixed 200 k-sample window
[20 000, 220 000) for LS. We rebuilt the Gram matrix on up to the
full N = 2,457,600 samples; the gain was nil, with the best variant at
6.91 dB. 200 k / 130 features is already 1500× overdetermined, and
because the scorer's own validity-check regression uses that same
window, any all-N coefficients that disagree with the window-LS
recovery leak into err.

C. FFT-based Wiener deconvolution of H

Substituting Richardson with a direct inverse

x = IFFT( FFT(target) · conj(H) / (|H|² + ε) )

failed (best 3.45 dB, worth than baseline). scipy.signal.convolve(mode="same")
is not exactly an FFT product even after a group-delay shift; boundary
distortions at the 2047-tap filter's edges are amplified by the
inverse. Richardson remains preferable because it only consumes
score_filter as a black-box operator.

D. Richardson on tx_pred

By symmetry with e_hat one could apply Richardson to tx_pred as well.
Tested configurations (steps ∈ {1, 2}, α ∈ {0.5, 0.7, 1.0}) degrade to
4.30–7.51 dB. The reason is that the scorer's own
fit_tx_prediction(rx - rx_hat) regresses using the band-once basis,
so any Richardson-driven deviation of our submitted tx_pred from the
LS-optimum over that same basis flows directly into err. The TX
stage must be emitted exactly as fit_tx_prediction(rx).

E. TX fit with a band² basis

Fitting h against shift(band_filter(band_filter(cubic))) so that the
scorer's second filter pass lands on the LS projection produced all
invalid runs. It is the same failure mode as in VI-D at the coefficient level.


References

[1] M. Z. Waheed, P. P. Campo, D. Korpi, A. Kiayani, L. Anttila, and M. Valkama, 
“Digital cancellation of passive intermodulation in FDD transceivers,” 
in 2018 52nd Asilomar Conference on Signals, Systems, and Computers, 2018.

[2] S. Krikunov, V. Zemlyakov, and A. Ivanov, 
“Physical modelling and cancellation of external passive intermodulation in FDD MIMO,” 
in 2024 IEEE International Multi-Conference on Engineering, Computer and Information Sciences (SIBIRCON), 2024, pp. 102–107.

[3] A. Leshem, A.-J. van der Veen, and A.-J. Boonstra, 
“Multichannel interference mitigation techniques in radio astronomy,” 
Astrophys. J. Suppl. Ser., vol. 131, no. 1, pp. 355–373, 2000.

[4] Y. Saad, Iterative methods for sparse linear systems, 2nd ed. Philadelphia, MS: SIAM, 2003.
