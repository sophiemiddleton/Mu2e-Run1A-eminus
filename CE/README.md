# CE: Conversion Electron Signal

The signal is a mono-energetic electron at about 104.97 MeV/c from $\mu^- \text{Al} \rightarrow e^- \text{Al}$. After energy loss in the target and detector, and after tracker resolution, the reconstructed momentum peak is shifted and has a low-side tail. Its time follows the muonic-aluminum lifetime (about 864 ns).

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/ce.py`
- Model definition: `model/physics_components*.py` (look for CE)
- Normalization: `CeM_mix` in the dataset registry. Expected yield = $N_{\text{stop}\,\mu}$ × captures per stopped muon × $R_{\mu e}$
- Background: `documentation/LikelihoodDefinition.md`

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Plots in pyfitter

`paper_plots.py` groups `signal` and `efficiency`, and the CE plots of `timing`. Run them as in [Plots from pyfitter](../README.md#plots-from-pyfitter). † marks a plot that is newer than the pinned pyfitter.

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `signal_mom` † | CE momentum in the fit window, with the fit-model DSCB and a free DSCB refit | Checks the CE momentum model |
| `signal_response` † | Momentum response $p_\text{reco} - p_\text{true}$ of flat e⁻ in slices of $p_\text{true}$, against CE | Separates the response from the radiative tail |
| `signal_time_slices` † | CE momentum shape in t0 slices | Checks that the 2D fit can factorize momentum and time |
| `signal_pileup` † | CE momentum shape in terciles of protons per pulse | Pileup dependence of the shape |
| `signal_eff_intensity` † | CE selection efficiency against protons per pulse | Pileup dependence of the efficiency, and the beam reweighting |
| `efficiency` | Flat e⁻ selection efficiency per true-momentum bin in the fit window (needs `--flat-dir`) | Acceptance × efficiency against momentum |
| `muon_time` † | CE t0 with an exponential fit, against the Muon time-model rate | Checks the CE time model |
| `muon_time_truth` † | CE true origin time before and after the cuts, and the efficiency against true time | Whether the cuts change the selected lifetime |
| `muon_time_cut_eff` † | Efficiency of each cut against CE true time, and its slope | Which cuts lengthen the selected lifetime |

The truth plots need `<sample>_truth_precut.parquet` from `process.py`.

## Open items

- Momentum response function (energy loss ⊗ resolution) and its parameterization
- Momentum scale and resolution systematics
- Acceptance × efficiency as a function of momentum and time
