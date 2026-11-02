# RMC: Radiative Muon Capture Background

In radiative muon capture, $\mu^- \text{Al} \rightarrow \nu_\mu \gamma\, \text{Mg}^*$, a photon converts to $e^+e^-$, externally or internally. The photon endpoint is below the CE energy, so RMC only enters the signal region through resolution tails. Its time follows the muonic-aluminum lifetime. Conversions give $e^+$ and $e^-$ in equal numbers, and DIO gives no $e^+$, so late-time $e^+$ can constrain RMC.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there (`main`, `d9be2e8`):

- Single-component fit: `component_fits/rmc.py`
- Model definition: `model/physics_components.py` (look for RMC). Momentum: the gamma-poly hybrid
  (`custom_models.rmc_gamma_poly_hybrid`) with $\alpha = \lambda = 0$, which leaves the endpoint power law
  $(k_\text{max} - p)^\beta$. The shape is fixed and $N_\text{RMC}$ floats. The $e^+$ control region has its own
  shape (`eplus`). Time: the CE/DIO exponential with a fixed slope (`momentum_pdf_builder.py`)
- Normalization: `RMCExt_mix` (`RMCPhaseSpace0NExternalMix1BB`) = $N_{\text{stop}\,\mu}$ × 0.609 × 1.41e-5 (RMC above
  57 MeV) × 0.099 (0N fraction) × the Plestid 0N fraction above 80 MeV relative to 57 MeV. `RMCInt_mix` is that
  × 0.0069 (internal/external). The 1N samples give nothing above the 80 MeV generator threshold
  (`documentation/DatasetRegistry.md`)
- Control region: late-time $e^+$ constraint on $N_\text{RMC}$, with `run_chain.py --rmc-control`
  (`documentation/ControlRegions.md`). It is off by default in `datasets_1BB.yaml`. `au_run1a_check` uses it
- Systematics: `lumi` (10%). With the control region on, its constraint replaces the yield
  uncertainty as `rmcN` on the input card

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Studies in pyfitter

Done in pyfitter at `d9be2e8`. The control-region numbers are from `au_run1a_check`, the Run-1A selection with the
d0 fix and TSDA veto (`<tag>/control/rmc/rmc_control.yaml`).

| Study | Result | Code |
| :--- | :--- | :--- |
| Momentum shape, e⁻ | RMCInt + RMCExt Mix1BB, 90–110 MeV/c, 475–1650 ns (46.7k MC events): $\beta$ = 2.826 ± 0.048, $k_\text{max}$ = 101.26 ± 0.10 MeV/c, χ²/ndf = 28.1/20. Replaces a Crystal Ball (χ² = 898/10) | `component_fits/rmc_shape.py`, output `component_fits/shapes/rmc_eminus_shape.yaml` |
| Momentum shape, e⁺ | Same form on the $e^+$ RMC MC, 95–110 MeV/c, 540–1650 ns (1361 MC events): $\beta$ = 2.92 ± 0.30, $k_\text{max}$ = 100.09 ± 0.27 MeV/c, χ²/ndf = 5.3/8 | `component_fits/rmc_shape.py`, output `component_fits/shapes/rmc_eplus_shape.yaml` |
| $e^+$ control region | $e^+$ with 95 < p < 110 MeV/c and 540 < t < 1650 ns. Fitted $N_\text{RMC}(e^+)$ = 129.5 ± 11.8 (129.1 expected), next to cosmics (23.5 ± 17.7) and RPC (≈ 0). Transfer factor 0.0122 ± 25% (16 RMC MC tracks in the fit window). Constraint: $N_\text{RMC}$ = 1.58 ± 0.42 in the fit window, `rmcN` = 1.269 | `control_region.py`, `run_chain.py --rmc-control` |
| $e^+$/$e^-$ charge asymmetry | Flat e⁻/e⁺ Mix1BB samples (7.5M each, 70–110 MeV/c). In the RMC window $R = \varepsilon(e^-)/\varepsilon(e^+)$ = 1.0277 ± 0.0021, so `charge_syst` = 0.0277 in `datasets_1BB.yaml`. It is conservative, because the transfer factor already takes $R$ from the MC | `charge_asymmetry.py`, outputs in `flat_charge/`, `documentation/Studies.md` |
| RMC in the signal window | $k_\text{max}$ is below the signal window, so RMC only reaches it through resolution. 1.57 events expected in the fit window (100–110 MeV/c, 475–1650 ns). 0 in the cut-and-count window (103.34–104.74 MeV/c, 540–1650 ns): 16 RMCExt MC tracks pass the cuts, none in the window | `au_run1a_check/input_card.yaml`, `cutcount/count_table.csv` |

### Plots in pyfitter

Run them as in [Plots from pyfitter](../README.md#plots-from-pyfitter). † marks a plot that is newer than the pinned pyfitter.

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `rmc_shape` | RMC momentum shape used in the fits (e⁻ and e⁺), over the MC it was fitted to (needs `--shapes-dir`) | Checks the fit from `component_fits/rmc_shape.py` |
| `rmc_int_ext` † | Internal and external RMC momentum, each with its own shape fit, over their sum | Whether one shape describes both |
| `rmc_time` † | RMC e⁻ t0 against the 2D fit's RMC time model (the fixed CE/DIO slope) | Checks the RMC time model |

The $e^+$ control-region plots (`eplus` group) cover RPC and RMC together:

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `eplus_fits` | RMC and RPC $e^+$ control fits, with the fitted components and pulls | The $N_\text{RMC}$ and $N_\text{RPC}$ constraints (`rmcN`, `rpcN`) |
| `eplus_time` | $e^+$ t0 by process, with the RPC and RMC windows and the cosmic sideband | Choice of the control windows |
| `eplus_budget` | Uncertainty of each constraint, by source | Size of `rmcN` and `rpcN` |
| `eplus_closure` | Fitted / expected yield of each component in the control fits | MC closure of the control fits |
| `charge_ratio` | $\varepsilon(e^-)/\varepsilon(e^+)$ per true-momentum bin from the flat samples (needs `--flat-dir`) | `charge_syst` = $|1 - R|$ in the registry |
| `charge_cuts` | e⁻/e⁺ ratio of each cut's efficiency (needs `--flat-dir`) | Which cuts are charge asymmetric |

## Open items

- Transfer-factor MC statistics: only 16 RMC MC tracks fall in the fit window, which gives 25% of the 27% `rmcN`.
  More RMC MC near the endpoint would shrink the constraint
- RMC rate uncertainty without the control region: the mu⁻ → e⁻ card has only `lumi` for RMC. The mu⁻ → e⁺
  card uses `rmcN` = 7.9% (rate 1.40 ± 0.11, TRIUMF 1999)
- Make the $e^+$ control region the default in `datasets_1BB.yaml`, or record which Run-1A results use it
- Charge asymmetry at 2BB: the 2BB registries have no RMC section, and would need a 2BB flat sample
- Generator-level spectrum: the 0N Plestid model above 80 MeV sets both the normalization and the MC shape. 1N
  is left out because it gives nothing above 80 MeV in mu⁻ → e⁻
