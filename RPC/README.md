# RPC: Radiative Pion Capture Background

Pions in the beam stop in the target and capture: $\pi^- N \rightarrow \gamma N^*$. The photon either converts in material (**external**) or produces an $e^+e^-$ pair directly (**internal** conversion). RPC is prompt and follows the pion lifetime, so the delayed signal time window suppresses it. Pion survival to late times sets the residual rate.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there (`main`, `d9be2e8`):

- Single-component fit: `component_fits/rpc.py`
- Model definition: `model/physics_components.py` (look for RPC). Momentum: a degree-2 Chebyshev (`poly2`) with
  `c1` and `c2` Gaussian-constrained to their MC fit. Time: an exponential (`piexp`) with the slope fixed in the 2D
  fit (`momentum_pdf_builder.DECAY_RATE_RPC`)
- Normalization: `RPCExt_mix` (`RPCExternalPhysicalMix1BB`) = POT × stopped π/POT × (Σ lifetime weights / N resampled)
  × 0.0215 (BR) × the fraction of the Bistirlich spectrum above 50 MeV. `RPCInt_mix` (`RPCInternalPhysicalMix1BB`) is
  that × 0.0069 (internal conversion) (`documentation/DatasetRegistry.md`)
- Control region: early-time $e^+$ constraint on $N_\text{RPC}$, with `run_chain.py --rpc-control`
  (`documentation/ControlRegions.md`). It is off by default in `datasets_1BB.yaml`. `au_run1a_check` uses it
- Systematics: `lumi` (10%) and `rpcsN` (27%, detection efficiency and acceptance). With the control region on,
  its constraint replaces both as `rpcN` on the input card

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Studies in pyfitter

Done in pyfitter at `d9be2e8`. The control-region numbers are from `au_run1a_check`, the Run-1A selection with the
d0 fix and TSDA veto (`<tag>/control/rpc/rpc_control.yaml`). The shape and time fits were made on `au_run1a`,
before the d0 fix.

| Study | Result | Code |
| :--- | :--- | :--- |
| Momentum shape | RPCInt + RPCExt Mix1BB, each weighted to its expected yield, 100–110 MeV/c, 475–1650 ns (12.0k MC events, 35.8 expected): `c1` = −0.2214 ± 0.0156, `c2` = −0.0235 ± 0.0152, χ²/ndf = 18.6/18. Replaces the Gaussian | `paper_plots.py rpc_int_ext` |
| Time slope | Same sample: τ = 17.93 ± 0.16 ns, fixed in the 2D fit | `paper_plots.py rpc_time` |
| Shape systematics | Time: an exponential plus a flat term fitted to the $e^+$ control MC gives τ = 17.6 ± 3.8 ns at the expected data size (e⁺ RPC MC alone: 17.89 ± 0.17 ns), so `rpcTau` shifts the slope by 0.012 /ns. Momentum: internal (`c1` = −0.255) and external (`c1` = −0.182) RPC have different shapes, and `rpcShape` moves (`c1`, `c2`) by half their difference, a stand-in for the unknown conversion mix | `component_fits/rpc_shape_systs.py`, output `component_fits/shapes/rpc_shape_systs.yaml` |
| $e^+$ control region | $e^+$ with 100 < p < 110 MeV/c (above the RMC endpoint) and 475 < t < 540 ns, about 96% RPC. Fitted $N_\text{RPC}(e^+)$ = 29.1 ± 5.5 (29.1 expected). The cosmics are constrained from the 650–1650 ns sideband (0.86 ± 0.24). Transfer factor 1.169 ± 1.4%. Constraint: $N_\text{RPC}$ = 34.0 ± 6.4 in the fit window, `rpcN` = 1.189, almost all from the control fit statistics | `control_region.py`, `run_chain.py --rpc-control` |
| $e^+$/$e^-$ charge asymmetry | Flat e⁻/e⁺ Mix1BB samples. In the RPC window $R$ = 1.0097 ± 0.0081, so `charge_syst` = 0.0097 in `datasets_1BB.yaml` | `charge_asymmetry.py`, outputs in `flat_charge/`, `documentation/Studies.md` |
| RPC in the signal window | 0.132 events in the cut-and-count window (103.34–104.74 MeV/c, 540–1650 ns), from 34 RPCInt and 9 RPCExt MC tracks | `cutcount/count_table.csv` |

### Plots in pyfitter

Run them as in [Plots from pyfitter](../README.md#plots-from-pyfitter).

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `rpc_int_ext` | Internal and external RPC e⁻ momentum in the fit window, each scaled to its yield and fitted with the model's Chebyshev | Source of the RPC momentum coefficients in `model/physics_components.py` |
| `rpc_time` | Internal and external RPC e⁻ t0, each with an exponential fit, against the RPC time model | Source of the RPC time slope |
| `pion_time` | RPC e⁻ t0 with an exponential fit, against the RPC time-model rate | Checks the RPC time model |

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

- Put `rpcTau` and `rpcShape` on the input card: `rpc_shape_systs.py` writes them, and the limit
  (`simple_limit_combine_12d.py`) reads `param` systematics, but `run_chain.py` does not add them
- Refit the momentum shape, time slope and shape systematics on the Run-1A selection with the d0 fix
  (`au_run1a_check`). They were made on `au_run1a`
- Out-of-time RPC from protons between pulses (extinction): flat in time, not in the RPC MC, and not
  constrained by the control region
- Pion survival to late times: the lifetime weights set the rate in the signal window, with no systematic yet
- Transfer factor: it assumes the simulated internal/external mix. The early-to-late time extrapolation has no
  systematic in mu⁻ → e⁻. The mu⁻ → e⁺ region uses `transfer_syst` = 0.33
- `rpcsN`: `uncertainties_config.py` has 1.27 with a comment saying 15%
- Make the $e^+$ control region the default in `datasets_1BB.yaml`, or record which Run-1A results use it
- Charge asymmetry at 2BB: the 2BB registries have no RPC section
