# RPC: Radiative Pion Capture Background

Pions in the beam stop in the target and capture: $\pi^- N \rightarrow \gamma N^*$. The photon either converts in material (**external**) or produces an $e^+e^-$ pair directly (**internal** conversion). RPC is prompt and follows the pion lifetime, so the delayed signal time window suppresses it. Pion survival to late times sets the residual rate.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/rpc.py`
- Normalization: `RPCExt_mix` and `RPCInt_mix`. Expected yield = POT × stopped π/POT × survival weights × RPC branching fraction × the fraction of the Bistirlich spectrum above 50 MeV, times 0.0069 for internal
- Control region: $e^+$ fit in `control_region.py` (`documentation/ControlRegions.md`)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Plots in pyfitter

Run them as in [Plots from pyfitter](../README.md#plots-from-pyfitter). † marks a plot that is newer than the pinned pyfitter.

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `rpc_int_ext` † | Internal and external RPC e⁻ momentum in the fit window, each scaled to its yield and fitted with the model's Chebyshev | Source of the RPC momentum coefficients in `model/physics_components.py` |
| `rpc_time` † | Internal and external RPC e⁻ t0, each with an exponential fit, against the RPC time model | Source of the RPC time slope |
| `pion_time` † | RPC e⁻ t0 with an exponential fit, against the RPC time-model rate | Checks the RPC time model |

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

- Momentum and time shape: replace the Gaussian with a data-driven model
- Pion survival and extinction systematics
- Use of the $e^+$ control region
