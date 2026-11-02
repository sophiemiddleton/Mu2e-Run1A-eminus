# Cosmics: Cosmic-Ray Background

Cosmic-ray muons that miss the cosmic ray veto can produce an electron in the tracker, or be misreconstructed as one, that looks like a downstream e⁻ near 105 MeV/c. The rate does not follow the beam, so the background is roughly flat in time across the window and smooth in momentum. The veto itself is covered in [`CRV/`](../CRV/).

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/cosmics.py`
- Model definition: `model/physics_components*.py` (look for Cosmic). Momentum: a degree-2 Chebyshev (`poly2`) with constrained `c1` and `c2`. Time: flat
- Normalization: `Cosmic_mix` (`CosmicSignalMix1BB`). The yield is the cut-flow count scaled by the run's on-spill time over the cosmic sample live time, `COSMIC_SAMPLE_LIVETIME_S` = 5.72e6 s in `run_chain.py`. Systematic: `ceHN`, a 20% lnN on the yield
- Control region: the OffSpill cosmic sample, selected with the nominal cuts minus the t0 window, can replace the cosmic momentum shape (`run_chain.py --control-region`, `documentation/ControlRegions.md`). OffSpill timing is not used
- Sideband: the late $e^+$ sideband of the RPC control region (`cosmic_sideband` in the registry) can constrain the cosmic yield

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Plots in pyfitter

`paper_plots.py` group `cosmic`. Run them as in [Plots from pyfitter](../README.md#plots-from-pyfitter).

| Plot | What it shows | Use in the fit |
| :--- | :--- | :--- |
| `cosmic_overlay` | OffSpill and OnSpill cosmic momentum at equal area, with their `poly2` fits and the ratio | Whether OffSpill can stand in for the on-spill shape |
| `cosmic_params` | Fitted (`c1`, `c2`) for OffSpill and OnSpill, with 68% and 95% contours, and the analysis model | The `c1`, `c2` constraints |
| `cosmic_band` | The `poly2` shapes with 68% bands from the fit covariance | Cosmic shape uncertainty |
| `cosmic_time` | OnSpill cosmic t0 after the cuts against the flat time model | Checks the flat time model |
| `cosmic_mixing` | Cosmic momentum with and without beam-background mixing, and OffSpill, each fitted | Effect of pileup on the shape |
| `cosmic_cutflow` | Relative efficiency of each cut for OffSpill and OnSpill, with and without mixing | Where the selections differ |

`component_fits/cosmic_time_shape.py` fits an early-time term to the OnSpill cosmic t0. The partner-track
vetoes are blind for the first ~100 ns of the readout. Nothing in the fit chain reads its output yet.

## Open items

- Cosmic sample live time: `COSMIC_SAMPLE_LIVETIME_S` assumes the full 5.72e6 s sample is processed
- MC statistics: about 40 Cosmic_mix tracks fall in the signal window
- Early-time term in the time model (`cosmic_time_shape.py`)
- Yield systematic: the 20% `ceHN` is a placeholder
