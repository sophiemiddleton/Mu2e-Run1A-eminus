# Cosmics: Cosmic-Ray Background

Cosmic-ray muons that miss the cosmic ray veto can produce an electron in the tracker, or be misreconstructed as one, that looks like a downstream e⁻ near 105 MeV/c. The rate does not follow the beam, so the background is roughly flat in time across the window and smooth in momentum. The veto itself is covered in [`CRV/`](../CRV/).

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there (`main`, `d9be2e8`):

- Single-component fit: `component_fits/cosmics.py`
- Model definition: `model/physics_components.py` (look for Cosmic). Momentum: a degree-2 Chebyshev (`poly2`),
  with `c1` = 0.020 and `c2` = −0.082 Gaussian-constrained. Time: flat (`uniform`)
- Normalization: `Cosmic_mix` (`CosmicSignalMix1BB`). The yield is the cut-flow count scaled by the run's on-spill
  time over the cosmic sample live time, `COSMIC_SAMPLE_LIVETIME_S` = 5.72e6 s in `run_chain.py`. Systematic:
  `ceHN`, a 20% lnN on the yield
- Control region: the OffSpill cosmic sample (`OffSpill` in the registry), selected with the nominal cuts minus
  the t0 window. Its `poly2` fit replaces the cosmic momentum shape, fixed or constrained
  (`run_chain.py --control-region --control-mode fixed|constrain`, `documentation/ControlRegions.md`). It is off
  by default in `datasets_1BB.yaml`. `au_run1a_check` uses it, with the shape fixed. OffSpill timing is not used,
  and the yield still comes from the on-spill MC
- Sideband: the late $e^+$ sideband of the RPC control region (650–1650 ns, `cosmic_sideband` in the registry)
  constrains the cosmics inside the RPC control fit

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Studies in pyfitter

Done in pyfitter at `d9be2e8`. The yields and the OffSpill fit are from `au_run1a_check`, the Run-1A selection with the
d0 fix and TSDA veto (`<tag>/control/cosmic_shape.yaml`). The time-shape and CRV diagnoses were made on `au_run1a`,
before the d0 fix.

| Study | Result | Code |
| :--- | :--- | :--- |
| OffSpill momentum shape | 377 OffSpill tracks, 100–110 MeV/c: `c1` = −0.123 ± 0.089, `c2` = −0.010 ± 0.087. `au_run1a_check` fixes the cosmic shape to these | `control_region.py`, `run_chain.py --control-region` |
| Early-time excess | The OnSpill cosmic t0 has an excess at early times. A flat-plus-exponential fit to 391 `Cosmic_mix` tracks (475–1650 ns) gives an early fraction $f_e$ = 0.054 ± 0.017 with $\tau_e$ = 23 ± 8 ns (about 21 tracks), preferred over flat (likelihood ratio 18.8, p ≈ 8e-5). Test only: nothing in the fit chain reads it | `component_fits/cosmic_time_shape.py`, output `component_fits/shapes/cosmic_time_shape_test.yaml` |
| Cause of the early excess | The partner-track vetoes (upstream, multi-track, CRV) are blind for the first ~100 ns of the readout. Selections without the veto and t0 cuts (`--cut-profile crvdiag`, `vetodiag`) give the vetoed fraction against t0 and the earliest CRV coincidence time | `cosmic_crv_diagnosis.py`, outputs in `au_run1a_crvdiag/crv_diagnosis/`, `au_run1a_vetodiag/crv_diagnosis/` |
| Cosmics in the windows | 42.1 events expected in the fit window (100–110 MeV/c, 475–1650 ns). 5.33 in the cut-and-count window (103.34–104.74 MeV/c, 540–1650 ns), from 39 MC tracks: the largest background there | `au_run1a_check/input_card.yaml`, `cutcount/count_table.csv` |

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

## Open items

- Cosmic sample live time: `COSMIC_SAMPLE_LIVETIME_S` assumes the full 5.72e6 s sample is processed (marked FIXME)
- Constraint widths of the default shape: the Cosmic `pars` look like (value, lower, upper) with ±0.061, but
  `constrain` reads them as (value, σ_lo, σ_hi) (`pdf_builder._create_constrained_param`), so σ is 0.081 for `c1`
  and 0.142 for `c2`. It only matters when the OffSpill control region is off
- Early-time term in the time model: the fit uses a flat time model. Add the `cosmic_time_shape.py` term, or move
  the t0 window past the blind period. Redo the test on `au_run1a_check`
- MC statistics: 39 `Cosmic_mix` tracks in the cut-and-count window give its largest background
- Yield systematic: the 20% `ceHN` is a placeholder
- Control-region mode: compare no control region, a fixed and a constrained OffSpill shape
  (`documentation/Studies.md`, "Cosmic control region", proposed)
