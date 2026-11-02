# StoppingTarget: CE and DIO for the target geometry options

Compares the CE signal and the DIO background after the Run-1A cuts (`run1a`, which includes
the TSDA veto) for the nominal stopping target and three alternative geometries. Each is scaled to a nominal run: 28 days of 1BB at 3.84 kW.
The study gives the signal acceptance and DIO yield in the signal window for each
geometry. It doesn't feed the main fit directly.

## Datasets

All are MDC2025ba, unmixed on-spill, reco `best_v1_5`, 200 files each:

| Geometry | CE | DIO |
| :--- | :--- | :--- |
| `nominal` | `nts.mu2e.CeEndpointOnSpill.MDC2025ba_best_v1_5.root` | `nts.mu2e.DIOtail95OnSpill.MDC2025ba_best_v1_5.root` |
| `st37plus6up6dn` | `nts.mu2e.CeEndpointOnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root` | `nts.mu2e.DIOtail95OnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root` |
| `st49hole1to5` | `nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole1to5_best_v1_5.root` | `nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole1to5_best_v1_5.root` |
| `st49hole9to41` | `nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole9to41_best_v1_5.root` | `nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole9to41_best_v1_5.root` |

From the stop positions (`plot_stops_st.py`, `stops_rz`):

- **`st37plus6up6dn`**: 6 upstream foils (r 20–112 mm), the 37 nominal foils (r 21.6–75 mm), and
  6 downstream foils (r 2–110 mm).
- **`st49hole1to5`**: 49 foils with central holes of 1–5 mm radius. The outer radius runs from
  80 mm to 112 mm and back to 96 mm.
- **`st49hole9to41`**: 49 foils with the hole radius growing from 9 mm upstream to 41 mm
  downstream, and the outer radius from 100 mm to 118 mm.

## Normalization

`make_registry_st.py` uses pyfitter's `make_registry.py` and `normalization.py`:

- **POT.** 28 days × 0.323 duty factor gives 7.81e5 s on spill. The 1BB intensity is
  scaled from its quoted 3.8 kW to 3.84 kW, which gives 7.36e18 POT.
- **Stopped muons.** The nominal-target SimEfficiencies2 rate (Sim_best v1_1): 5.65e15 stopped μ⁻,
  scaled by each geometry's `stop_rate_rel` (below).
- **CE.** Stopped μ × 0.609 captures × R_μe = 1e-15, which gives 3.44 events
  generated. The sample is `CeEndpoint`: mono-energetic, with no radiative tail.
- **DIO.** Stopped μ × DIO fraction × the fraction of the Czarnecki/Szafron spectrum
  above 95 MeV, which gives 8.03e4 events generated.
- **Acceptance.** mcs events / generated events, from metacat.

Each track after the cuts gets the weight `expected × acceptance / N_generated`, as in
`run_chain.py`.

The muon-stop rate depends on the geometry. The registry records it as `stop_rate_rel`, stopped μ⁻
per POT relative to the nominal target:

- **Nominal.** SimEfficiencies2, as Production `normalizations.py`: `MuBeamCat` × `MuminusStopsCat` × 1000
  = 7.67e-4 per POT. The 1000 undoes the stops stage's `MuminusTargetStopPrescale`.
- **Geometries.** Their stops stage (`cnf.mu2e.TargetStops.MDC2025ba_<geometry>.tar`) resampled the same
  `MuBeamCat.MDC2025ab` with `TargetStopPrescaleFilter.nPrescale: 1`, so stops per POT are
  `MuBeamCat` × `MuminusStopsCat` with no 1000.

| Geometry | μ⁻ stops per `MuBeamCat` event | Stops per POT | `stop_rate_rel` |
| :--- | ---: | ---: | ---: |
| `nominal` | 0.359 | 7.67e-4 | 1 |
| `st37plus6up6dn` | 0.586 | 1.25e-3 | 1.63 |
| `st49hole1to5` | 0.653 | 1.40e-3 | 1.82 |
| `st49hole9to41` | 0.583 | 1.25e-3 | 1.62 |

The nominal rate is from the MDC2025ac stops stage and the geometries' from MDC2025ba, on the same
`MuBeamCat` input. The sensitivity table always applies `stop_rate_rel`; `plot_st.py --scale-stops` also
applies it to the plots.
The ratio plots are to the first geometry, `nominal`.

## Sensitivity

`plot_st.py` writes `sensitivity_st.csv`. The table below is from the `run1a` cuts (signed d0 and the
TSDA veto), the signal window 103.4 < p < 104.8 MeV/c and 540 < t < 1650 ns, for 28 days at 3.84 kW:

| Geometry | N<sub>stops</sub> | N<sub>CE</sub> (R<sub>μe</sub> = 1e-15) | N<sub>DIO</sub> | N<sub>cosmic</sub> | CE reco eff | CE sel eff | SES | Limit (events) | R<sub>μe</sub> (90% CL) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `nominal` | 5.65e15 | 0.464 | 0.010 ± 0.010 (1 track) | 5.19 ± 0.84 | 0.392 | 0.344 | 2.15e-15 | 5.00 | 1.08e-14 |
| `st37plus6up6dn` | 9.22e15 | 0.458 | 0 (0 tracks) | 5.19 ± 0.84 | 0.367 | 0.223 | 2.18e-15 | 4.99 | 1.09e-14 |
| `st49hole1to5` | 1.03e16 | 0.434 | 0 (0 tracks) | 5.19 ± 0.84 | 0.364 | 0.191 | 2.31e-15 | 4.99 | 1.15e-14 |
| `st49hole9to41` | 9.17e15 | 0.429 | 0.016 ± 0.016 (1 track) | 5.19 ± 0.84 | 0.363 | 0.212 | 2.33e-15 | 5.00 | 1.16e-14 |

The limit also includes RPC, 0.128 events in this window, the same for every geometry (`N_RPC` in the csv). It was made with the
cosmics and RPC of `au_run1a_tsda` (the main analysis rerun with these cuts) and the systematics of the
older `au_run1a` card, since `au_run1a_tsda` has no count card yet. 

```bash
python StoppingTarget/plot_st.py --tag $OUT/run1a_tsda --main-tag $PYFITTER_DIR/au_run1a_tsda \
    --card $PYFITTER_DIR/au_run1a/cutcount/count_card.yaml --p-window 103.4 104.8 --t-window 540 1650 \
    --optimize-p 102.0 103.8 104.2 105.2 --scale-stops
```

`--p-window` and `--t-window` set the signal window (default: the cut set's `signal_region`).

### Momentum window per geometry

The extra material moves the CE peak down (median p<sub>reco</sub> − p<sub>gen</sub> −1.3 to −1.5 MeV/c, against
−0.9 for nominal), so the fixed window loses CE in the alternatives. `--optimize-p LO_MIN LO_MAX HI_MIN HI_MAX`
scans each geometry's window edges (`--optimize-step`, default 0.1 MeV/c) for the lowest expected limit and
writes `sensitivity_st_opt.csv`. The scan takes the cosmics as flat in momentum (the rate per MeV/c over
100–110 MeV/c, 3.7 per MeV/c here), so it doesn't follow the fluctuations of the 38 MC tracks in the window:

| Geometry | Window (MeV/c) | N<sub>CE</sub> | N<sub>DIO</sub> | N<sub>cosmic</sub> | N<sub>RPC</sub> | CE sel eff | SES | Limit (events) | R<sub>μe</sub> (90% CL) |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `nominal` | 103.4–104.8 | 0.464 | 0.010 (1 track) | 5.20 | 0.128 | 0.344 | 2.15e-15 | 5.00 | 1.08e-14 |
| `st37plus6up6dn` | 102.7–104.8 | 0.598 | 0.066 (4 tracks) | 7.80 | 0.167 | 0.290 | 1.67e-15 | 6.05 | 1.01e-14 (−6%) |
| `st49hole1to5` | 102.5–104.7 | 0.646 | 0.091 (5 tracks) | 8.18 | 0.187 | 0.284 | 1.55e-15 | 6.19 | 0.96e-14 (−11%) |
| `st49hole9to41` | 102.6–104.7 | 0.586 | 0.065 (4 tracks) | 7.80 | 0.183 | 0.289 | 1.71e-15 | 6.05 | 1.03e-14 (−4%) |

The fixed window is the best one for `nominal`. With their own windows the alternatives do 4–11% better than
nominal, as the Set 79 v14 comparison (cut-and-count 10.68, 10.01, 9.49, 10.24 × 1e-15). The windows are
chosen on the same MC they are evaluated on, so the gain is slightly optimistic, and the DIO rests on 4–5 tracks.

One shared window, 102.7 < p < 104.8 MeV/c, comes within 0–2% of each alternative's own best limit, and costs
nominal 6% (`window_scan_plo`, `window_scan_2d`). `ce_eff_compare` shows where the alternatives lose CE: about 7%
in acceptance and reconstruction, then `within_d0` (5–9%), the 100–110 MeV/c range, and the fixed window
(down to 0.52–0.60 of nominal). Their own windows bring that back to 0.77–0.79 of nominal.

- **N<sub>stops</sub>.** The run's stopped μ⁻ × `stop_rate_rel`. The table always scales the yields by
  the stop rate, with or without `--scale-stops`.
- **N<sub>cosmic</sub>.** The same for every geometry: the main-analysis `Cosmic_mix` sample
  (CosmicSignalMix1BB, nominal target) after the Run-1A cuts, counted in the same signal window
  (38 tracks × 0.137, on-spill time / livetime). `--main-tag` sets the main-analysis tag
  (default `$PYFITTER_DIR/au_run1a`). The cosmics are also stacked in the `mom_` and `time_`
  plots. That tag starts at 100 MeV/c, so there are no cosmics drawn below that.
- **CE reco eff.** A reconstructed track per generated CE: acceptance × `has_a_track` / all events.
- **CE sel eff.** CE tracks in the signal window after all the cuts, per reconstructed track. Reco ×
  selection is the total CE efficiency.
- **SES.** 1 / (N<sub>stops</sub> × 0.609 × total CE efficiency), which is 1e-15 / N<sub>CE</sub>.
- **R<sub>μe</sub>.** The Combine AsymptoticLimits expected (Asimov) limit from
  `sensitivity/cutcount/counting_stats.py`. It runs on the main-analysis `cutcount/count_card.yaml`
  (`--card` to use another), with CE, DIO and cosmics replaced by the geometry's yields, and RPC and
  RMC recounted from the main-analysis tag in the window: the same processes and lnN systematics
  (lumi, ceHN, dioN, rpcsN).

All four samples are unmixed best_v1_5, so compare the geometries with `nominal`, not with the
main analysis (mixed MDC2025au). `st37plus6up6dn` is the nominal 37 foils plus 6 upstream and 6
downstream: the extra foils give more energy loss, which moves CE out of the momentum window.

The DIO in the window is limited by MC statistics: one track weighs about 0.01 events.

## Stops and the CE response

`plot_st.py` also plots the CE momentum response, p<sub>reco</sub> − p<sub>gen</sub>. The generated
momentum is the production momentum, so the response includes the energy loss in the target as
well as the tracker resolution.

`dump_stops_st.py` runs Offline's `StoppedParticlesDumper` (`dump_stops.fcl`) on the
`MuminusStopsCat` files. uproot can't read `SimParticle` from art files. The stops come out in the
detector frame (Mu2e x + 3904 mm, z − 10171 mm), which is the frame of the EventNtuple truth.
The nominal stops come from `MuminusStopsCat.MDC2025ac` (the SimEfficiencies2 sample, prescale 1000) and the
geometries' from MDC2025ba (prescale 1). The plots are per `MuBeamCat` event (stops × prescale / `gen.count`),
the same input for all four: 0.359 (nominal), 0.587, 0.654 and 0.584 (+63%, +82%, +63%). These count stops, not
events, so they differ from `stop_rate_rel` by 0.2%.

The CE electron starts where its muon stopped. `process.py` writes that position to the parquet
(`pos_x_true`, `pos_y_true`, `pos_z_true`). `plot_stops_st.py` uses it to give the CE efficiency and
signal for each foil.

## Running

```bash
source setup.sh
getToken
OUT=/exp/mu2e/data/users/$USER/stoppingtarget

# 1. Registry and file lists
python StoppingTarget/make_registry_st.py --outdir $OUT

# 2. Run-1A selection with the TSDA veto (95-110 MeV/c, so the DIO tail shows)
cd $PYFITTER_DIR
python run_chain.py --registry $OUT/datasets_st.yaml --tag $OUT/run1a_tsda --version run1a --group all \
    --skip-fit --loc tape --jobs 3 --fitrange_low 95 475 --fitrange_hi 110 1650
cd -

# 3. Plots and the sensitivity table
python StoppingTarget/plot_st.py --tag $OUT/run1a_tsda
```

```bash
# 4. Muon stops (art dump of MuminusStopsCat; ~2 min and ~1 GB per file, runs on a gpvm)
python StoppingTarget/dump_stops_st.py --outdir $OUT --max-files 3 --jobs 3   # nominal: one 7 GB file, ~15 min

# 5. Stops per foil, and CE per foil (the per-foil CE needs pos_z_true: run_chain with the current process.py)
python StoppingTarget/plot_stops_st.py --outdir $OUT --tag $OUT/run1a_tsda --p-window 103.4 104.8 --t-window 540 1650
```

`make_registry_st.py` takes `--days`, `--beam-power` and `--rmue` to change the run.
Keep `--jobs` at about one worker per 3 GB of memory. On a 4-core, 11 GB node, 16
workers ran out of memory.

## Outputs

Everything is under `/exp/mu2e/data/users/sophie/stoppingtarget/`:

| Path | Contents |
| :--- | :--- |
| `datasets_st.yaml` | Registry: acceptances, expected yields, run normalization, `stop_rate_rel` |
| `run1a_tsda/manifest.yaml`, `run1a_tsda/cutflows/` | Selection output from `run_chain.py` |
| `run1a_tsda/plots/mom_<geometry>.{png,pdf}` | Momentum in the signal time window, DIO and cosmics stacked, CE line |
| `run1a_tsda/plots/time_<geometry>.{png,pdf}` | Time in the signal momentum window, DIO and cosmics stacked, CE line |
| `run1a_tsda/plots/mom_compare.{png,pdf}` | All geometries overlaid, ratio to `nominal` |
| `run1a_tsda/plots/mom_<geometry>_opt`, `time_<geometry>_opt` | As above, in the geometry's optimized momentum window (`--optimize-p`) |
| `run1a_tsda/plots/ce_eff_compare.{png,pdf,csv}` | CE efficiency per generated event after each cut, then in the fixed and own windows (540–1650 ns), ratio to `nominal` |
| `run1a_tsda/plots/window_scan_plo.{png,pdf}` | Limit against the lower window edge, each geometry at its best upper edge (`--optimize-p`) |
| `run1a_tsda/plots/window_scan_2d.{png,pdf}`, `window_scan.csv` | Limit over both window edges per geometry, optimum marked; the csv has every window's yields (`--optimize-p`) |
| `run1a_tsda/plots/window_yields.{png,pdf}` | CE and backgrounds against the lower edge, at the best shared upper edge (`--optimize-p`) |
| `run1a_tsda/plots/sensitivity_st_opt.csv` | The sensitivity table in each geometry's optimized momentum window (`--optimize-p`) |
| `run1a_tsda/plots/sensitivity_st.csv` | Per geometry: stops, CE, DIO and cosmics in the signal window, CE reco and selection efficiency, SES, R<sub>μe</sub> limit |
| `run1a_tsda/plots/resolution_compare.{png,pdf}` | CE p<sub>reco</sub> − p<sub>gen</sub>, unit area, signal window marked, ratio to `nominal` |
| `run1a_tsda/plots/resolution_st.csv` | CE response: peak, median, FWHM, core RMS, fraction in the signal window |
| `stops/stops_<geometry>.parquet`, `stops/stops_summary.yaml` | Stop x, y, z (detector frame), t; stops per generated event |
| `stops/plots/stops_z`, `stops_r`, `stops_rz` | Stops per foil against z, stop radius, r against z, per 10<sup>3</sup> `MuBeamCat`; ratios to `nominal` |
| `stops/plots/stops_foils.csv` | Per foil: z, inner and outer stop radius, stops per generated event |
| `stops/plots/ce_foils.{png,pdf,csv}` | Per foil: CE signal-window efficiency, median response, CE signal for the run |

## Contact

Sophie Middleton (sophie@fnal.gov)
