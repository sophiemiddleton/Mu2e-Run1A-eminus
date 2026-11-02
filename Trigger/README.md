# Trigger: Online Track Trigger

Events are kept online by the track triggers. The analysis requires either the calorimeter-seeded (CPR) or the tracker-seeded (APR) downstream-electron path. The trigger efficiency falls with pileup, and the beam flash is concentrated early in the microbunch, so the loss can depend on both intensity and t0.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Cut: `or_trigger` in `analyze.py`, `trig_cpr_TrkDe_80m70p` OR `trig_apr_TrkDe_80m70p`
- Efficiency: taken from the MC, inside the CE selection efficiency. There is no separate trigger systematic yet

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Plots in pyfitter

The trigger plots of `paper_plots.py` group `signal`. Run them as in
[Plots from pyfitter](../README.md#plots-from-pyfitter). † marks a plot that is newer than the pinned pyfitter.
All of them need `--trigdiag-dir`: a tag made with `process.py --cut-profile trigdiag`, which applies
every nominal cut except the trigger and stores each trigger decision.

| Plot | What it shows |
| :--- | :--- |
| `trigger_eff` † | Absolute CE efficiency with each trigger path against protons per pulse, `run1a` selection (t0 > 475 ns) |
| `trigger_eff_run1a_1D` † | As `trigger_eff`, in the 1D-fit window (t0 > 540 ns) |
| `trigger_eff_time` † | Trigger efficiency of offline-selected CE against protons per pulse, in t0 chunks, per path. A slope that steepens at early t0 is the flash |
| `trigger_eff_flash_scale` † | Change of the trigger efficiency with the beam flash scaled by k = 1 to 3, per path and t0 chunk |

## Open items

- Trigger efficiency systematic for the fit
- Dependence on the beam-flash level
