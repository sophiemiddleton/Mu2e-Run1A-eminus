# Sensitivity: Running pyfitter

The sensitivity (90% CL upper limit on $R_{\mu e}$) comes from
[pyfitter](https://github.com/Mu2e/pyfitter). One script, `run_chain.py`, runs
the whole chain:

1. **Selection.** Apply the analysis cuts to each dataset.
2. **Fit.** Build the per-process PDFs and fit them.
3. **Limit.** Run the CLs scan and compute the upper limit.

For full details see `documentation/RunChain.md` in pyfitter.

## 1. Set up

From the top of this repository, in every new shell:

```bash
source setup.sh
getToken            # only needed to read files from tape or disk
cd $PYFITTER_DIR    # for run_chain.py
```

`setup.sh` does the following:

1. Runs `mu2einit` if you haven't already.
2. Activates `pyenv ana 2.8.0`.
3. Checks out the `pyfitter` submodule at the top of this repository, at the commit
   this repository pins, and warns if the checkout has moved off that commit.
4. Exports `PYFITTER_DIR` and `RUN1A_DIR` and checks that every package and pyfitter
   file the scripts need is present.

Mu2e/pyfitter is private, so the checkout uses ssh and needs a GitHub ssh key: test it
with `ssh -T git@github.com`. To change the defaults, set these before sourcing:

| Variable | Default | Use |
| :--- | :--- | :--- |
| `PYFITTER_DIR` | `pyfitter` (the submodule) | Use another pyfitter checkout, for development. Results are then not tied to the pinned version |
| `MU2E_ANA_VERSION` | `2.8.0` | Another `pyenv ana` version |

To move to a new pyfitter version, check out the tag in the submodule and commit the change:

```bash
git -C pyfitter fetch --tags && git -C pyfitter checkout <tag>
git add pyfitter && git commit -m "Use pyfitter <tag>"
```

## 2. Run the chain

```bash
python run_chain.py --tag run1a --version run1a --group mix --loc tape --jobs 16
```

| Option | Meaning |
| :--- | :--- |
| `--tag` | Name of the output directory |
| `--version` | Cut set from `VERSION_CUTS` in `config.py`. Use `run1a` for the Run-1A analysis cuts, which include the TSDA veto (`no_tsda`) |
| `--group mix` | Run all the mixed signal and background samples in `datasets.yaml` |
| `--loc` | Where the input files are: `tape`, `disk` or `local` |
| `--jobs` | Number of parallel workers |

The datasets, acceptances, expected yields and fit type (`1d` or `2d`) come from
`datasets.yaml`. Don't edit it by hand: regenerate it with
`python make_registry.py`.

Always pass `--version run1a`, including when you rerun part of the chain, so
every stage uses the same cuts.

Add `--dry-run` to print the commands without running them.

## 3. Find the results

Everything is written under `run1a/` (the `--tag` name):

| Path | Contents |
| :--- | :--- |
| `run1a/cutflows/` | Cut flow per dataset |
| `run1a/manifest.yaml` | Acceptances, cut efficiencies and expected yields |
| `run1a/fit/plots/` | Fit plots |
| `run1a/fit/output_card.yaml` | Fitted values passed to the limit |
| `run1a/limits/` | The upper limit |
| `run1a/logs/` | Logs for each stage, check here if something fails |

## Rerunning part of the chain

Selection is the slow step. Once it has run, you can redo later stages on their own:

```bash
# Redo the fit and limit, reusing the selection
python run_chain.py --tag run1a --version run1a --skip-selection

# Redo only the limit
python run_chain.py --tag run1a --version run1a --limit-only

# Selection and fit only, no limit
python run_chain.py --tag run1a --version run1a --group mix --skip-limit
```

Other useful options:

| Option | Meaning |
| :--- | :--- |
| `--fit-type 1d` / `2d` | Override the fit type in the registry |
| `--toys N` | Number of toys for the limit |
| `--asimov` | Fast asymptotic limit (Asimov) instead of toys |
| `--control-region` | Take the cosmic momentum shape from the OffSpill sample |
| `--rmc-control` | Constrain the RMC yield with late-time $e^+$ |

## Cut-and-count cross-check

`cutcount/` gives the limit from a simple count in a fixed momentum and time window,
instead of the full shape fit. It uses the same selection and the same statistics:

1. **Selection.** Calls pyfitter's `run_chain.py --skip-fit`, so `process.py` and
   `analyze.py` apply the cuts of `--version`.
2. **Count.** Counts the post-cut tracks of each dataset in the window. It scales them
   to the run as `run_chain.py` does: `expected × acceptance × N_window / N_generated`
   for MC, and `N_window × (on-spill time / cosmic livetime)` for cosmics.
3. **Limit.** `counting_stats.py` is a one-bin version of pyfitter's
   `simple_limit_combine_card.py`. It uses the same conventions: POI μ = N_CE / card
   CE yield, Combine lnN nuisances, and a one-sided q_μ. It reports Asimov expected
   limits and bands, a toy cross-check, or an empirical toy CLs scan
   (`--method empirical`).

By default the window is the `signal_region` of the cut set: 103.34–104.74 MeV/c and
640–1650 ns for `run1a`. Use `--p-window` and `--t-window` to change it.

```bash
source setup.sh   # sets PYFITTER_DIR

# Full chain: selection, count, limit
python sensitivity/cutcount/run_cutcount.py --tag run1a_count --version run1a --group mix --loc tape --jobs 16

# Reuse an existing run_chain tag (no selection), try another window
python sensitivity/cutcount/run_cutcount.py --tag $PYFITTER_DIR/run1a --version run1a --skip-selection \
    --p-window 103.6 104.9 --t-window 700 1650

# Limit only, from an existing card
python sensitivity/cutcount/counting_stats.py --card run1a_count/cutcount/count_card.yaml --method empirical
```

The window has to lie inside the preselection: 100–110 MeV/c and the
`within_t0_475` window (475–1650 ns for `run1a`). The script warns if it does not.

Everything is written to `<tag>/cutcount/`:

| File | Contents |
| :--- | :--- |
| `count_table.csv` | Counts, efficiencies and yields in the window, per dataset |
| `count_card.yaml` | One-bin card: yields, lnN systematics, R_μe normalisation |
| `count_limits_<method>.yaml` | Limits in events and R_μe |

Other options:

| Option | Meaning |
| :--- | :--- |
| `--mc-stat` | Add a 1/√N lnN for each dataset's MC statistics in the window |
| `--n-obs N` | Also compute an observed limit for N events in the window |
| `--freeze-nuisances` | No systematics |
| `--skip-limit` | Stop after writing the card |

## Summary web page

`make_page.py` writes a single page with no external dependencies. It shows the momentum
and time distributions per process, scaled to the run and with the counting window shaded.
It also shows the expected R_μe limits from the cut and count and from the shape fit.

The page lives in the `sensitivity/display` submodule
([senisitivity_display](https://github.com/sophiemiddleton/senisitivity_display)). After
cloning this repository, fetch the submodule once:

```bash
git submodule update --init sensitivity/display
```

Run `make_page.py` after `run_chain.py` and `cutcount/run_cutcount.py` on the same tag:

```bash
python sensitivity/make_page.py --tag $PYFITTER_DIR/run1a     # writes sensitivity/display/index.html
python sensitivity/make_page.py --tag $PYFITTER_DIR/run1a --out run1a.html --cutcount-dir run1a_count/cutcount
```

Two optional panels compare projections. Each appears only when you pass its input:

| Option | Panel |
| :--- | :--- |
| `--runtime-summary <pyfitter>/studies/runtime/summary.csv` | Limit vs run time, from `run_studies.py --study runtime` |
| `--beam 'LABEL=TAG' ...` | Shape-fit limit of each tag, e.g. 1.5 kW, 1BB and 2BB. A tag without a limit shows as "Not run" |

Both limit tiles show the median expected limit of the toys and its 68% band. Pass `--limits asimov`
for the Asimov median instead.

`--shape-rmue 1.09e-14` shows that shape-fit limit in place of `<tag>/limits`, for a tag whose
shape fit was run elsewhere. The tile then shows only the median, with no band.

```bash
python sensitivity/make_page.py --tag $PYFITTER_DIR/run1a \
    --runtime-summary $PYFITTER_DIR/studies/runtime/summary.csv \
    --beam "1.5 kW=$PYFITTER_DIR/studies/beampower/1BB1500W" "1BB=$PYFITTER_DIR/run1a" "2BB=$PYFITTER_DIR/run1a_2BB"
```

Both panels are also saved as PNG and PDF next to the page (`limit_vs_runtime.{png,pdf}`,
`limit_by_beam.{png,pdf}`) and linked from it, for talks and notes.

The 1.5 kW tag comes from `python run_studies.py --study beampower --points 1BB1500W`.

To publish a new version, commit and push the page in the submodule, then commit the
updated submodule pointer here:

```bash
git -C sensitivity/display commit -am "Update for <tag>" && git -C sensitivity/display push
git add sensitivity/display && git commit -m "Update sensitivity display"
```

The display repository is public. The page shows file paths relative to the tag, not
your area, but it does show the limits and yields.

## Contact

Sophie Middleton (sophie@fnal.gov)
