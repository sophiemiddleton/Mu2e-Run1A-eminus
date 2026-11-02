# DIO: Decay-in-Orbit Background

A muon bound in aluminum can decay in orbit. Nuclear recoil gives the electron spectrum a tail that reaches the CE endpoint, falling roughly as $(E_{\text{end}} - E)^5$. This makes DIO the irreducible background just below the signal peak. Its time distribution matches CE.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there (`main`, `d9be2e8`):

- Single-component fit: `component_fits/dio.py`
- Model definition: `model/physics_components.py` (look for DIO). Momentum, by default: `poly58`, the 2011
  polynomial $a_5\delta^5 + \dots + a_8\delta^8$ in $\delta = m_\mu - p - p^2/2m_\text{Al}$, with the coefficients and
  $N_\text{DIO}$ floating (`custom_models.py`). Time: the muonic-aluminum exponential, 864 ns, fixed in the 2D fit
  (`momentum_pdf_builder.DECAY_RATE_CE_DIO`)
- Alternative momentum model: `dio_tabulated_morph`, the production DIO table with two Gaussian-constrained
  endpoint-shape directions (`dio_shape_model.py`, `documentation/DIOShapeModel.md`). Implemented and validated,
  but not the default
- Normalization: `DIO_mix` (`DIOtail95Mix1BB`) = $N_{\text{stop}\,\mu}$ × 0.391 (DIO per stop) × the fraction of the
  Czarnecki/Szafron spectrum above 95 MeV (`documentation/DatasetRegistry.md`)
- Systematics: `lumi` (10%) and `dioN` (2.5%, spectrum shape and rate) on the input card
- Background: `documentation/LikelihoodDefinition.md`

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

### Studies in pyfitter

Done in pyfitter at `d9be2e8`. The yields are from `au_run1a_check`, the Run-1A selection with the d0 fix and TSDA
veto.

| Study | Result | Code |
| :--- | :--- | :--- |
| Endpoint shape uncertainty | The production table `heeck_finer_binning_2016_szafron.tbl` (10 keV bins, SHA256-pinned), deformed as $S_\text{nom}\,e^{\theta_1 v_1 + \theta_2 v_2}$. The directions $v_1$, $v_2$ are the two leading SVD components of endpoint-theory alternatives (radiative corrections, finite nuclear size), and are zero more than 1.96 MeV below the endpoint. Adopted constraints: $\sigma(\theta_1)$ = 0.25, $\sigma(\theta_2)$ = 0.01. Closure and signal-injection toys recover the injected CE without bias (worst case about 0.03 events), and are insensitive to doubling the widths or to ±0.9 keV endpoint shifts. The directions must not float unconstrained: free, they absorb injected signal | `dio_shape_model.py`, `common/dio_shape_morph_basis.csv`, `unit_test/test_dio_shape_model.py`. The construction is in `github.com/zwl0331/DIO_Shape_Systematic` (Finding 011) |
| Shape systematics in the profile | `DIO_Shape_Theta1` and `DIO_Shape_Theta2` profile each direction at ±1σ. They replace the planned `DIO_PDF_Variant` template swap, and only apply with the tabulated-morph model | `uncertainties/model/sysunc_components.py`, `uncertainties/run_profile_systematics.py` |
| DIO in the windows | 42.6 events expected in the fit window (100–110 MeV/c, 475–1650 ns). 0.032 events in the cut-and-count window (103.34–104.74 MeV/c, 540–1650 ns), from 10 MC tracks | `au_run1a_check/input_card.yaml`, `cutcount/count_table.csv` |

### Plots in pyfitter

`paper_plots.py` has no DIO-specific plot. DIO is in the stacked MC of `data_mc` ([sensitivity](../sensitivity/README.md#plots-in-pyfitter)).

## Open items

- Make `dio_tabulated_morph` the default: a one-line change at the end of `model/physics_components.py`
  (`documentation/DIOShapeModel.md`). The default `poly58` floats its coefficients, so the endpoint shape
  is set by the fit, not by theory
- Correlation with the CE momentum response: the morph model enters in reconstructed momentum. Folding the
  deformed truth spectrum through the response (`theo_exp`, or a measured response matrix) is the planned follow-up
- DIO yield: `DIO_Theory` (2.5% on $N_\text{DIO}$) is on hold in `sysunc_components.py`. `dioN` on the card
  is 1.025, with a comment in `uncertainties_config.py` saying 10%
- MC statistics in the signal window: 10 DIO MC tracks, about 0.003 events each
