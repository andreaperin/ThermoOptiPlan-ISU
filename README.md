# ThermoOptiPlan

Uncertainty quantification of an OpenGeoSys (OGS) 3-layer thermal doublet model.
Output: `crossing_year`, the year the mixed production temperature drops 1 K below its initial value.

## Layout

| Folder | Content |
|---|---|
| `src/` | model definition: random inputs, OGS model, post-processing, paths |
| `scripts/` | OGS runs (`run_samples`, `run_nominal`); `surrogates/` fits (`pce`, `pce_wafp`, `ipm`, `lbfm`, shared `settings.jl`); `auxiliary/` dataset and figures |
| `tools/` | `dataset.jl` (inputs + output DataFrame for GP training), `domain_analyzer.jl` |
| `models/` | OGS inputs: `multilayer_irz` (template), `multilayer_irz_nominal` |
| `data/` | generated, not in git: `runs/`, `surrogates/`, `datasets/` |
| `archive/`, `docs/` | 1-layer model, original models, Kim's post-processing |

## Setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python3 -m venv .venv && .venv/bin/pip install -r py_requirements.txt
```

## Run

From the repository root. On the cluster: `sbatch scripts/surrogates/pce.sh` (data goes to `$THERMOOPTIPLAN_DATA`, set in the `.sh`).

## Notes

- Layer naming: `sandstone1/2/3` = `Sandstone_2/3/4` = top/middle/bottom.
- Data in `data/surrogates/` (128/256/512 runs) uses the old input distributions (shared κ/porosity window).
- Future work: pumping rates per layer are fixed in the `.prj`; make them κ-dependent (fixed total, split by κ·h).
