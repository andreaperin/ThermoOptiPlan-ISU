# Models

| Folder | Content |
|---|---|
| `multilayer_irz/` | 3-layer OGS model, template: the 15 random inputs are `{{{ :name }}}` placeholders in `MULTI_BW_line_IRZ.prj` |
| `multilayer_irz_nominal/` | same model with the nominal values (μ below) filled in |

## Layers

| Layer | OGS medium | Mesh group | z [m] | Pumping share |
|---|---|---|---|---|
| top | 1 | `Sand_2` | −1240.2 … −1265.4 | 0.245 |
| middle | 2 | `Sand_3` | −1309.1 … −1338.0 | 0.430 |
| bottom | 3 | `Sand_4` | −1364.6 … −1374.8 | 0.325 |

Solid properties are numbered `sandstone1/2/3`, OGS parameters `_2/_3/_4`; both mean top/middle/bottom.

## Random inputs

All are truncated normal: TN(μ, σ) on [lower, upper]. The name is both the placeholder in the `.prj` and the column name in the datasets. Defined in `src/model_3layers.jl`.

| Name | Property | Layer | μ | σ | Bounds |
|---|---|---|---|---|---|
| `thermal_conductivity_sandstone1` | solid thermal conductivity [W/(m·K)] | top | 2.38 | 0.25 | [1.785, 2.975] |
| `thermal_conductivity_sandstone2` | | middle | 2.38 | 0.25 | [1.785, 2.975] |
| `thermal_conductivity_sandstone3` | | bottom | 2.38 | 0.25 | [1.785, 2.975] |
| `specific_heat_capacity_sandstone1` | solid specific heat capacity [J/(kg·K)] | top | 820 | 80 | [615, 1025] |
| `specific_heat_capacity_sandstone2` | | middle | 820 | 80 | [615, 1025] |
| `specific_heat_capacity_sandstone3` | | bottom | 820 | 80 | [615, 1025] |
| `density_sandstone1` | solid density [kg/m³] | top | 2690 | 100 | [2400, 2800] |
| `density_sandstone2` | | middle | 2690 | 100 | [2400, 2800] |
| `density_sandstone3` | | bottom | 2690 | 100 | [2400, 2800] |
| `sandstone_porosity_parameter_2` | porosity [–] | top | 0.19 | 0.019 | [0.152, 0.228] |
| `sandstone_porosity_parameter_3` | | middle | 0.20 | 0.020 | [0.160, 0.240] |
| `sandstone_porosity_parameter_4` | | bottom | 0.22 | 0.022 | [0.176, 0.264] |
| `kappa_Sandstone_2` | permeability [m²] | top | 1.56e-13 | 1.56e-14 | [1.248e-13, 1.872e-13] |
| `kappa_Sandstone_3` | | middle | 2.39e-13 | 2.39e-14 | [1.912e-13, 2.868e-13] |
| `kappa_Sandstone_4` | | bottom | 4.97e-13 | 4.97e-14 | [3.976e-13, 5.964e-13] |

Porosity and permeability: σ = 10% of μ, bounds ±20%. Thermal conductivity and heat capacity: bounds ±25%.

## Previous distributions (v1)

The runs in `data/surrogates/` (128/256/512) were sampled before 2026-10 with:

| Inputs | v1 distribution |
|---|---|
| thermal conductivity, all layers | TN(2.38, 0.25) on [1.5675, 2.6125] |
| porosity top / middle / bottom | TN(0.19 / 0.20 / 0.22, 0.05) on [0.18, 0.23] |
| permeability top / middle / bottom | TN(1.56 / 2.39 / 4.97 e-13, 2.868e-13) on [1.912e-13, 2.868e-13] |
| heat capacity, density | as above |
