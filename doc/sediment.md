# The `&sediment` Block — Sediment Zones and Bed Heat Exchange

*GLM Science Manual — configuration reference*

---

## 1. Overview

The lake bed is an active boundary of a lake. It exchanges heat with the overlying
water, shades and reflects light, exerts drag, and — when a water-quality module is
coupled — drives benthic biogeochemistry. The **`&sediment`** namelist block configures
this boundary.

It plays **two distinct roles**:

1. **Bed geometry** — whether the bed is treated as a single uniform surface or
   subdivided into a small number of depth-banded **sediment zones**, and the optical /
   roughness properties of those zones.
2. **Sediment temperature model** — how the temperature of the bed is determined and how
   heat is exchanged with the water column. Two approaches are available: a lightweight
   **static** annual cycle, and a physically based **dynamic** 1-D soil-heat solver.

These two roles are independent: you choose a bed geometry with `benthic_mode`, and a
heat-exchange approach with `sed_heat_model`.

> **At a glance**
> | I want… | Set |
> |---|---|
> | No sediment heating | omit `&sediment` (or `sed_heat_model = 0`) |
> | Simple seasonal bed temperature | `sed_heat_model = 1` |
> | One bed for the whole lake | `benthic_mode = 1` |
> | Depth-banded bed zones | `benthic_mode = 2` + `n_zones`, `zone_heights` |
> | A process-based soil column | `sed_heat_model = 2` (requires `benthic_mode = 2` + a WQ module) |

If the `&sediment` block is absent, sediment heating is switched off entirely.

---

## 2. Quick start

**Uniform bed, simple seasonal heating** — the most common starting point:

```fortran
&sediment
   benthic_mode       = 1
   sed_heat_model     = 1
   sed_temp_mean      = 12.0      ! annual-mean bed temperature (degC)
   sed_temp_amplitude = 6.0       ! seasonal swing (degC)
   sed_temp_peak_doy  = 60        ! day-of-year of warmest bed
   sed_heat_Ksoil     = 5.0       ! bed-water coupling conductivity
   sed_temp_depth     = 0.1       ! conductive depth scale (m)
/
```

**Two depth zones, zoned seasonal heating:**

```fortran
&sediment
   benthic_mode       = 2
   n_zones            = 2
   zone_heights       = 5.0, 9.5          ! top of each zone (m above bed datum)
   sed_temp_mean      = 11.4, 14.9        ! per-zone
   sed_temp_amplitude = 3.1, 7.1
   sed_temp_peak_doy  = 278, 277
   sed_heat_model     = 1
/
```

**Dynamic soil column** (requires `benthic_mode = 2` and an active WQ module):

```fortran
&sediment
   benthic_mode   = 2
   n_zones        = 2
   zone_heights   = 5.0, 9.5
   sed_heat_model = 2
   n_sed_layers   = 8                      ! soil nodes (>= 3)
   sed_temp_depth = 1.5                    ! soil-column depth (m)
   sed_temp_mean  = 11.4, 14.9             ! deep-soil temperature per zone
   sed_vwc        = 0.4                    ! volumetric water content
/
```

---

## 3. Bed geometry — `benthic_mode`

`benthic_mode` selects how the bed is represented.

| `benthic_mode` | Bed representation |
|:---:|---|
| **1** | **Uniform bed.** A single set of sediment parameters is applied to the whole lake. Every water layer exchanges heat with the bed through the sediment area it overlies. |
| **2** | **Sediment zones.** The bed is divided into `n_zones` depth bands. Each water layer is matched to the zone it sits above, and is heated using *that zone's* parameters. **Required** for the dynamic soil model (`sed_heat_model = 2`). |
| **3** | **AED zonation variant.** As for mode 2 for water-quality zone-averaging, with an additional biogeochemical variant on the AED side. Has no separate effect on GLM's core sediment-heat calculation. |

Notes:

- Modes **2 and 3 require sediment zones** to be defined (`n_zones > 0` and
  `zone_heights`). If zones are missing, GLM reverts to a uniform bed (mode 1) and warns.
- Parameters such as `sed_temp_mean` are read as **lists, one value per zone**, in
  mode 2/3. In mode 1 only the first value is used.

### Defining zones

Sediment zones are depth bands of the bed, defined bottom-up by their **top heights**:

| Key | Type | Units | Description |
|---|---|---|---|
| `n_zones` | integer | – | Number of sediment zones. |
| `zone_heights` | list | m | Height of the **top** of each zone above the bed datum, in ascending order. Zone *z* spans from `zone_heights[z-1]` up to `zone_heights[z]`. |
| `sed_reflectivity` | list | – (0–1) | Short-wave reflectivity of each zone's bed surface. |
| `sed_roughness` | list | m | Hydraulic roughness length of each zone's bed. |

At every step GLM maps each water layer to a zone by comparing the layer height to the
zone top heights, then applies that zone's sediment properties to the layer.

---

## 4. Sediment temperature models — `sed_heat_model`

`sed_heat_model` selects how bed temperature and bed→water heat flux are computed.

| `sed_heat_model` | Approach |
|:---:|---|
| **0** | Off — no sediment heating. |
| **1** | **Static**: bed temperature follows a prescribed annual sinusoid; heat is exchanged with the water through a fixed conductive layer. Cheap, robust, no extra state. |
| **2** | **Dynamic**: a 1-D heat-conduction equation is solved through a layered soil column beneath the bed, driven by the overlying water above and a deep soil temperature below. Process-based; resolves the sediment thermal store and phase lag. |

### 4.1 Static model (`sed_heat_model = 1`)

The bed temperature is prescribed as an **annual cosine cycle**:

```
T_sed(d) = sed_temp_mean
         + sed_temp_amplitude · cos( 2π (d − sed_temp_peak_doy) / 365 )
```

where `d` is the day of year. Heat is exchanged with each overlying water layer as a
conductive flux through a sediment layer of effective conductivity `sed_heat_Ksoil`
(`KSED`) and thickness `sed_temp_depth` (`ZSED`):

```
q = KSED · ( T_sed − T_water ) / ZSED          [W m⁻²]
```

This flux, integrated over the sediment area each layer overlies, warms or cools that
layer. In zoned mode each zone carries its own `sed_temp_mean`, `sed_temp_amplitude` and
`sed_temp_peak_doy`.

| Key | Units | Default | Per-zone | Description |
|---|---|:---:|:---:|---|
| `sed_temp_mean` | °C | – | yes | Annual-mean bed temperature. *(Also the deep-soil temperature for model 2.)* |
| `sed_temp_amplitude` | °C | – | yes | Amplitude of the annual cycle (half peak-to-trough). |
| `sed_temp_peak_doy` | day | – | yes | Day of year of the warmest bed temperature. |
| `sed_heat_Ksoil` | W m⁻¹ K⁻¹ | 5.0 | no | Effective conductivity of the bed–water coupling layer (`KSED`). |
| `sed_temp_depth` | m | 0.1 | no | Conductive depth scale of the coupling layer (`ZSED`). |

**Use the static model when** you want a simple, well-behaved seasonal bed signal, lack
data to parameterise a full soil column, or are not running a water-quality module.

### 4.2 Dynamic model (`sed_heat_model = 2`)

The dynamic model solves the **1-D heat-conduction equation** through a soil column
beneath each sediment zone:

```
C(θ) ∂T/∂t = ∂/∂z [ λ(θ) ∂T/∂z ] ,        0 ≤ z ≤ D
```

where `z` is depth into the bed, `D = sed_temp_depth` is the column depth, `θ` is the
volumetric water content, `λ` is the soil thermal conductivity and `C` is the volumetric
heat capacity. The boundaries are:

- **Top (`z = 0`):** the temperature of the overlying water above that zone (Dirichlet).
- **Bottom (`z = D`):** a fixed deep-soil temperature, `sed_deep_temp` (Dirichlet); if
  unset it defaults to the zone's `sed_temp_mean`.

The heat flux across the bed surface is returned to the water column each step, closing
the bed–water energy exchange.

**Numerics.** The column is discretised as a **cell-centred finite-volume** grid and
advanced with a **fully implicit (backward-Euler)** step, which is unconditionally
stable. Interface conductivities use a harmonic mean of the adjacent cells, and the
resulting tridiagonal system is solved directly. The model time step is the GLM
hydrodynamic time step.

**Soil thermal properties.** Conductivity follows a **Johansen-style** moisture
dependence and heat capacity a volumetric mixing rule:

```
porosity      φ     = 1 − ρ_bulk / ρ_mineral
saturation    Sr    = θ / φ
Kersten no.   Ke    = clip( log₁₀(Sr) + 1 , 0 , 1 )
dry cond.     k_dry = (0.135 ρ_bulk + 64.7) / (ρ_mineral − 0.947 ρ_bulk)   [ρ in kg m⁻³]
sat. cond.    k_sat = k_mineral^(1−φ) · k_water^φ
conductivity  λ(θ)  = k_dry + (k_sat − k_dry) · Ke
heat capacity C(θ)  = (1−φ) c_mineral + θ c_water + (φ−θ) c_air
```

**The soil grid.** You set the number of nodes (`n_sed_layers`) and the column depth
(`sed_temp_depth`); GLM builds a **geometric grid that is fine near the surface and
coarsens with depth** (a 5 mm surface cell, growing to span the full column). This
resolves the steep near-surface thermal gradient where the diurnal/seasonal signal is
strongest. `n_sed_layers` is the **total** node count *N* (including the surface and deep
boundary nodes); the solver advances the *N − 2* interior cells. To prescribe the grid
explicitly instead, provide `sed_layer_depth` (a list of *N* node depths); it overrides
the automatic grid.

**Spin-up.** Before the run, the profile is relaxed toward quasi-equilibrium for
`sed_spinup_days` so the simulation does not start from an arbitrary state.

| Key | Units | Default | Description |
|---|---|:---:|---|
| `n_sed_layers` | – | – | Total soil nodes *N* (≥ 3), including both boundaries. |
| `sed_temp_depth` | m | 0.1 | Soil-column depth *D*. *(Set deeper — e.g. 1–2 m — for the dynamic model.)* |
| `sed_layer_depth` | list, m | auto | *Optional.* Explicit node depths (*N* values). Omit to auto-build the geometric grid. |
| `sed_vwc` | m³ m⁻³ | 0.4 | Volumetric water content; one value, or one per node. |
| `sed_spinup_days` | days | 365 | Spin-up duration before the run. |
| `sed_deep_temp` | °C | `sed_temp_mean` | Deep-boundary soil temperature. |
| `sed_k_mineral` | W m⁻¹ K⁻¹ | 2.5 | Thermal conductivity of the mineral fraction. |
| `sed_k_water` | W m⁻¹ K⁻¹ | 0.57 | Thermal conductivity of water. |
| `sed_k_air` | W m⁻¹ K⁻¹ | 0.025 | Thermal conductivity of air. |
| `sed_c_mineral` | J m⁻³ K⁻¹ | 2.0×10⁶ | Volumetric heat capacity of the mineral fraction. |
| `sed_c_water` | J m⁻³ K⁻¹ | 4.18×10⁶ | Volumetric heat capacity of water. |
| `sed_c_air` | J m⁻³ K⁻¹ | 1.25×10³ | Volumetric heat capacity of air. |
| `sed_bulk_density` | Mg m⁻³ | 1.5 | Dry bulk density of the soil. |
| `sed_mineral_density` | Mg m⁻³ | 2.6 | Density of the mineral grains. |
| `sed_porosity` | – (0–1) | derived | Porosity; if unset, computed as `1 − bulk/mineral density`. |

All thermal properties are **optional** — the defaults reproduce a typical mineral soil
and match the reference implementation (see §8). Most users need only set
`n_sed_layers`, `sed_temp_depth`, `sed_vwc` and the deep temperature.

**Requirements.** The dynamic model:

- requires **`benthic_mode = 2`** (it operates per sediment zone), and
- requires an **active water-quality module** (`&wq_setup wq_lib = 'aed'` or `'api'`),
  which provides the soil solver. Running `sed_heat_model = 2` without a WQ module is a
  configuration error.

**Use the dynamic model when** you need a realistic sediment heat store — for example to
capture the thermal inertia and seasonal phase lag of the bed, near-bed temperatures, or
to couple bed temperature to benthic biogeochemistry.

---

## 5. Outputs

- **Static model (`= 1`):** the effect appears in the simulated water temperatures
  (`output.nc`); the model carries no additional state.
- **Dynamic model (`= 2`):** the per-zone soil column is written to the output NetCDF as
  time series, for inspection and validation:
  - `sed_temp(time, zone, layer)` — soil temperature profile (°C)
  - `sed_heatflux(time, zone)` — bed-surface heat flux (W m⁻²)
  - `zone_ztemp(time, zone)` — overlying-water temperature used as the top boundary (°C)
  - `sed_layer_depth(layer)` — the soil node depths (m)

---

## 6. Choosing and setting up

- **Start simple.** If you only need a seasonal bed signal, use `sed_heat_model = 1`. It
  needs only three temperatures (`sed_temp_mean`, `sed_temp_amplitude`,
  `sed_temp_peak_doy`) and is hard to misconfigure.
- **Move to the dynamic model** when the bed's thermal *memory* matters — strong seasonal
  phase lag, near-bed temperature dynamics, or coupling to sediment biogeochemistry.
- **For the dynamic model**, set a realistic column depth (`sed_temp_depth`, ~1–2 m so the
  deep boundary is below the seasonal penetration depth), a node count that resolves the
  surface (8–20 nodes is usually ample given the surface-refined grid), and a deep
  temperature (`sed_deep_temp`, or rely on `sed_temp_mean`). Increase `sed_spinup_days` if
  the column is deep.
- **Zones** let bed properties vary with depth (e.g. littoral vs profundal). Define them
  once with `n_zones` / `zone_heights`; all per-zone lists must then have `n_zones` values.

---

## 7. Full parameter summary

| Key | Role | Units | Default | Notes |
|---|---|---|:---:|---|
| `benthic_mode` | geometry | – | 1 | 1 uniform · 2 zones · 3 AED variant |
| `n_zones` | geometry | – | 0 | number of sediment zones |
| `zone_heights` | geometry | m | – | per-zone top heights (ascending) |
| `sed_reflectivity` | geometry | – | – | per-zone bed reflectivity |
| `sed_roughness` | geometry | m | – | per-zone bed roughness |
| `sed_heat_model` | selector | – | 0 | 0 off · 1 static · 2 dynamic |
| `sed_temp_mean` | static + dynamic | °C | – | per-zone mean; deep-soil temp for model 2 |
| `sed_temp_amplitude` | static | °C | – | per-zone annual amplitude |
| `sed_temp_peak_doy` | static | day | – | per-zone peak day-of-year |
| `sed_heat_Ksoil` | static | W m⁻¹ K⁻¹ | 5.0 | coupling conductivity `KSED` |
| `sed_temp_depth` | static + dynamic | m | 0.1 | `ZSED` (model 1) / column depth (model 2) |
| `n_sed_layers` | dynamic | – | – | total soil nodes *N* (≥ 3) |
| `sed_layer_depth` | dynamic | m | auto | optional explicit node depths |
| `sed_vwc` | dynamic | m³ m⁻³ | 0.4 | volumetric water content |
| `sed_spinup_days` | dynamic | days | 365 | spin-up duration |
| `sed_deep_temp` | dynamic | °C | `sed_temp_mean` | deep boundary temperature |
| `sed_k_mineral` / `sed_k_water` / `sed_k_air` | dynamic | W m⁻¹ K⁻¹ | 2.5 / 0.57 / 0.025 | phase conductivities |
| `sed_c_mineral` / `sed_c_water` / `sed_c_air` | dynamic | J m⁻³ K⁻¹ | 2.0e6 / 4.18e6 / 1.25e3 | phase heat capacities |
| `sed_bulk_density` / `sed_mineral_density` | dynamic | Mg m⁻³ | 1.5 / 2.6 | soil densities |
| `sed_porosity` | dynamic | – | derived | `1 − bulk/mineral` if unset |

*Keys marked “per-zone” are read as one value per zone under `benthic_mode = 2/3`, and as
a single value under `benthic_mode = 1`.*

---

## 8. Theoretical basis — the soil heat-conduction model

The dynamic model (§4.2) is GLM's implementation of the **subsurface heat-conduction**
component of the AED Wetland Soil Model. This section sets out its theoretical basis.
In GLM the lake bed always lies beneath the water column, so the surface boundary is
always the **overlying water temperature** (the *submerged* case). The wetland-soil model's
surface energy balance for an *exposed*, drying bed, and its dynamic soil-moisture solver,
are summarised in §8.7 as context and planned extensions; they are not active in GLM
today (moisture content is prescribed via `sed_vwc`).

### 8.1 Governing equation

Vertical heat conduction in the bed follows Fourier's law. The one-dimensional heat
equation for the soil column is

$$
C(\theta)\,\frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\!\left[\lambda(\theta)\,\frac{\partial T}{\partial z}\right]
$$

where $T$ is temperature (°C), $z$ is depth below the bed surface (m), $C(\theta)$ is the
volumetric heat capacity (J m⁻³ K⁻¹), $\lambda(\theta)$ is the thermal conductivity
(W m⁻¹ K⁻¹) and $\theta$ is the volumetric water content (m³ m⁻³). Because both $C$ and
$\lambda$ depend on moisture, the thermal state is tied to the soil's water content.

### 8.2 Spatial discretisation (cell-centred finite volume)

The column is divided into $n$ cells with geometric spacing — fine near the surface,
coarsening with depth (§4.2). The layer boundaries are the depth nodes $d[0..n+1]$, with
$d[0]=0$ at the surface and $d[n+1]$ at the column base. Cell thicknesses and centre
depths are

$$
\Delta z_i = d[i+1] - d[i], \qquad z_{c,i} = \frac{d[i] + d[i+1]}{2}, \qquad i = 0,\ldots,n-1 .
$$

Inter-cell distances and interface conductivities use the **harmonic mean** for internal
interfaces (the correct series-resistance average for conduction across cells of differing
conductivity):

$$
\Delta z_{c,i} = z_{c,i} - z_{c,i-1}, \qquad
\lambda_{i+\frac12} = \frac{2\,\lambda_i\,\lambda_{i-1}}{\lambda_i + \lambda_{i-1}} .
$$

The top interface takes the top-cell conductivity, $\lambda_{-\frac12}=\lambda_0$, with the
distance from the surface to the first cell centre $\Delta z_{c,0}=z_{c,0}$.

### 8.3 Time integration (implicit backward Euler)

An implicit backward-Euler step is used, giving a tridiagonal system that is
**unconditionally stable** for any time step or cell thickness:

$$
C_i\,\Delta z_i\,\frac{T_i^{\,n+1} - T_i^{\,n}}{\Delta t}
= \frac{\lambda_{i-\frac12}}{\Delta z_{c,i}}\left(T_{i-1}^{\,n+1} - T_i^{\,n+1}\right)
+ \frac{\lambda_{i+\frac12}}{\Delta z_{c,i+1}}\left(T_{i+1}^{\,n+1} - T_i^{\,n+1}\right) .
$$

The system is solved directly by the Thomas (tridiagonal) algorithm each step. The time
step $\Delta t$ is the GLM hydrodynamic time step.

### 8.4 Boundary conditions (as used in GLM)

- **Top — Dirichlet.** The surface temperature is prescribed as the temperature of the
  water overlying that sediment zone, $T_{-1}^{\,n+1} = T_s = T_{\text{water}}$. Because the
  hydrodynamic model has already solved the energy balance of the water column, the bed
  simply conducts heat to or from the saturated soil beneath it. *(In the fuller wetland-soil
  model, when the bed is exposed to air, $T_s$ instead comes from a surface energy balance —
  see §8.7.)*
- **Bottom — Dirichlet.** A fixed deep-soil temperature $T_{\text{deep}}$ is held at the
  column base (`sed_deep_temp`, defaulting to the zone's `sed_temp_mean`), representing the
  quasi-constant temperature below the seasonal penetration depth.

The bed-surface heat flux returned to the water column is computed diagnostically from the
updated profile,

$$
G = \frac{\lambda_0}{\Delta z_{c,0}}\left(T_s - T_0^{\,n+1}\right) \qquad \text{[W m⁻², positive into the soil]} .
$$

### 8.5 Thermal properties

**Thermal conductivity — Johansen model.** The effective conductivity interpolates between
dry and saturated values with the Kersten number $K_e$:

$$
\lambda_{\text{eff}} = \lambda_{\text{dry}} + \left(\lambda_{\text{sat}} - \lambda_{\text{dry}}\right) K_e ,
\qquad
K_e = \log_{10}(S_r) + 1, \quad K_e \in [0,1], \quad S_r = \theta/\phi ,
$$

where $S_r$ is the degree of saturation and $\phi$ the porosity. The dry and saturated
end-members (Johansen, 1975; geometric-mean saturated value) are

$$
\lambda_{\text{dry}} = \frac{0.135\,\rho_b + 64.7}{\rho_s - 0.947\,\rho_b},
\qquad
\lambda_{\text{sat}} = \lambda_m^{\,(1-\phi)}\,\lambda_w^{\,\phi},
$$

with $\rho_b$ the bulk density and $\rho_s$ the mineral density (both kg m⁻³).

**Volumetric heat capacity — volume-weighted mixture:**

$$
C = (1-\phi)\,C_m + \theta\,C_w + (\phi-\theta)\,C_a ,
\qquad \phi = 1 - \frac{\rho_b}{\rho_s} .
$$

| Property | Symbol | Default | Units | Key |
|---|---|---|---|---|
| Mineral conductivity | $\lambda_m$ | 2.5 | W m⁻¹ K⁻¹ | `sed_k_mineral` |
| Water conductivity | $\lambda_w$ | 0.57 | W m⁻¹ K⁻¹ | `sed_k_water` |
| Mineral heat capacity | $C_m$ | 2.0×10⁶ | J m⁻³ K⁻¹ | `sed_c_mineral` |
| Water heat capacity | $C_w$ | 4.18×10⁶ | J m⁻³ K⁻¹ | `sed_c_water` |
| Air heat capacity | $C_a$ | 1.25×10³ | J m⁻³ K⁻¹ | `sed_c_air` |
| Mineral density | $\rho_s$ | 2.6 | Mg m⁻³ | `sed_mineral_density` |
| Bulk density | $\rho_b$ | 1.5 | Mg m⁻³ | `sed_bulk_density` |

**Physical implication.** Both $\lambda$ and $C$ rise strongly with moisture: a saturated
mineral soil conducts heat roughly five times faster and stores about twice the heat of the
same soil when dry. Wetter beds therefore transmit the surface (water) signal deeper and
respond more sluggishly; drier beds behave as surface-confined insulators. This moisture
control is why `sed_vwc` matters even though it is currently prescribed.

### 8.6 Initialisation (spin-up)

Before the run, the profile is relaxed toward quasi-equilibrium by integrating the heat
equation with constant surface and deep temperatures for `sed_spinup_days`, starting from a
linear interpolation between the surface and deep temperatures. This avoids transients from
an arbitrary initial profile.

**Numerical note.** Backward Euler is unconditionally stable but first-order in time; at an
hourly step the solution is accurate to ~0.1 °C over a diurnal cycle. A Dirichlet
(temperature) surface condition is used rather than a prescribed-flux (Neumann) condition,
because the very thin surface cell (~5 mm) makes a flux boundary numerically oscillatory.

### 8.7 The fuller wetland-soil model and planned GLM extensions

GLM currently implements the heat-conduction core above with **prescribed** moisture and an
always-submerged surface. The complete AED Wetland Soil Model adds two coupled parts that
are candidates for future GLM coupling:

- **Surface energy balance (SEB).** When the bed is *exposed* to air, the surface
  temperature is not simply the water temperature but the solution of an energy balance,
  $SW_{net} + LW_{down} - LW_{out} - H - LE = G$, combining net short- and long-wave
  radiation, sensible ($H$) and latent ($LE$) turbulent fluxes and the ground flux $G$,
  solved for $T_s$ by Newton iteration. Evaporative cooling is moisture-limited by a factor
  $\beta = \min(\theta_{top}/\theta_{fc},\,1)$. This is what lets an exposed flat heat to
  50–60 °C. In GLM this pathway is inactive because the bed is always submerged.
- **Soil moisture.** A Campbell (1985) unsaturated-flow solver, water-table tracking and an
  optional macropore domain make $\theta$ prognostic (rather than the fixed `sed_vwc`),
  which then feeds back on $\lambda$ and $C$. Coupling this — so the top bed layer can dry
  above the water line and zones can dry as a lake recedes — is planned work.

Both are documented in full in the AED Wetland Soil Model Technical Manual (Parts A and B).

---

## 9. Background and references

GLM's dynamic sediment-temperature model is the heat-conduction component of the
**AED Wetland Soil Model**, a 1-D finite-volume soil heat-and-moisture solver. The
conductivity and heat-capacity formulations, the numerical scheme and the default
properties follow that reference implementation, and GLM reproduces it to numerical
precision on a matched configuration.

- AED Wetland Soil Model — <https://github.com/AquaticEcoDynamics/intertidal-soil-py>
  (see `docs/` — *AED Wetland Soil Model Technical Manual*, Part A: Soil Moisture,
  Part B: Soil Temperature and Surface Energy Balance).
- Johansen, O. (1975) *Thermal conductivity of soils.* PhD thesis, University of Trondheim.
  (dry/saturated effective-conductivity scheme.)
- Campbell, G. S. (1985) *Soil Physics with BASIC: Transport Models for Soil-Plant Systems.*
  Elsevier. (soil thermal and hydraulic background.)
