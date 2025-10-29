# PPRO Operator Configuration Guide

This document describes how to configure the PPRO (Parameterized Polarimetric Radar Operator) in JEDI-UFO for different microphysics schemes and polarimetric operators.

## Overview

The PPRO operator supports:
- **Polarimetric Operators**: 
  - Zhang21 (polynomial functions with coefficients from T-matrix simulations)
  - TCWA2 (analytical gamma distribution formulation)
- **Microphysics Schemes**: Thompson, WSM6, NSSL, TCWA2
- **Radar Bands**: S-band (iband=1), C-band (iband=2)
- **Output Variables**: Zhh, ZDR, KDP, ρhv

## Basic Configuration Structure

```yaml
obs operator:
  name: PPRO
  polarimetric operator: <Zhang21|TCWA2>    # Choose operator type
  microphysics option: <Thompson|WSM6|NSSL|TCWA2>
  VertCoord: height
  # Variable mappings (depends on microphysics scheme)
```

---

## Configuration Examples by Microphysics Scheme

### 1. WSM6 (Single-Moment) with Zhang21 Operator

**Required GeoVaLs: 7 variables**
- Basic: air_density, air_temperature, air_pressure, specific_humidity
- Hydrometeors: qr, qs, qg

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21        # Polynomial-based from T-matrix (default)
  microphysics option: WSM6
  VertCoord: height
  var_rain_mixing_ratio: qr              # rain water mixing ratio (kg/kg)
  var_snow_mixing_ratio: qs              # snow mixing ratio (kg/kg)
  var_graupel_mixing_ratio: qg           # graupel mixing ratio (kg/kg)
```

**Notes:**
- Simplest configuration with 7 GeoVaLs
- Uses fixed intercept parameters N₀
- Zhang21 operator uses polynomial functions of Dm and fx fitted from T-matrix simulations
- Temperature-based melting layer QC applied internally

---

### 2. Thompson (Double-Moment for Rain) with Zhang21 Operator

**Required GeoVaLs: 8 variables**
- WSM6 variables (7) + nr

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21         # Polynomial-based from T-matrix (default)
  microphysics option: Thompson
  VertCoord: height
  var_rain_mixing_ratio: qr
  var_snow_mixing_ratio: qs
  var_graupel_mixing_ratio: qg
  var_rain_number_concentration: nr      # rain number concentration (#/kg)
```

**Notes:**
- Adds rain number concentration for improved rain representation
- 8 total GeoVaLs required
- Zhang21 operator handles both single-moment (WSM6) and double-moment (Thompson) schemes

---

### 3. NSSL (Double-Moment with Hail) with Zhang21 Operator

**Required GeoVaLs: 14 variables**
- Basic (4) + mass mixing ratios (4) + number concentrations (5) + volume mixing ratios (2)

```yaml
obs operator:
  name: PPRO
  polarimetric operator: Zhang21         # Polynomial-based from T-matrix (default)
  microphysics option: NSSL
  VertCoord: height
  # Mass mixing ratios
  var_rain_mixing_ratio: qr
  var_snow_mixing_ratio: qs
  var_graupel_mixing_ratio: qg
  var_hail_mixing_ratio: qh              # hail mixing ratio (kg/kg)
  # Number concentrations
  var_rain_number_concentration: nr       # #/kg
  var_snow_number_concentration: ns       # #/kg
  var_graupel_number_concentration: ng    # #/kg
  var_hail_number_concentration: nh       # #/kg
  # Volume mixing ratios
  var_graupel_vol_mixing_ratio: vg       # m³/kg
  var_hail_vol_mixing_ratio: vh          # m³/kg
```

**Notes:**
- Most comprehensive scheme with 14 GeoVaLs
- Includes hail category for severe weather
- Requires volume mixing ratios for graupel and hail

---

### 4. TCWA2 (Taiwan CWA Analytical Scheme)

**Required GeoVaLs: 15 variables**
- Basic (4) + mass mixing ratios (5) + number concentrations (4) + melted fractions (2)

```yaml
obs operator:
  name: PPRO
  polarimetric operator: TCWA2           # MUST use TCWA2 operator!
  microphysics option: TCWA2             # MUST use TCWA2 microphysics!
  VertCoord: height
  # Mass mixing ratios
  var_rain_mixing_ratio: qr              # kg/kg
  var_snow_mixing_ratio: qs              # kg/kg
  var_graupel_mixing_ratio: qg           # kg/kg
  var_cloud_mixing_ratio: qc             # cloud water (kg/kg)
  var_ice_mixing_ratio: qi               # ice crystal (kg/kg)
  # Number concentrations
  var_rain_number_concentration: nr       # #/kg
  var_snow_number_concentration: ns       # #/kg
  var_graupel_number_concentration: ng    # #/kg
  var_ice_number_concentration: ni        # #/kg
  # Melted fractions (unique to TCWA2)
  var_snow_melted_fraction: smlf         # snow melting fraction [0-1]
  var_graupel_melted_fraction: gmlf      # graupel melting fraction [0-1]
```

**Notes:**
- **CRITICAL**: Must use `polarimetric operator: TCWA2` AND `microphysics option: TCWA2`
- Uses analytical gamma distribution formulation (no polynomial coefficients needed)
- Requires melted fraction variables for melting layer treatment
- Total 15 GeoVaLs required

**TCWA2 Features:**
- No coefficient files needed (analytical formulation)
- Requires melting layer information (smlf, gmlf)
- Designed for Taiwan CWA cloud microphysics scheme
- Complete ice crystal treatment (qi, ni)

---

## GeoVaLs Summary Table

| Scheme    | # GeoVaLs | Basic | qr,qs,qg | qh | qc,qi | nr,ns,ng | nh | ni | vg,vh | smlf,gmlf |
|-----------|-----------|-------|----------|----|----|----------|----|----|-------|-----------|
| **WSM6**     | 7         | 4     | ✓        |    |    |          |    |    |       |           |
| **Thompson** | 8         | 4     | ✓        |    |    | nr only  |    |    |       |           |
| **NSSL**     | 14        | 4     | ✓        | ✓  |    | ✓        | ✓  |    | ✓     |           |
| **TCWA2**    | 15        | 4     | ✓        |    | ✓  | ✓        |    | ✓  |       | ✓         |

**Legend:**
- Basic (4): air_density, air_temperature, air_pressure, specific_humidity
- qr/qs/qg/qh/qc/qi: mass mixing ratios (kg/kg)
- nr/ns/ng/nh/ni: number concentrations (#/kg)
- vg/vh: volume mixing ratios (m³/kg)
- smlf/gmlf: melted fractions [0-1]

---

## Operator Type Selection

### Zhang21 (Polynomial-Based from T-Matrix Simulations)
- **Method**: Polynomial functions of Dm (mean diameter) and fx (water fraction) fitted from T-matrix simulations
- **Use with**: WSM6, Thompson, NSSL
- **Pros**: Validated, includes melting layer QC, efficient polynomial evaluation
- **Cons**: Requires coefficient files for polynomials
- **Config**: `polarimetric operator: Zhang21` (or omit, it's default)

### TCWA2 (Analytical Gamma Distribution)
- **Method**: Analytical formulation based on gamma distribution PSD parameters
- **Use with**: TCWA2 microphysics **ONLY**
- **Pros**: No polynomial coefficients needed, direct gamma PSD calculation
- **Cons**: Requires melted fraction fields explicitly
- **Config**: `polarimetric operator: TCWA2` (REQUIRED!)

---

## Common Issues and Solutions

### 1. Missing GeoVaLs
**Error**: "Variable xxx not found in GeoVaLs"
**Solution**: Ensure your model output includes all required variables for your chosen scheme (see table above)

### 2. Wrong Operator/Scheme Combination
**Error**: "TCWA2 operator requires qi, ni, qc, smlf, gmlf, nr, ns, ng"
**Solution**: 
- If using TCWA2 microphysics → Must set `polarimetric operator: TCWA2`
- If using Thompson/WSM6/NSSL → Use `polarimetric operator: Zhang21` (default)

### 3. Missing Melted Fraction Fields
**Error**: Variable not found for smlf/gmlf
**Solution**: 
- These are required for TCWA2 only
- Other schemes (Thompson/WSM6/NSSL) don't need these fields

---

## Complete Example: TCWA2 with All Settings

```yaml
observations:
  observers:
  - obs space:
      name: Radar
      obsdatain:
        engine:
          type: H5File
          obsfile: radar_obs.nc4
      simulated variables: [equivalentReflectivityFactor, differentialReflectivity, 
                            specificDifferentialPhase, copolarCorrelationCoefficient]
    obs operator:
      name: PPRO
      polarimetric operator: TCWA2
      microphysics option: TCWA2
      VertCoord: height
      # Variable mappings
      var_rain_mixing_ratio: qr
      var_snow_mixing_ratio: qs
      var_graupel_mixing_ratio: qg
      var_cloud_mixing_ratio: qc
      var_ice_mixing_ratio: qi
      var_rain_number_concentration: nr
      var_snow_number_concentration: ns
      var_graupel_number_concentration: ng
      var_ice_number_concentration: ni
      var_snow_melted_fraction: smlf
      var_graupel_melted_fraction: gmlf
    obs error:
      covariance model: diagonal
```

---

## References

1. **Zhang21 Operator**:  
   Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. **TCWA2 Operator**:  
   - Chen, J.-P., and D. Lamb, 1994: The theoretical basis for the parameterization of ice crystal habits: Growth by vapor deposition. *J. Atmos. Sci.*, **51**, 1206–1221.
   - Reisner, J., R. M. Rasmussen, and R. T. Bruintjes, 1998: Explicit forecasting of supercooled liquid water in winter storms using the MM5 forecast model. *Quart. J. Roy. Meteor. Soc.*, **124**, 1071–1107.
   - Cheng, C.-T., W.-C. Wang, and J.-P. Chen, 2010: Simulation of the effects of increasing cloud condensation nuclei on mixed-phase clouds and precipitation of a front system. *Atmos. Res.*, **96**, 461-476.
   - Chen, J.-P., I-C. Tsai and Y.-C. Lin, 2013: A statistical-numerical aerosol parameterization scheme. *Atmos. Chem. Phys.*, **13**, 10483–10504.
   - Chen, J.-P., and T.-C. Tsai, 2016: Triple-moment modal parameterization for the adaptive growth habit of pristine ice crystals. *J. Atmos. Sci.*, **73**, 2105–2122.

3. **PPRO Library Documentation**:  
   See `code/ppro/README.md` for library-level documentation

---

## Version History

- **2025-01**: Added TCWA2 support
- **2024**: Initial PPRO implementation with Zhang21, Thompson, WSM6, NSSL support

