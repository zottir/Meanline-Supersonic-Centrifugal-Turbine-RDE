# Radial Supersonic Turbine Meanline Solver

This repository contains several MATLAB scripts and supporting function folders.

The main script is:

Meanline.m

---

## Overview

This MATLAB code implements a meanline solver for a radial supersonic turbine, including:

- Perfect gas model with temperature-dependent properties  
- Profile losses (Stewart method with Stratford–Beavers boundary layer quantities)  
- Bow shock losses  
- Reflected shock losses  
- Kantrowitz limit verification  
- Conservation checks (mass and rothalpy)  
- Computation of:
  - Total-to-total efficiency  
  - Total-to-static efficiency  
  - Degree of reaction  
  - Power output  

### Visualization Features

The solver can generate:

- T–s diagram  
- Velocity triangles  
- Kantrowitz limit plot  
- Shock pattern visualization  
- Stator–rotor configuration  

The code is intended for preliminary design and analysis of turbomachinery operating in transonic and supersonic regimes.

---

## Model Features

- Thermodynamic properties are computed iteratively from a reference temperature:
  - γ(T)
  - R(T)
  - cp(T)
  - μ(T) (via Sutherland’s law)

The reference temperature is updated until convergence is achieved.

- The number of blades is calculated using the Zweifel criterion.

---

## Physical Checks

The solver verifies:

- Mass flow conservation  
- Rothalpy conservation  
- Non-physical solutions (e.g., radial Mach number inconsistencies)  
- Kantrowitz limit compliance  
- Presence of collective shocks  
- Rotor geometric classification (convergent/divergent)  

---

## Main Input Parameters

### Inlet Conditions

- T0 — Inlet static temperature  
- M0 — Inlet Mach number  
- p0t — Inlet total pressure  
- alpha0 — Inlet flow angle  
- m — Mass flow rate  

### Design Parameters

- phi — Flow coefficient  
- lambda — Work coefficient  
- C_Ft_s, C_Ft_r — Zweifel coefficients  
- Umax — Maximum peripheral speed  
- r0–r3 — Machine radii  

### Flags

- flagPlot → Enables plots  
- flagPlot_Converg → Displays iterative errors  
- flagISO → Forces isentropic flow  
- flagMaps → Disables reflected shocks  

---

## Output Quantities

The solver returns:

- Total-to-total efficiency  
- Total-to-static efficiency  
- Machine power  
- Degree of reaction  
- Stator and rotor losses  
- Mach numbers and velocity triangles  
- Kantrowitz limit verification  

---

## Additional Scripts

### Meanline_Nb_imposed.m

Analogous to Meanline.m, but allows the number of blades to be imposed directly instead of using Zweifel coefficients.

### Meanline_Nb_cp_imposed.m

Similar to the previous script, but also allows the fluid properties (cp, R, γ, μ, etc.) to be imposed.

This removes the iterative loop on the reference temperature, significantly improving convergence speed.

---

## Requirements

- MATLAB  
- Optimization Toolbox (required for fzero)  
- Fluid property database (called via FluidPropInitializer.m)  

---

## Author Notes

The code is structured to:

- Facilitate parametric studies on phi and lambda  
- Analyze the effects of shock waves and viscous losses  
- Support supersonic turbine configurations  

### Important Notes

- If input parameters are poorly tuned, the code may crash due to violation of conservation equations.  
  This occurs particularly when:
  - The mass flow rate is too small  
  - Flow and work coefficients fall outside operational limits  

- When collective shock formation occurs, machine performance may be significantly overestimated.  
  In such conditions, the shock system and reflections are not accurately resolved.

---

## Structure

All supporting functions are organized within the corresponding folders in the repository.
