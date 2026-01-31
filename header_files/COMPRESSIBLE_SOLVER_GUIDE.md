# Compressible Navier-Stokes Solver for Low Mach Number Flows
## Implementation Guide

## 1. ADDITIONAL FIELD VARIABLES REQUIRED

Add these to your `FieldVariables` structure:

### Primary Thermodynamic Variables
```c
// Density
double *rho;           // Current density
double *rho_old;       // Previous time step density  
double *rho_new;       // Intermediate density

// Temperature
double *T;             // Current temperature
double *T_old;         // Previous time step temperature
double *T_new;         // Intermediate temperature

// Internal Energy
double *e;             // Internal energy (e = cv*T)
double *e_old;         // Previous time step internal energy
```

### Gradient Variables
```c
// Density gradients
double *drhodx, *drhody, *drhodz;

// Temperature gradients  
double *dTdx, *dTdy, *dTdz;

// Internal energy gradients
double *dedx, *dedy, *dedz;
```

### Viscous Stress Tensor Components
```c
// Dynamic viscosity (temperature-dependent)
double *mu;

// Thermal conductivity
double *kappa;

// Normal stresses
double *tau_xx, *tau_yy, *tau_zz;

// Shear stresses
double *tau_xy, *tau_xz, *tau_yz;

// Divergence of stress tensor
double *div_tau_x, *div_tau_y, *div_tau_z;
```

### Source Terms
```c
double *Q_visc;        // Viscous dissipation (heating)
double *Q_source;      // External heat source term
```

**Total new arrays needed: ~30 additional double arrays**

---

## 2. ADDITIONAL PARAMETERS

Add these to your `parameters` structure:

```c
// Gas properties
double gamma;          // Ratio of specific heats (1.4 for air)
double R_gas;          // Gas constant (287 J/(kg·K) for air)
double Pr;             // Prandtl number (0.71 for air)
double cv;             // Specific heat at constant volume
double cp;             // Specific heat at constant pressure (cp = gamma*cv/(gamma-1))

// Reference conditions
double T_ref;          // Reference temperature (e.g., 300 K)
double rho_ref;        // Reference density (e.g., 1.225 kg/m³)
double p_ref;          // Reference pressure (e.g., 101325 Pa)
double mu_ref;         // Reference viscosity (e.g., 1.81e-5 Pa·s)

// Viscosity model
double T_sutherland;   // Sutherland's constant (110.4 K for air)
int viscosity_model;   // 0=constant, 1=Sutherland, 2=power-law

// Mach number
double Mach;           // Reference Mach number

// Solver options
int energy_equation;   // 0=isothermal, 1=solve energy equation
```

---

## 3. MAIN DIFFERENCES FROM INCOMPRESSIBLE SOLVER

### Algorithm Structure
**Incompressible:**
1. Solve momentum for u*
2. Solve pressure Poisson
3. Correct velocity

**Compressible (Low Mach):**
1. Solve continuity for ρ
2. Solve momentum for u* (with variable ρ)
3. Solve energy for T
4. Solve pressure correction (low Mach approximation)
5. Correct velocity

### Key Equations

#### Continuity Equation
```
∂ρ/∂t + ∇·(ρu) = 0
```
- Solves for density evolution
- Density depends on temperature via ideal gas law

#### Momentum Equations
```
∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·τ
```
- Similar to incompressible but with variable density
- Stress tensor depends on temperature-dependent viscosity

#### Energy Equation
```
∂(ρe)/∂t + ∇·(ρeu) = -p∇·u + ∇·(k∇T) + Φ + Q
```
Where:
- e = internal energy = cv*T
- k = thermal conductivity
- Φ = viscous dissipation
- Q = heat source

#### Pressure (Low Mach Approximation)
```
p = p_ref + p'    where p' << p_ref
∇²p' = ∇·u/Δt
```
- Thermodynamic pressure p_ref is nearly constant
- Only solve for dynamic pressure correction p'

---

## 4. VISCOSITY MODELS

### Constant Viscosity
```c
μ = μ_ref
```

### Sutherland's Law (Most Accurate for Air)
```c
μ(T) = μ_ref * (T/T_ref)^(3/2) * (T_ref + T_s)/(T + T_s)
```
Where T_s = 110.4 K for air

### Power Law (Simple Approximation)
```c
μ(T) = μ_ref * (T/T_ref)^n
```
Where n ≈ 0.7 for air

---

## 5. BOUNDARY CONDITIONS

### Additional BC Types Needed

```c
typedef enum {
    BC_VELOCITY_INLET,      // Specify: u, v, w, T, (p or ρ)
    BC_PRESSURE_OUTLET,     // Specify: p, extrapolate T and u
    BC_WALL,                // No-slip: u=v=w=0
    BC_ISOTHERMAL_WALL,     // Wall with fixed T
    BC_ADIABATIC_WALL,      // Wall with ∂T/∂n = 0
    BC_SYMMETRY,            // Zero normal velocity and gradients
} BoundaryConditionType;
```

### BC Structure Updates
```c
typedef struct {
    BoundaryConditionType type;
    double u, v, w;         // Velocity components
    double p;               // Pressure
    double T;               // Temperature
    double rho;             // Density (alternative to p)
} BoundaryCondition;
```

---

## 6. INITIALIZATION

### Initial Conditions Required
```c
// Set initial values for all nodes
for (int i = 0; i < num_nodes; i++) {
    field->T[i] = T_initial;                    // e.g., 300 K
    field->p[i] = p_ref;                        // e.g., 101325 Pa
    field->rho[i] = p_ref / (R_gas * T[i]);    // Ideal gas law
    field->e[i] = cv * T[i];                    // Internal energy
    field->u[i] = u_initial;                    // Initial velocity
    field->v[i] = v_initial;
    field->w[i] = 0.0;
}
```

---

## 7. TIME STEP CALCULATION

Compressible flow has **three CFL conditions**:

```c
// Convective CFL
dt_conv = Δx / (|u| + |v| + |w| + a)

// Diffusive CFL  
dt_diff = Δx² / (2ν)

// Acoustic CFL
dt_acoustic = Δx / a

dt = min(dt_conv, dt_diff, dt_acoustic) * CFL
```

Where `a = sqrt(γRT)` is the speed of sound.

**Important:** For low Mach flows, acoustic CFL can be relaxed since we're using pressure splitting.

---

## 8. TYPICAL PARAMETER VALUES FOR AIR AT STANDARD CONDITIONS

```c
parameters.gamma = 1.4;              // Ratio of specific heats
parameters.R_gas = 287.0;            // J/(kg·K)
parameters.Pr = 0.71;                // Prandtl number
parameters.T_ref = 300.0;            // K
parameters.p_ref = 101325.0;         // Pa
parameters.rho_ref = 1.177;          // kg/m³
parameters.mu_ref = 1.85e-5;         // Pa·s
parameters.T_sutherland = 110.4;     // K
parameters.cv = 717.5;               // J/(kg·K)
parameters.cp = 1004.5;              // J/(kg·K)
```

---

## 9. MEMORY REQUIREMENTS COMPARISON

**Incompressible solver:**
- ~15 arrays per field variable (u, v, w, p + helpers)
- Memory ≈ 15 × N × sizeof(double)

**Compressible solver:**
- ~45 arrays (u, v, w, p, ρ, T, e + stresses + gradients)
- Memory ≈ 45 × N × sizeof(double)

**Memory increase: ~3× compared to incompressible**

---

## 10. VALIDATION TEST CASES

### Test 1: Isothermal Low Mach Channel Flow
- Set `energy_equation = 0` (skip energy solver)
- Fix T = constant everywhere
- Should match incompressible results at low Re

### Test 2: Natural Convection (Boussinesq)
- Add buoyancy force: `F = -ρg(1 - β(T-T_ref))`
- Solve energy equation with fixed walls
- Classic test: differentially heated cavity

### Test 3: Heated Channel Flow
- Hot wall at inlet, cold wall elsewhere
- Develop temperature profile
- Check energy balance

### Test 4: Compressible Boundary Layer
- High-speed flow over plate
- Temperature rise due to viscous dissipation
- Compare with analytical solutions

---

## 11. DEBUGGING TIPS

### Check Conservation
```c
// Mass conservation
double mass_total = 0.0;
for (int i = 0; i < N; i++) {
    mass_total += field->rho[i] * volume[i];
}
printf("Total mass: %e (should be constant)\n", mass_total);

// Energy conservation (if adiabatic)
double energy_total = 0.0;
for (int i = 0; i < N; i++) {
    double KE = 0.5 * (u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
    energy_total += (field->e[i] + KE) * field->rho[i] * volume[i];
}
printf("Total energy: %e (should be constant)\n", energy_total);
```

### Check Thermodynamic Consistency
```c
// Verify ideal gas law: p = ρRT
for (int i = 0; i < N; i++) {
    double p_calc = field->rho[i] * R_gas * field->T[i];
    if (fabs(p_calc - p_ref) / p_ref > 0.01) {
        printf("WARNING: Thermodynamic inconsistency at node %d\n", i);
    }
}
```

---

## 12. PERFORMANCE CONSIDERATIONS

### GPU Optimization
- Reuse arrays when possible (e.g., dpdx, dpdy for multiple purposes)
- Fuse loops where dependencies allow
- Group related operations in same kernel

### Convergence Acceleration
- Use multigrid for pressure (already implemented)
- Consider implicit treatment of diffusion terms
- Use higher-order time integration for smooth problems

---

## 13. EXTENSIONS

### Adding Real Gas Effects
- Replace ideal gas law with real gas EOS (van der Waals, Redlich-Kwong)
- Modify pressure calculation

### Adding Turbulence
- Add RANS models (k-ε, k-ω)
- Additional equations for turbulent kinetic energy and dissipation
- Modify viscosity: μ_total = μ + μ_turbulent

### Adding Chemical Reactions
- Add species transport equations
- Include reaction source terms in energy equation
- Variable molecular weight and specific heats

---

## 14. QUICK START GUIDE

1. **Extend structures**
   - Add 30 new double arrays to FieldVariables
   - Add 15 parameters to parameters structure

2. **Initialize**
   - Set T, ρ, p to physical values
   - Calculate derived quantities (e, μ, κ)

3. **Replace solver**
   - Use `compressible_solver_explicit()` instead of `fractional_step_explicit()`
   - Keep same multigrid pressure solver

4. **Test incrementally**
   - Start with isothermal (energy_equation = 0)
   - Add energy equation once isothermal works
   - Add temperature-dependent viscosity last

5. **Validate**
   - Check mass, momentum, energy conservation
   - Compare with known solutions
   - Verify CFL stability limits

---

## SUMMARY

The compressible solver adds:
- ✓ Density as primary variable (via continuity equation)
- ✓ Temperature evolution (via energy equation)  
- ✓ Variable viscosity (temperature-dependent)
- ✓ Viscous dissipation (converts kinetic to thermal energy)
- ✓ Thermal conductivity (heat diffusion)
- ✓ Compressibility effects (via ideal gas law)

While maintaining:
- ✓ Same meshless RBF framework
- ✓ Same multigrid Poisson solver
- ✓ Same GPU acceleration approach
- ✓ Similar code structure and style
