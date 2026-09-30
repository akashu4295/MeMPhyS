"""Solver configuration schema and CSV helpers used by the Qt interface."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ParameterSpec:
    key: str
    label: str
    default: Any
    kind: str
    section: str
    tooltip: str = ""
    choices: tuple[tuple[str, int], ...] = ()


PARAMETERS = (
    ParameterSpec("domain_dimensions", "Spatial dimensions", 2, "choice", "Meshless discretization", "Choose the dimensionality used by the solver.", (("2D", 2), ("3D", 3))),
    ParameterSpec("poly_deg", "Polynomial degree", 3, "int", "Meshless discretization", "Degree of the appended polynomial basis."),
    ParameterSpec("phs_deg", "PHS degree", 3, "int", "Meshless discretization", "Odd degree of the polyharmonic spline basis."),
    ParameterSpec("cloud_size_multiplier", "Cloud size multiplier", 2, "int", "Meshless discretization", "Scales the local RBF-FD stencil size."),
    ParameterSpec("test_derivative", "Derivative test", 0, "int", "Meshless discretization", "Run the derivative verification test before solving."),
    ParameterSpec("num_time_steps", "Number of time steps", 10000, "int", "Flow and time", "Maximum number of solver time steps."),
    ParameterSpec("write_interval", "Write interval", 10, "int", "Flow and time", "Write solution output every N steps."),
    ParameterSpec("Re", "Reynolds number", 10.0, "float", "Flow and time"),
    ParameterSpec("time_step", "Maximum time step (dt)", 0.1, "float", "Flow and time", "Caps the mesh-based automatic time step."),
    ParameterSpec("courant_number", "Courant number", 0.1, "float", "Flow and time", "Controls the mesh-based automatic time step."),
    ParameterSpec("steady_tolerance", "Steady-state tolerance", 1e-8, "float", "Flow and time"),
    ParameterSpec("compressible_flow", "Enable compressible flow", 0, "bool", "Flow and time"),
    ParameterSpec("fractional_step", "Navier-Stokes method", 1, "choice", "Time integration", choices=(("Fractional step", 1), ("Time implicit", 0))),
    ParameterSpec("iter_momentum", "Momentum iterations", 5, "int", "Time integration"),
    ParameterSpec("iter_timple", "Time-implicit iterations", 1, "int", "Time integration"),
    ParameterSpec("time_scheme", "Time integration scheme", 0, "choice", "Time integration", choices=(("Explicit", 0), ("Implicit", 1))),
    ParameterSpec("theta", "Theta", 0.5, "float", "Time integration", "Implicit time integration weighting; 0.5 is Crank-Nicolson."),
    ParameterSpec("Poisson_solver_type", "Poisson solver", 1, "choice", "Pressure solver", choices=(("Jacobi", 1), ("Gauss-Seidel", 2), ("BiCGStab", 3))),
    ParameterSpec("poisson_solver_tolerance", "Poisson tolerance", 1e-8, "float", "Pressure solver"),
    ParameterSpec("sor_parameter", "SOR parameter", 1.0, "float", "Pressure solver"),
    ParameterSpec("num_vcycles", "Multigrid V-cycles", 10, "int", "Pressure solver"),
    ParameterSpec("num_relax", "Relaxation steps", 100, "int", "Pressure solver"),
    ParameterSpec("num_colors", "Color count", 1, "int", "Pressure solver", "Used for colored Gauss-Seidel."),
    ParameterSpec("facRe", "Reynolds correction factor", 1.0, "float", "Advanced numerics"),
    ParameterSpec("facdt", "Time-step correction factor", 1.0, "float", "Advanced numerics"),
    ParameterSpec("use_hyperviscosity", "Enable hyperviscosity", 0, "bool", "Advanced numerics"),
    ParameterSpec("gamma_hyper", "Hyperviscosity coefficient", 1e-3, "float", "Advanced numerics"),
    ParameterSpec("write_processed_grid_data", "Write processed grid data", 0, "bool", "Advanced numerics"),
    ParameterSpec("restart", "Restart from a previous solution", 0, "bool", "Initial and restart"),
    ParameterSpec("restart_filename", "Restart VTK file", "Field_000029.vtk", "str", "Initial and restart"),
    ParameterSpec("gamma", "Specific heat ratio", 1.4, "float", "Compressible properties"),
    ParameterSpec("R_gas", "Gas constant", 287.0, "float", "Compressible properties"),
    ParameterSpec("Pr", "Prandtl number", 0.71, "float", "Compressible properties"),
    ParameterSpec("Cv", "Specific heat at constant volume", 718.0, "float", "Compressible properties"),
    ParameterSpec("Cp", "Specific heat at constant pressure", 1004.5, "float", "Compressible properties"),
    ParameterSpec("T_ref", "Reference temperature", 300.0, "float", "Compressible properties"),
    ParameterSpec("rho_ref", "Reference density", 1.225, "float", "Compressible properties"),
    ParameterSpec("p_ref", "Reference pressure", 101325.0, "float", "Compressible properties"),
    ParameterSpec("mu_ref", "Reference viscosity", 1.81e-5, "float", "Compressible properties"),
    ParameterSpec("T_sutherland", "Sutherland temperature", 110.4, "float", "Compressible properties"),
    ParameterSpec("viscosity_model", "Viscosity model", 0, "choice", "Compressible properties", choices=(("Constant", 0), ("Sutherland", 1), ("Power law", 2))),
    ParameterSpec("Mach", "Mach number", 0.1, "float", "Compressible properties"),
    ParameterSpec("energy_equation", "Energy equation", 0, "choice", "Compressible properties", choices=(("Isothermal", 0), ("Solve energy", 1))),
)


PARAMETER_BY_KEY = {spec.key: spec for spec in PARAMETERS}


def default_parameters() -> dict[str, Any]:
    return {spec.key: spec.default for spec in PARAMETERS}


def read_parameter_csv(path: str | Path = "flow_parameters.csv") -> dict[str, Any]:
    values = default_parameters()
    path = Path(path)
    if not path.is_file():
        return values

    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.reader(stream):
            if len(row) < 2:
                continue
            key, raw_value = row[0].strip(), row[1].strip()
            spec = PARAMETER_BY_KEY.get(key)
            if spec is None:
                continue
            try:
                if spec.kind == "bool":
                    values[key] = int(raw_value != "0")
                elif spec.kind in {"int", "choice"}:
                    values[key] = int(raw_value)
                elif spec.kind == "float":
                    values[key] = float(raw_value)
                else:
                    values[key] = raw_value
            except ValueError:
                continue
    return values


def write_parameter_csv(values: dict[str, Any], path: str | Path = "flow_parameters.csv") -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        for spec in PARAMETERS:
            if spec.section == "Compressible properties" and not bool(values.get("compressible_flow", 0)):
                continue
            value = values.get(spec.key, spec.default)
            if spec.kind == "bool":
                value = int(bool(value))
            elif spec.kind == "float":
                value = format(float(value), ".12g")
            writer.writerow((spec.key, value))


def read_grid_csv(path: str | Path = "grid_filenames.csv") -> list[str]:
    path = Path(path)
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.reader(stream))
    if not rows:
        return []
    if rows[0] and rows[0][0].strip() == "num_levels":
        rows = rows[1:]
    return [row[0].strip() for row in rows if row and row[0].strip()]


def write_grid_csv(mesh_files: list[str], path: str | Path = "grid_filenames.csv") -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("num_levels", len(mesh_files)))
        writer.writerows((mesh_file,) for mesh_file in mesh_files)