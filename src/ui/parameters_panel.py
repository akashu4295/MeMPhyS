"""
Parameters panel construction for MeMPhyS GUI

Creates the left panel with solver parameters, mesh file selection,
and run controls.
"""

import dearpygui.dearpygui as dpg

from src.config import (
    BASE_PARAMETERS,
    IMPLICIT_PARAMETERS,
    MULTIGRID_PARAMETERS,
    SOLVER_METHODS,
    Poisson_SOLVER_METHODS,
    DEFAULT_SOLVER_METHOD,
    MAX_MESH_LEVELS,
    LEFT_PANEL_WIDTH,
    PARAMETER_INPUT_WIDTH,
    MESH_PATH_INPUT_WIDTH,
    RUN_BUTTON_WIDTH,
    RUN_BUTTON_HEIGHT,
    COLORS,
    FILE_EXTENSIONS,
    FILE_DIALOG_WIDTH,
    FILE_DIALOG_HEIGHT,
)
from src.callbacks import (
    run_solver_callback,
    validate_numeric_input,
    show_implicit_callback,
    show_multigrid_callback,
    update_mesh_inputs_callback,
    open_file_dialog_callback,
    select_mesh_file_callback,
    launch_gmsh_callback,
    set_mesh_from_geometry_callback,
    browse_geometry_file_callback,
    select_geometry_file_callback,
    show_bc_window_callback,
)


# def create_parameters_panel(themes: dict) -> int:
#     """
#     Create the left panel with parameters and controls
    
#     Args:
#         themes: Dictionary of theme IDs
    
#     Returns:
#         Child window tag/ID
#     """
#     with dpg.child_window(
#         width=LEFT_PANEL_WIDTH,
#         border=True,
#         tag="parameters_panel"
#     ) as panel:
                

#         # Header
#         dpg.add_text("Input & Parameters", color=COLORS["header"])
#         dpg.add_spacer(height=5)
#         dpg.add_separator()
        
#         # Solver Method Selection
#         _create_solver_method_section()
        
#         dpg.add_separator()
        
#         # Flow Parameters
#         _create_flow_parameters_section()
#         _create_poisson_solver_method_section()
#         # Implicit Parameters (hidden by default)
#         _create_implicit_parameters_section()
        
#         dpg.add_separator()
        
#         # Header
#         dpg.add_text("Geometry and Grids", color=COLORS["header"])

#         # Multigrid Section
#         _create_multigrid_section()
        
#         # Solver Method Selection
#         _create_geometry_section(themes)

#         # Mesh File Selection
#         _create_mesh_files_section(themes)
        
#         dpg.add_separator()
        
#         # Initialization File Selection
#         _create_init_file_section()
#         dpg.add_spacer(height=5)
#         _create_restart_file_section()
#         dpg.add_spacer(height=5)
#         dpg.add_separator()
#         # Run Button
#         _create_run_button(themes)
    
#     return panel


# def _create_solver_method_section():
#     """Create solver method selection combo box"""
#     with dpg.group(horizontal=True):
#         dpg.add_text("Solver Method", color=COLORS["header"])
#         dpg.add_combo(
#             items=SOLVER_METHODS,
#             default_value=DEFAULT_SOLVER_METHOD,
#             tag="solver_method",
#             width=PARAMETER_INPUT_WIDTH+80,
#             callback=show_implicit_callback
#         )

#     with dpg.tooltip("solver_method"):
#         dpg.add_text(
#             "Select the solver method for the Navier-Stokes equations",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("Implicit solvers are often more stable for larger time steps, but require solving a linear system at each step")
#         dpg.add_text("Explicit solvers can be faster for smaller meshes and time steps, but may require more careful tuning of parameters for stability")

# def _create_poisson_solver_method_section():
#     """Create solver method selection combo box"""
#     with dpg.group(
#         horizontal=True,
#         tag="poisson_solver_method_group"
#     ):
#         dpg.add_combo(
#             items=Poisson_SOLVER_METHODS,
#             default_value=Poisson_SOLVER_METHODS[0],
#             tag="poisson_solver_method",
#             width=PARAMETER_INPUT_WIDTH,
#         )
#         dpg.add_text("Poisson Solver Type")
    
#     with dpg.tooltip("poisson_solver_method_group"):
#         dpg.add_text(
#             "Select the solver type for the Poisson equation",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("The choice of Poisson solver can affect convergence and performance, especially for larger meshes")
#         dpg.add_text("Multigrid solvers are often more efficient for large problems, while direct solvers can be faster for smaller meshes")
#         dpg.add_text("If using multigrid, make sure to enable it in the multigrid section below and provide the appropriate mesh files for each level")

# def _create_flow_parameters_section():
#     """Create flow parameters input section"""
#     dpg.add_text("Flow Parameters", color=COLORS["header"])
    
#     for pname, pval in BASE_PARAMETERS.items():
#         tag = f"param_{pname}"
#         dpg.add_input_text(
#             label=pname,
#             tag=tag,
#             default_value=str(pval),
#             width=PARAMETER_INPUT_WIDTH,
#             callback=validate_numeric_input,
#             on_enter=True
#         )


# def _create_implicit_parameters_section():
#     """Create implicit solver parameters (hidden by default)"""
#     for pname, pval in IMPLICIT_PARAMETERS.items():
#         tag = f"param_{pname}"
#         dpg.add_input_text(
#             label=pname,
#             tag=tag,
#             default_value=str(pval),
#             width=PARAMETER_INPUT_WIDTH,
#             show=False,
#             callback=validate_numeric_input,
#             on_enter=True
#         )


# def _create_multigrid_section():
#     """Create multigrid enable checkbox and parameters"""
#     # Multigrid toggle
#     dpg.add_checkbox(
#         label="Enable Multigrid",
#         tag="multigrid_toggle",
#         callback=show_multigrid_callback
#     )
    
#     with dpg.tooltip("multigrid_toggle"):
#         dpg.add_text(
#             "Check to enable multigrid solver with multiple mesh levels",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("When enabled, you can specify parameters for each mesh level and provide multiple mesh files (finest to coarsest)")
#         dpg.add_text("Make sure to provide the correct number of mesh files corresponding to the number of mesh levels specified in the multigrid parameters below")


#     # Multigrid parameters group (hidden by default)
#     with dpg.group(
#         horizontal=False,
#         tag="multigrid_parameters_section",
#         show=False
#     ):
#         dpg.add_text("Multigrid Parameters", color=COLORS["header"])
        
#         # Multigrid parameter inputs
#         for pname, pval in MULTIGRID_PARAMETERS.items():
#             tag = f"param_{pname}"
#             dpg.add_input_text(
#                 label=pname,
#                 tag=tag,
#                 default_value=str(pval),
#                 width=PARAMETER_INPUT_WIDTH,
#                 callback=validate_numeric_input,
#                 on_enter=True
#             )
        
#         # Number of mesh levels
#         dpg.add_input_int(
#             label="Number of mesh levels",
#             default_value=1,
#             min_value=1,
#             max_value=MAX_MESH_LEVELS,
#             width=PARAMETER_INPUT_WIDTH,
#             tag="num_mesh_levels",
#             callback=update_mesh_inputs_callback
#         )


# def _create_geometry_section(themes: dict):
#     """Create geometry/Gmsh section"""
#     dpg.add_text("Create geometry in Gmsh, use a .geo file or browse for existing .msh file", color=COLORS["success"], wrap=LEFT_PANEL_WIDTH - 20)
    
#     with dpg.group(horizontal=True):
#         new_geo_btn = dpg.add_button(
#             label="Open Gmsh GUI",
#             callback=launch_gmsh_callback,
#             width=120
#         )
        
#         open_geo_btn = dpg.add_button(
#             label="Use .geo File",
#             callback=browse_geometry_file_callback,
#             width=120
#         )

#         # Add tooltip
#         with dpg.tooltip(new_geo_btn):
#             dpg.add_text(
#                 "Click to open the Gmsh GUI for geometry creation",
#                 color=COLORS["info"]
#             )
#             dpg.add_separator()
#             dpg.add_text("Make sure to:")
#             dpg.add_text("  1. Name the boundaries in Gmsh for correct BC association")
#             dpg.add_text("  2. Save the mesh file with .msh (ASCII version 2) extension")
#             dpg.add_text("  3. Load the .msh file in the mesh section below")

#         with dpg.tooltip(open_geo_btn):
#             dpg.add_text(
#                 "Click to select a .geo file to generate the mesh from",
#                 color=COLORS["info"]
#             )
#             dpg.add_separator()
#             dpg.add_text("Make sure to:")
#             dpg.add_text("  1. Name the boundaries in Gmsh for correct BC association")
#             dpg.add_text("  2. Save the .geo file and select it here")
#             dpg.add_text("  3. The mesh will be generated automatically and loaded in the mesh section below")

#         # Apply themes
#         if "button_secondary" in themes:
#             dpg.bind_item_theme(new_geo_btn, themes["button_secondary"])
#             dpg.bind_item_theme(open_geo_btn, themes["button_secondary"])
        
    
#     # File dialog for .geo files
#     dpg.add_file_dialog(
#         directory_selector=False,
#         tag="file_dialog_geometry",
#         callback=select_geometry_file_callback,
#         show=False,
#         width=FILE_DIALOG_WIDTH,
#         height=FILE_DIALOG_HEIGHT
#     )
#     dpg.add_file_extension(".geo", parent="file_dialog_geometry", color=(150, 255, 150, 255))
#     dpg.add_file_extension(".geo_unrolled", parent="file_dialog_geometry", color=(150, 255, 150, 255))
    
#     dpg.add_spacer(height=5)


# def _create_mesh_files_section(themes: dict):
#     """Create mesh file selection section"""
#     with dpg.group(tag="multigrid_mesh_section"):
#         # dpg.add_text("Choose Mesh Files", color=COLORS["success"])
#         dpg.add_text("(finest to coarsest if multigrid)", color=COLORS["success"])
        
#         # Create mesh file inputs for all levels
#         for i in range(MAX_MESH_LEVELS):
#             idx = i + 1
#             show_initial = (i == 0)  # Only show first one initially
            
#             # Group for each mesh level
#             with dpg.group(
#                 horizontal=True,
#                 show=show_initial,
#                 tag=f"mesh_group_{idx}"
#             ):
#                 dpg.add_input_text(
#                     tag=f"mesh_file_{idx}",
#                     hint=f"Mesh file {idx} (choose or type path)",
#                     width=MESH_PATH_INPUT_WIDTH,
#                     show=show_initial
#                 )
                
#                 dpg.add_button(
#                     label="Browse",
#                     tag=f"browse_{idx}",
#                     callback=open_file_dialog_callback,
#                     user_data=idx,
#                     show=show_initial
#                 )
        
#         # Create file dialogs for mesh files
#         for i in range(MAX_MESH_LEVELS):
#             idx = i + 1
#             dialog_tag = f"file_dialog_{idx}"
#             input_tag = f"mesh_file_{idx}"
            
#             dpg.add_file_dialog(
#                 directory_selector=False,
#                 tag=dialog_tag,
#                 user_data=input_tag,
#                 callback=select_mesh_file_callback,
#                 show=False,
#                 width=FILE_DIALOG_WIDTH,
#                 height=FILE_DIALOG_HEIGHT
#             )
            
#             # Add file extension filter
#             if ".msh" in FILE_EXTENSIONS:
#                 dpg.add_file_extension(
#                     ".msh",
#                     parent=dialog_tag,
#                     color=FILE_EXTENSIONS[".msh"]
#                 )
    
#     dpg.add_spacer(height=5)
#     associate_bc_btn = dpg.add_button(
#         label="Set Boundary Conditions",
#         callback=show_bc_window_callback,
#         width=250
#     )

#     with dpg.tooltip(associate_bc_btn):
#         dpg.add_text(
#             "Click to open the boundary condition association window",
#             color=COLORS["info"]
#         )
        
#     # Apply themes
#     if "button_secondary" in themes:
#         dpg.bind_item_theme(associate_bc_btn, themes["button_secondary"])
    


# def _create_init_file_section():
#     """Create initialization file selection section"""
#     # dpg.add_text("Choose the initialisation file if you need non-zero initial conditions", color=COLORS["subheader"], wrap=LEFT_PANEL_WIDTH - 20)
#     dpg.add_checkbox(
#         label="Initialise with user defined conditions",
#         tag="init_toggle",
#         callback=show_init_callback
#     )

#     with dpg.tooltip("init_toggle"):
#         dpg.add_text(
#             "Check to provide an initialization file with user defined initial conditions",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("The initialization file should be a .c file that defines the function:")
#         dpg.add_text("An example template can be found in the repository named as 'init.c'")
    
#     with dpg.group(horizontal=True, show=False, tag="init_group"):
#         dpg.add_input_text(
#             hint="Path to Initialisation file",
#             tag="init_path",
#             width=MESH_PATH_INPUT_WIDTH,
#             show=True
#         )
        
#         dpg.add_button(
#             label="Browse",
#             tag="init_browse",
#             callback=open_file_dialog_callback,
#             user_data="init",
#             show=True
#         )
    
#     # File dialog for init file
#     dpg.add_file_dialog(
#         directory_selector=False,
#         tag="file_dialog_init",
#         user_data="init_path",
#         callback=select_mesh_file_callback,
#         show=False,
#         width=FILE_DIALOG_WIDTH,
#         height=FILE_DIALOG_HEIGHT
#     )
    
#     # Add .c file extension filter
#     if ".c" in FILE_EXTENSIONS:
#         dpg.add_file_extension(
#             ".c",
#             parent="file_dialog_init",
#             color=FILE_EXTENSIONS[".c"]
#         )

# def _create_restart_file_section():
#     """Create Restart file selection section"""
#     # dpg.add_text("Choose the restart file to continue from a previous run", color=COLORS["subheader"], wrap=LEFT_PANEL_WIDTH - 20)
#     dpg.add_checkbox(
#         label="Restart from previous run",
#         tag="restart_toggle",
#         callback=show_restart_callback
#     )

#     with dpg.tooltip("restart_toggle"):
#         dpg.add_text(
#             "Check to provide a restart file from a previous run",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("The restart file should be a .csv file generated by a previous run of the solver with the same mesh file")

#     with dpg.group(horizontal=True, show=False, tag="restart_group"):
#         dpg.add_input_text(
#             hint="Path to Restart file",
#             tag="restart_path",
#             width=MESH_PATH_INPUT_WIDTH,
#             show=True
#         )
        
#         dpg.add_button(
#             label="Browse",
#             tag="restart_browse",
#             callback=open_file_dialog_callback,
#             user_data="restart",
#             show=True
#         )
    
#     # File dialog for restart file
#     dpg.add_file_dialog(
#         directory_selector=False,
#         tag="file_dialog_restart",
#         user_data="restart_path",
#         callback=select_mesh_file_callback,
#         show=False,
#         width=FILE_DIALOG_WIDTH,
#         height=FILE_DIALOG_HEIGHT
#     )
    
#     # Add .csv file extension filter
#     if ".csv" in FILE_EXTENSIONS:
#         dpg.add_file_extension(
#             ".csv",
#             parent="file_dialog_restart",
#             color=FILE_EXTENSIONS[".csv"]
#         )

# def show_restart_callback(sender, app_data, user_data):
#     """Show or hide restart file input based on checkbox"""
#     show = app_data  # Checkbox value (True/False)
#     dpg.configure_item("restart_group", show=show)

# def show_init_callback(sender, app_data, user_data):
#     """Show or hide initialization file input based on checkbox"""
#     show = app_data  # Checkbox value (True/False)
#     dpg.configure_item("init_group", show=show)

# def _create_run_button(themes: dict):
#     """Create the main run button"""
#     dpg.add_spacer(height=5)
#     dpg.add_text("SOLVE!", color=COLORS["header"], wrap=LEFT_PANEL_WIDTH - 20)
#     compile_gpu_checkbox = dpg.add_checkbox(
#         label="Run on GPU (if supported)",
#         tag="gpu_toggle"
#     )
#     run_btn = dpg.add_button(
#         label="Compile and Run Solver",
#         width=RUN_BUTTON_WIDTH,
#         height=RUN_BUTTON_HEIGHT,
#         callback=run_solver_callback,
#         tag="run_solver_button"
#     )
    
#     # Add tooltip
#     with dpg.tooltip(run_btn):
#         dpg.add_text(
#             "Click to compile and execute the solver",
#             color=COLORS["info"]
#         )
#         dpg.add_separator()
#         dpg.add_text("This will:")
#         dpg.add_text("  1. Write configuration files")
#         dpg.add_text("  2. Compile all C source files with gcc (CPU only)")
#         dpg.add_text("  3. Run the solver executable")
    
#     # Apply button theme
#     if "button" in themes:
#         dpg.bind_item_theme(run_btn, themes["button"])

def create_parameters_panel(themes: dict) -> int:
    with dpg.child_window(
        width=LEFT_PANEL_WIDTH,
        border=True,
        tag="parameters_panel"
    ) as panel:

        _panel_header("Geometry and Grids")
        _create_geometry_section(themes)
        _create_multigrid_section()
        _create_mesh_files_section(themes)   # BC button lives at bottom of this

        _section_gap()
        _panel_header("Solver Configuration")
        _create_solver_method_section(themes)
        _create_poisson_solver_method_section()

        _section_gap()
        _panel_header("Flow Parameters")
        _create_flow_parameters_section()
        _create_implicit_parameters_section()

        _section_gap()
        _panel_header("Initial and Restart Conditions")
        _create_init_file_section()
        dpg.add_spacer(height=6)
        _create_restart_file_section()

        _section_gap()
        _create_run_button(themes)

    return panel

# ── Visual helpers ─────────────────────────────────────────────────────────────

def _panel_header(title: str):
    """Visually distinct section header with accent bar effect."""
    dpg.add_spacer(height=2)
    dpg.add_text(title, color=COLORS["header"])
    dpg.add_separator()
    dpg.add_spacer(height=2)


def _section_gap():
    """Consistent vertical breathing room between sections."""
    dpg.add_spacer(height=8)


def _labeled_row(label: str, widget_fn, *args, **kwargs):
    """Helper to put a fixed-width label beside any widget."""
    with dpg.group(horizontal=True):
        dpg.add_text(f"{label:<22}")   # pad to align widgets
        return widget_fn(*args, **kwargs)


# ── Sections ───────────────────────────────────────────────────────────────────

def _create_solver_method_section(themes: dict):
    with dpg.group(horizontal=True):
        dpg.add_text("Method", color=COLORS["subheader"])
        dpg.add_combo(
            items=SOLVER_METHODS,
            default_value=DEFAULT_SOLVER_METHOD,
            tag="solver_method",
            width=PARAMETER_INPUT_WIDTH + 80,
            callback=show_implicit_callback
        )
    with dpg.tooltip("solver_method"):
        dpg.add_text("Select the solver method for the Navier-Stokes equations",
                     color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("Implicit: more stable for large time steps, solves a linear system each step.")
        dpg.add_text("Explicit: faster for small meshes/time steps, needs careful stability tuning.")


def _create_poisson_solver_method_section():
    with dpg.group(horizontal=True, tag="poisson_solver_method_group"):
        dpg.add_text("Poisson Solver", color=COLORS["subheader"])
        dpg.add_combo(
            items=Poisson_SOLVER_METHODS,
            default_value=Poisson_SOLVER_METHODS[0],
            tag="poisson_solver_method",
            width=PARAMETER_INPUT_WIDTH+25,
        )
    with dpg.tooltip("poisson_solver_method_group"):
        dpg.add_text("Select the solver type for the Poisson equation",
                     color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("Multigrid: efficient for large problems.")
        dpg.add_text("Direct: faster for smaller meshes.")


def _create_flow_parameters_section():
    for pname, pval in BASE_PARAMETERS.items():
        tag = f"param_{pname}"
        dpg.add_input_text(
            label=pname,
            tag=tag,
            default_value=str(pval),
            width=PARAMETER_INPUT_WIDTH,
            callback=validate_numeric_input,
            on_enter=True
        )


def _create_implicit_parameters_section():
    for pname, pval in IMPLICIT_PARAMETERS.items():
        tag = f"param_{pname}"
        dpg.add_input_text(
            label=pname,
            tag=tag,
            default_value=str(pval),
            width=PARAMETER_INPUT_WIDTH,
            show=False,
            callback=validate_numeric_input,
            on_enter=True
        )


def _create_multigrid_section():
    dpg.add_checkbox(
        label="Enable Multigrid",
        tag="multigrid_toggle",
        callback=show_multigrid_callback
    )
    with dpg.tooltip("multigrid_toggle"):
        dpg.add_text("Enable multigrid solver with multiple mesh levels",
                     color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("Provide mesh files finest -> coarsest.")
        dpg.add_text("Number of files must match the mesh levels specified below.")

    with dpg.group(horizontal=False, tag="multigrid_parameters_section", show=False):
        dpg.add_spacer(height=4)
        dpg.add_text("Multigrid Parameters", color=COLORS["subheader"])
        for pname, pval in MULTIGRID_PARAMETERS.items():
            dpg.add_input_text(
                label=pname,
                tag=f"param_{pname}",
                default_value=str(pval),
                width=PARAMETER_INPUT_WIDTH,
                callback=validate_numeric_input,
                on_enter=True
            )
        dpg.add_input_int(
            label="Mesh levels",
            default_value=1,
            min_value=1,
            max_value=MAX_MESH_LEVELS,
            width=PARAMETER_INPUT_WIDTH,
            tag="num_mesh_levels",
            callback=update_mesh_inputs_callback
        )
        dpg.add_spacer(height=4)


def _create_geometry_section(themes: dict):
    dpg.add_text(
        "Create geometry in Gmsh, use a .geo file, or browse for an existing .msh file",
        color=COLORS["success"],
        wrap=LEFT_PANEL_WIDTH - 20
    )
    dpg.add_spacer(height=4)

    with dpg.group(horizontal=True):
        new_geo_btn = dpg.add_button(
            label="Open Gmsh GUI",
            callback=launch_gmsh_callback,
            width=130
        )
        open_geo_btn = dpg.add_button(
            label="Use .geo File",
            callback=browse_geometry_file_callback,
            width=130
        )

    with dpg.tooltip(new_geo_btn):
        dpg.add_text("Open the Gmsh GUI for geometry creation", color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("1. Name boundaries for correct BC association")
        dpg.add_text("2. Save mesh as .msh (ASCII v2)")
        dpg.add_text("3. Load the .msh file below")

    with dpg.tooltip(open_geo_btn):
        dpg.add_text("Select a .geo file to generate the mesh from", color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("1. Name boundaries for correct BC association")
        dpg.add_text("2. Save and select the .geo file here")
        dpg.add_text("3. Mesh is generated and loaded automatically")

    if "button_secondary" in themes:
        dpg.bind_item_theme(new_geo_btn, themes["button_secondary"])
        dpg.bind_item_theme(open_geo_btn, themes["button_secondary"])

    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_geometry",
        callback=select_geometry_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT
    )
    dpg.add_file_extension(".geo",          parent="file_dialog_geometry", color=(150, 255, 150, 255))
    dpg.add_file_extension(".geo_unrolled", parent="file_dialog_geometry", color=(150, 255, 150, 255))
    dpg.add_spacer(height=4)


def _create_mesh_files_section(themes: dict):
    with dpg.group(tag="multigrid_mesh_section"):
        dpg.add_text("Mesh files  (finest to coarsest)",
                     color=COLORS["success"], tag="multigrid_hint_text", show=False)
        dpg.add_spacer(height=3)

        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            show_initial = (i == 0)
            with dpg.group(horizontal=True, show=show_initial, tag=f"mesh_group_{idx}"):
                dpg.add_input_text(
                    tag=f"mesh_file_{idx}",
                    hint=f"Mesh {idx}",
                    width=MESH_PATH_INPUT_WIDTH,
                    show=show_initial
                )
                dpg.add_button(
                    label="Browse",
                    tag=f"browse_{idx}",
                    callback=open_file_dialog_callback,
                    user_data=idx,
                    show=show_initial
                )

        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            dialog_tag = f"file_dialog_{idx}"
            dpg.add_file_dialog(
                directory_selector=False,
                tag=dialog_tag,
                user_data=f"mesh_file_{idx}",
                callback=select_mesh_file_callback,
                show=False,
                width=FILE_DIALOG_WIDTH,
                height=FILE_DIALOG_HEIGHT
            )
            if ".msh" in FILE_EXTENSIONS:
                dpg.add_file_extension(".msh", parent=dialog_tag, color=FILE_EXTENSIONS[".msh"])

    dpg.add_spacer(height=6)
    associate_bc_btn = dpg.add_button(
        label="Set Boundary Conditions",
        callback=show_bc_window_callback,
        width=270
    )
    with dpg.tooltip(associate_bc_btn):
        dpg.add_text("Open the boundary condition association window",
                     color=COLORS["info"])
    if "button_secondary" in themes:
        dpg.bind_item_theme(associate_bc_btn, themes["button_secondary"])


def _create_init_file_section():
    dpg.add_checkbox(
        label="Initialise with user-defined conditions",
        tag="init_toggle",
        callback=show_init_callback
    )
    with dpg.tooltip("init_toggle"):
        dpg.add_text("Provide a .c file defining custom initial conditions",
                     color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("See 'init.c' in the repository for a template.")

    with dpg.group(horizontal=True, show=False, tag="init_group"):
        dpg.add_input_text(hint="Path to initialisation file", tag="init_path",
                           width=MESH_PATH_INPUT_WIDTH)
        dpg.add_button(label="Browse", tag="init_browse",
                       callback=open_file_dialog_callback, user_data="init")

    dpg.add_file_dialog(
        directory_selector=False, tag="file_dialog_init",
        user_data="init_path", callback=select_mesh_file_callback,
        show=False, width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT
    )
    if ".c" in FILE_EXTENSIONS:
        dpg.add_file_extension(".c", parent="file_dialog_init",
                               color=FILE_EXTENSIONS[".c"])


def _create_restart_file_section():
    dpg.add_checkbox(
        label="Restart from previous run",
        tag="restart_toggle",
        callback=show_restart_callback
    )
    with dpg.tooltip("restart_toggle"):
        dpg.add_text("Provide a .csv restart file from a previous run",
                     color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("Must use the same mesh file as the original run.")

    with dpg.group(horizontal=True, show=False, tag="restart_group"):
        dpg.add_input_text(hint="Path to restart file", tag="restart_path",
                           width=MESH_PATH_INPUT_WIDTH)
        dpg.add_button(label="Browse", tag="restart_browse",
                       callback=open_file_dialog_callback, user_data="restart")

    dpg.add_file_dialog(
        directory_selector=False, tag="file_dialog_restart",
        user_data="restart_path", callback=select_mesh_file_callback,
        show=False, width=FILE_DIALOG_WIDTH, height=FILE_DIALOG_HEIGHT
    )
    if ".csv" in FILE_EXTENSIONS:
        dpg.add_file_extension(".csv", parent="file_dialog_restart",
                               color=FILE_EXTENSIONS[".csv"])


def _create_run_button(themes: dict):
    dpg.add_checkbox(label="Run on GPU (if supported)", tag="gpu_toggle")
    dpg.add_spacer(height=6)
    run_btn = dpg.add_button(
        label="Compile and Run Solver",
        width=270,
        height=RUN_BUTTON_HEIGHT,
        callback=run_solver_callback,
        tag="run_solver_button"
    )
    with dpg.tooltip(run_btn):
        dpg.add_text("Compile and execute the solver", color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text("1. Writes configuration files")
        dpg.add_text("2. Compiles C sources with gcc")
        dpg.add_text("3. Runs the solver executable")

    if "button" in themes:
        dpg.bind_item_theme(run_btn, themes["button"])


# ── Callbacks (unchanged) ──────────────────────────────────────────────────────

def show_restart_callback(sender, app_data, user_data):
    dpg.configure_item("restart_group", show=app_data)

def show_init_callback(sender, app_data, user_data):
    dpg.configure_item("init_group", show=app_data)