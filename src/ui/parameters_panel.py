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
)


def create_parameters_panel(themes: dict) -> int:
    """
    Create the left panel with parameters and controls
    
    Args:
        themes: Dictionary of theme IDs
    
    Returns:
        Child window tag/ID
    """
    with dpg.child_window(
        width=LEFT_PANEL_WIDTH,
        border=True,
        tag="parameters_panel"
    ) as panel:
                

        # Header
        dpg.add_text("Input & Parameters", color=COLORS["header"])
        dpg.add_spacer(height=5)
        
        # Solver Method Selection
        _create_solver_method_section()
        
        dpg.add_separator()
        
        # Flow Parameters
        _create_flow_parameters_section()
        
        # Implicit Parameters (hidden by default)
        _create_implicit_parameters_section()
        
        dpg.add_separator()
        
        # Multigrid Section
        _create_multigrid_section()
        
        # Header
        dpg.add_text("Create Geometry or Browse", color=COLORS["header"])
        
        # Solver Method Selection
        _create_geometry_section(themes)

        # Mesh File Selection
        _create_mesh_files_section()
        
        dpg.add_separator()
        
        # Initialization File Selection
        _create_init_file_section()
        
        dpg.add_spacer(height=10)
        
        # Run Button
        _create_run_button(themes)
    
    return panel


def _create_solver_method_section():
    """Create solver method selection combo box"""
    dpg.add_text("Solver Method")
    dpg.add_combo(
        items=SOLVER_METHODS,
        default_value=DEFAULT_SOLVER_METHOD,
        tag="solver_method",
        width=200,
        callback=show_implicit_callback
    )


def _create_flow_parameters_section():
    """Create flow parameters input section"""
    dpg.add_text("Flow Parameters", color=COLORS["subheader"])
    
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
    """Create implicit solver parameters (hidden by default)"""
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
    """Create multigrid enable checkbox and parameters"""
    # Multigrid toggle
    dpg.add_checkbox(
        label="Enable Multigrid",
        tag="multigrid_toggle",
        callback=show_multigrid_callback
    )
    
    # Multigrid parameters group (hidden by default)
    with dpg.group(
        horizontal=False,
        tag="multigrid_parameters_section",
        show=False
    ):
        dpg.add_text("Multigrid Parameters", color=COLORS["success"])
        
        # Multigrid parameter inputs
        for pname, pval in MULTIGRID_PARAMETERS.items():
            tag = f"param_{pname}"
            dpg.add_input_text(
                label=pname,
                tag=tag,
                default_value=str(pval),
                width=PARAMETER_INPUT_WIDTH,
                callback=validate_numeric_input,
                on_enter=True
            )
        
        # Number of mesh levels
        dpg.add_input_int(
            label="Number of mesh levels",
            default_value=1,
            min_value=1,
            max_value=MAX_MESH_LEVELS,
            width=PARAMETER_INPUT_WIDTH,
            tag="num_mesh_levels",
            callback=update_mesh_inputs_callback
        )


def _create_geometry_section(themes: dict):
    """Create geometry/Gmsh section"""
    dpg.add_text("Geometry (Gmsh)", color=COLORS["success"])
    
    with dpg.group(horizontal=True):
        new_geo_btn = dpg.add_button(
            label="Open Gmsh",
            callback=launch_gmsh_callback,
            width=120
        )
        
        open_geo_btn = dpg.add_button(
            label="Open Geometry",
            callback=browse_geometry_file_callback,
            width=120
        )

        # Apply themes
        if "button_secondary" in themes:
            dpg.bind_item_theme(new_geo_btn, themes["button_secondary"])
            dpg.bind_item_theme(open_geo_btn, themes["button_secondary"])
        
    # set_mesh_btn = dpg.add_button(
    #     label="Set Mesh from Geo",
    #     callback=set_mesh_from_geometry_callback,
    #     width=140
    # )
        
    # # Apply themes
    # if "button_secondary" in themes:
    #     dpg.bind_item_theme(set_mesh_btn, themes["button_secondary"])
    
    # File dialog for .geo files
    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_geometry",
        callback=select_geometry_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT
    )
    dpg.add_file_extension(".geo", parent="file_dialog_geometry", color=(150, 255, 150, 255))
    dpg.add_file_extension(".geo_unrolled", parent="file_dialog_geometry", color=(150, 255, 150, 255))
    
    dpg.add_spacer(height=5)


def _create_mesh_files_section():
    """Create mesh file selection section"""
    with dpg.group(tag="multigrid_mesh_section"):
        dpg.add_text("Choose Mesh Files", color=COLORS["success"])
        dpg.add_text("(finest to coarsest if multigrid)", color=COLORS["success"])
        
        # Create mesh file inputs for all levels
        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            show_initial = (i == 0)  # Only show first one initially
            
            # Group for each mesh level
            with dpg.group(
                horizontal=True,
                show=show_initial,
                tag=f"mesh_group_{idx}"
            ):
                dpg.add_input_text(
                    tag=f"mesh_file_{idx}",
                    hint=f"Mesh file {idx} (choose or type path)",
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
        
        # Create file dialogs for mesh files
        for i in range(MAX_MESH_LEVELS):
            idx = i + 1
            dialog_tag = f"file_dialog_{idx}"
            input_tag = f"mesh_file_{idx}"
            
            dpg.add_file_dialog(
                directory_selector=False,
                tag=dialog_tag,
                user_data=input_tag,
                callback=select_mesh_file_callback,
                show=False,
                width=FILE_DIALOG_WIDTH,
                height=FILE_DIALOG_HEIGHT
            )
            
            # Add file extension filter
            if ".msh" in FILE_EXTENSIONS:
                dpg.add_file_extension(
                    ".msh",
                    parent=dialog_tag,
                    color=FILE_EXTENSIONS[".msh"]
                )


def _create_init_file_section():
    """Create initialization file selection section"""
    dpg.add_text("Choose the initialisation file", color=COLORS["subheader"])
    
    with dpg.group(horizontal=True, show=True, tag="init_group"):
        dpg.add_input_text(
            hint="Path to Initialisation file",
            tag="init_path",
            width=MESH_PATH_INPUT_WIDTH,
            show=True
        )
        
        dpg.add_button(
            label="Browse",
            tag="init_browse",
            callback=open_file_dialog_callback,
            user_data="init",
            show=True
        )
    
    # File dialog for init file
    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_init",
        user_data="init_path",
        callback=select_mesh_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT
    )
    
    # Add .c file extension filter
    if ".c" in FILE_EXTENSIONS:
        dpg.add_file_extension(
            ".c",
            parent="file_dialog_init",
            color=FILE_EXTENSIONS[".c"]
        )


def _create_run_button(themes: dict):
    """Create the main run button"""
    compile_btn = dpg.add_button(
        label="Compile and Run Solver",
        width=RUN_BUTTON_WIDTH,
        height=RUN_BUTTON_HEIGHT,
        callback=run_solver_callback,
        tag="run_solver_button"
    )
    
    # Add tooltip
    with dpg.tooltip(compile_btn):
        dpg.add_text(
            "Click to compile and execute the solver",
            color=COLORS["info"]
        )
        dpg.add_separator()
        dpg.add_text("This will:")
        dpg.add_text("  1. Write configuration files")
        dpg.add_text("  2. Compile all C source files")
        dpg.add_text("  3. Run the solver executable")
    
    # Apply button theme
    if "button" in themes:
        dpg.bind_item_theme(compile_btn, themes["button"])