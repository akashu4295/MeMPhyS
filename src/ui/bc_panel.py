"""
Boundary Conditions panel for MeMPhyS GUI

Allows users to assign boundary conditions to physical entities
from the mesh file.
"""

import dearpygui.dearpygui as dpg

from src.config import BC_TYPES, BC_VARIABLES, COLORS
from src.utils.gmsh_bc_manager import (
    get_physical_names,
    get_boundary_condition,
    set_boundary_condition,
)


def create_bc_panel(themes: dict) -> int:
    """
    Create boundary conditions configuration panel
    
    Args:
        themes: Dictionary of theme IDs
        
    Returns:
        Panel tag/ID
    """
    with dpg.child_window(
        width=-1,
        height=300,
        border=True,
        tag="bc_panel"
    ) as panel:
        
        dpg.add_text("Boundary Conditions", color=COLORS["success"])
        dpg.add_separator()
        
        # Instructions
        dpg.add_text(
            "Load a mesh file to configure boundary conditions",
            tag="bc_instructions",
            color=COLORS["info"]
        )
        
        # Container for BC widgets (will be populated when mesh is loaded)
        with dpg.group(tag="bc_widgets_container"):
            pass
        
        dpg.add_spacer(height=10)
        
        # Buttons
        with dpg.group(horizontal=True):
            refresh_btn = dpg.add_button(
                label="Refresh from Mesh",
                callback=lambda: refresh_bc_panel(),
                tag="refresh_bc_button"
            )
            
            write_btn = dpg.add_button(
                label="Write bc.csv",
                callback=lambda: write_bc_csv_callback(),
                tag="write_bc_button"
            )
            
            load_btn = dpg.add_button(
                label="Load bc.csv",
                callback=lambda: load_bc_csv_callback(),
                tag="load_bc_button"
            )
            
            # Apply themes
            if "button_secondary" in themes:
                dpg.bind_item_theme(refresh_btn, themes["button_secondary"])
                dpg.bind_item_theme(write_btn, themes["button_secondary"])
                dpg.bind_item_theme(load_btn, themes["button_secondary"])
    
    return panel


def refresh_bc_panel():
    """
    Refresh the BC panel with physical names from loaded mesh
    """
    from src.core import logger
    from src.utils.gmsh_bc_manager import read_physical_names_from_msh
    import dearpygui.dearpygui as dpg
    
    # Get current mesh file
    mesh_file = None
    if dpg.does_item_exist("mesh_file_1"):
        mesh_file = dpg.get_value("mesh_file_1")
    
    if not mesh_file:
        logger.warning("No mesh file selected. Please select a mesh file first.")
        return
    
    # Read physical names
    physical_names = read_physical_names_from_msh(mesh_file)
    
    if not physical_names:
        logger.warning("No physical names found in mesh file")
        return
    
    # Clear existing widgets
    if dpg.does_item_exist("bc_widgets_container"):
        dpg.delete_item("bc_widgets_container", children_only=True)
    
    # Update instructions
    if dpg.does_item_exist("bc_instructions"):
        dpg.set_value("bc_instructions", f"Configure boundary conditions for {len(physical_names)} physical entities:")
    
    # Create BC widgets for each physical name
    with dpg.group(parent="bc_widgets_container"):
        for physical_name in physical_names:
            create_bc_widget(physical_name)


def create_bc_widget(physical_name: str):
    """
    Create widgets for configuring BC for one physical entity
    
    Args:
        physical_name: Name of the physical entity
    """
    # Get existing BC if any
    bc_data = get_boundary_condition(physical_name)
    default_type = bc_data['type'] if bc_data else BC_TYPES[0]
    default_vars = bc_data['variables'] if bc_data else {}
    
    with dpg.group(horizontal=True):
        # Physical name label
        dpg.add_text(f"{physical_name}:", color=COLORS["subheader"])
        
        # BC type combo
        combo_tag = f"bc_type_{physical_name}"
        dpg.add_combo(
            items=BC_TYPES,
            default_value=default_type,
            tag=combo_tag,
            width=150,
            callback=lambda s, a, u: on_bc_type_changed(physical_name, a)
        )
        
        # Variable inputs container
        vars_container_tag = f"bc_vars_{physical_name}"
        with dpg.group(horizontal=True, tag=vars_container_tag):
            create_variable_inputs(physical_name, default_type, default_vars)


def create_variable_inputs(physical_name: str, bc_type: str, current_values: dict):
    """
    Create input fields for variables based on BC type
    
    Args:
        physical_name: Name of physical entity
        bc_type: Type of boundary condition
        current_values: Current variable values
    """
    variables = BC_VARIABLES.get(bc_type, [])
    
    for var_name in variables:
        default_val = current_values.get(var_name, 0.0)
        
        dpg.add_text(f"{var_name}=")
        dpg.add_input_float(
            default_value=default_val,
            width=80,
            tag=f"bc_var_{physical_name}_{var_name}",
            format="%.2f",
            on_enter=True,
            callback=lambda s, a, u: update_bc_from_widgets(physical_name)
        )


def on_bc_type_changed(physical_name: str, new_type: str):
    """
    Callback when BC type is changed
    
    Args:
        physical_name: Name of physical entity
        new_type: New BC type
    """
    # Clear and recreate variable inputs
    vars_container_tag = f"bc_vars_{physical_name}"
    
    if dpg.does_item_exist(vars_container_tag):
        dpg.delete_item(vars_container_tag, children_only=True)
        
        with dpg.group(horizontal=True, parent=vars_container_tag):
            create_variable_inputs(physical_name, new_type, {})
    
    # Update BC in memory
    update_bc_from_widgets(physical_name)


def update_bc_from_widgets(physical_name: str):
    """
    Update boundary condition from widget values
    
    Args:
        physical_name: Name of physical entity
    """
    # Get BC type
    combo_tag = f"bc_type_{physical_name}"
    if not dpg.does_item_exist(combo_tag):
        return
    
    bc_type = dpg.get_value(combo_tag)
    
    # Get variable values
    variables = {}
    for var_name in BC_VARIABLES.get(bc_type, []):
        var_tag = f"bc_var_{physical_name}_{var_name}"
        if dpg.does_item_exist(var_tag):
            variables[var_name] = dpg.get_value(var_tag)
    
    # Set BC
    set_boundary_condition(physical_name, bc_type, variables)


def write_bc_csv_callback():
    """Write boundary conditions to bc.csv"""
    from src.core import logger
    from src.utils.gmsh_bc_manager import write_bc_csv, validate_bc_assignment
    
    # Update all BCs from widgets first
    physical_names = get_physical_names()
    for name in physical_names:
        update_bc_from_widgets(name)
    
    # Validate
    all_assigned, missing = validate_bc_assignment()
    
    if not all_assigned:
        logger.warning(f"Not all physical names have BCs assigned. Missing: {', '.join(missing)}")
        logger.info("Writing bc.csv anyway with assigned BCs only")
    
    # Write
    success = write_bc_csv()
    
    if success:
        logger.success("Boundary conditions saved to bc.csv")


def load_bc_csv_callback():
    """Load boundary conditions from bc.csv"""
    from src.core import logger
    from src.utils.gmsh_bc_manager import read_bc_csv
    
    success = read_bc_csv()
    
    if success:
        logger.success("Boundary conditions loaded from bc.csv")
        # Refresh panel to show loaded BCs
        refresh_bc_panel()
    else:
        logger.warning("Could not load bc.csv")