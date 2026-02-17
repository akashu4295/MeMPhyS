"""
Menu bar callbacks for MeMPhyS GUI

Handles callbacks for:
- File menu (open, save, exit)
- Edit menu (preferences)
- Help menu (help, about)
"""

import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils import open_folder, open_url, change_font, get_platform
from src.config import (
    LOG_DIR,
    HELP_URL,
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    FONT_SIZES,
    FONT_PREFERENCES,
    COLORS,
    BC_TYPES,
    BC_VARIABLES
)

from src.utils.gmsh_bc_manager import (
    get_physical_names,
    get_boundary_condition,
    set_boundary_condition,
)


def show_bc_window_callback(sender, app_data, user_data):
    """
    Show the Boundary Conditions dialog
    
    Args:
        sender: button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # If window already exists, just show it
    if dpg.does_item_exist("bc_window"):
        dpg.configure_item("bc_window", show=True)
        dpg.focus_item("bc_window")
        return
    
    # Get current font settings
    current_font_name = app_state.current_font_name
    current_font_size = app_state.current_font_size
    
    # Get available fonts for current platform
    platform = get_platform()
    available_fonts = FONT_PREFERENCES.get(platform, [])
    
    # Create Preferences window
    with dpg.window(
        label="Boundary Conditions",
        tag="bc_window",
        modal=True,
        width=650,
        height=400,
        no_resize=True,
        pos=(415, 275)
    ):
        dpg.add_text("Boundary Conditions", color=(200, 220, 255))
        dpg.add_separator()
        

        # Instructions
        dpg.add_text(
            "Load a mesh file to configure boundary conditions",
            tag="bc_dialog_instructions",
            color=COLORS["info"]
        )
        
        # Container for BC widgets (will be populated when mesh is loaded)
        with dpg.group(tag="bc_widgets_dialog_container"):
            pass
        
        dpg.add_spacer(height=10)
        
        # Buttons
        with dpg.group(horizontal=True):
            refresh_btn = dpg.add_button(
                label="Read from Mesh",
                callback=lambda: refresh_bc_callback(),
                tag="refresh_bc_dialog_button"
            )
            
            load_btn = dpg.add_button(
                label="Load bc.csv",
                callback=lambda: load_bc_csv_dialog_callback(),
                tag="load_bc_dialog_button"
            )
            
            # Apply themes
            # if "button_secondary" in themes:
            #     dpg.bind_item_theme(refresh_bc_btn, themes["button_secondary"])
            #     dpg.bind_item_theme(write_bc_btn, themes["button_secondary"])
            #     dpg.bind_item_theme(load_bc_btn, themes["button_secondary"])
    

        # Buttons
        with dpg.group(horizontal=True):
            dpg.add_button(
                label="Apply",
                width=80,
                callback=apply_bc_callback
            )
            
            dpg.add_button(
                label="Close",
                width=80,
                callback=lambda: dpg.configure_item("bc_window", show=False)
            )

def apply_bc_callback(sender, app_data, user_data):
    """
    Apply boundary condition changes
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """    
    from src.core import logger
    from src.utils.gmsh_bc_manager import write_bc_csv, validate_bc_assignment
    
    # Update all BCs from widgets first
    physical_names = get_physical_names()
    for name in physical_names:
        update_bc_from_widgets_dialog(name)
    
    # Validate
    all_assigned, missing = validate_bc_assignment()
    
    if not all_assigned:
        logger.warning(f"Not all physical names have BCs assigned. Missing: {', '.join(missing)}")
        logger.info("Writing bc.csv anyway with assigned BCs only")
    
    # Write
    success = write_bc_csv()
    
    if success:
        logger.success("Boundary conditions saved to bc.csv")

def refresh_bc_callback():
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
    if dpg.does_item_exist("bc_widgets_dialog_container"):
        dpg.delete_item("bc_widgets_dialog_container", children_only=True)
    
    # Update instructions
    if dpg.does_item_exist("bc_dialog_instructions"):
        dpg.set_value("bc_dialog_instructions", f"Configure boundary conditions for {len(physical_names)} physical entities:")
    
    # Create BC widgets for each physical name
    with dpg.group(parent="bc_widgets_dialog_container"):
        for physical_name in physical_names:
            create_bc_widget(physical_name)

def load_bc_csv_dialog_callback():
    """Load boundary conditions from bc.csv"""
    from src.core import logger
    from src.utils.gmsh_bc_manager import read_bc_csv
    
    success = read_bc_csv()
    
    if success:
        logger.success("Boundary conditions loaded from bc.csv")
        # Refresh panel to show loaded BCs
        refresh_bc_callback()
    else:
        logger.warning("Could not load bc.csv")

def update_bc_from_widgets_dialog(physical_name: str):
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
            callback=lambda s, a, u: on_bc_type_changed_dialog(physical_name, a)
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
            callback=lambda s, a, u: update_bc_from_widgets_dialog(physical_name)
        )

def on_bc_type_changed_dialog(physical_name: str, new_type: str):
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
    update_bc_from_widgets_dialog(physical_name)