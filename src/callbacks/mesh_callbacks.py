"""
Mesh and multigrid-related callbacks for MeMPhyS GUI

Handles callbacks for:
- Multigrid enable/disable
- Mesh level selection
- File browsing
- File validation
"""

import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils import validate_mesh_file
from src.config import MAX_MESH_LEVELS, FILE_DIALOG_WIDTH, FILE_DIALOG_HEIGHT


def show_multigrid_callback(sender, app_data, user_data):
    """
    Show/hide multigrid-related UI elements
    
    Args:
        sender: Checkbox tag
        app_data: Checkbox state
        user_data: User data (unused)
    """
    multigrid_on = dpg.get_value("multigrid_toggle")
    
    # Update app state
    app_state.multigrid_enabled = multigrid_on
    
    # Show/hide entire multigrid parameter group
    if dpg.does_item_exist("multigrid_parameters_section"):
        dpg.configure_item("multigrid_parameters_section", show=multigrid_on)
        dpg.configure_item("multigrid_hint_text", show=multigrid_on)  # show hint when multigrid ON
    # Get number of mesh levels
    if dpg.does_item_exist("num_mesh_levels"):
        n = int(dpg.get_value("num_mesh_levels"))
    else:
        n = 1
    
    # Show/hide mesh input groups
    if multigrid_on:
        # Show first n mesh groups
        for i in range(1, MAX_MESH_LEVELS + 1):
            show_i = (i <= n)
            _configure_mesh_group_visibility(i, show_i)
        
        logger.info(f"Multigrid enabled with {n} levels")
    else:
        # Hide all except the first
        for i in range(1, MAX_MESH_LEVELS + 1):
            show_i = (i == 1)
            _configure_mesh_group_visibility(i, show_i)
        
        logger.info("Multigrid disabled")
    
    # Update mesh inputs visibility
    if multigrid_on:
        update_mesh_inputs_callback(sender, app_data, user_data)


def update_mesh_inputs_callback(sender, app_data, user_data):
    """
    Update visibility of mesh input groups based on number of levels
    
    Args:
        sender: Input widget tag
        app_data: New value
        user_data: User data (unused)
    """
    num_levels = int(dpg.get_value("num_mesh_levels"))
    multigrid_on = dpg.get_value("multigrid_toggle")
    
    # Update app state
    app_state.num_mesh_levels = num_levels
    
    # Show/hide mesh groups
    for i in range(1, MAX_MESH_LEVELS + 1):
        show_i = multigrid_on and (i <= num_levels)
        _configure_mesh_group_visibility(i, show_i)
    
    logger.debug(f"Mesh levels set to: {num_levels}")


def _configure_mesh_group_visibility(level: int, show: bool):
    """
    Helper function to show/hide mesh group elements
    
    Args:
        level: Mesh level (1-indexed)
        show: Whether to show the elements
    """
    tags = [
        f"mesh_group_{level}",
        f"mesh_file_{level}",
        f"browse_{level}",
    ]
    
    for tag in tags:
        if dpg.does_item_exist(tag):
            dpg.configure_item(tag, show=show)
    
    # Keep file dialogs hidden (only shown when Browse is clicked)
    dialog_tag = f"file_dialog_{level}"
    if dpg.does_item_exist(dialog_tag):
        dpg.configure_item(dialog_tag, show=False)


def open_file_dialog_callback(sender, app_data, user_data):
    """
    Open file dialog for mesh file selection
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: Mesh level number (int) or file type (str)
    """
    if isinstance(user_data, int):
        # Mesh file dialog
        dialog_tag = f"file_dialog_{user_data}"
    else:
        # Other file dialogs (init, vtk, config, etc.)
        dialog_tag = f"file_dialog_{user_data}"
    
    if dpg.does_item_exist(dialog_tag):
        dpg.configure_item(dialog_tag, show=True)
    else:
        logger.error(f"File dialog not found: {dialog_tag}")


def select_mesh_file_callback(sender, app_data, user_data):
    """
    Callback when a mesh file is selected from file dialog
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: Target input text tag to update
    """
    if 'file_path_name' not in app_data:
        logger.error("No file path in file dialog data")
        return
    
    file_path = app_data['file_path_name']
    
    # Update the input field
    if dpg.does_item_exist(user_data):
        dpg.set_value(user_data, file_path)
        logger.info(f"Selected file: {file_path}")
        
        # If it's a mesh file, validate it
        if "mesh_file_" in user_data:
            validate_selected_mesh_file(file_path, user_data)
    else:
        logger.error(f"Target input field not found: {user_data}")


def validate_selected_mesh_file(file_path: str, input_tag: str):
    """
    Validate a selected mesh file and provide feedback
    
    Args:
        file_path: Path to mesh file
        input_tag: Tag of the input field
    """
    is_valid, error_msg = validate_mesh_file(file_path)
    
    if is_valid:
        logger.success(f"Mesh file validated: {file_path}")
    else:
        logger.error(f"Mesh file validation failed: {error_msg}")


def validate_all_mesh_files_callback(sender, app_data, user_data):
    """
    Validate all selected mesh files
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Validating all mesh files...")
    
    multigrid_on = dpg.get_value("multigrid_toggle")
    num_levels = 1
    
    if multigrid_on and dpg.does_item_exist("num_mesh_levels"):
        num_levels = int(dpg.get_value("num_mesh_levels"))
    
    all_valid = True
    validated_count = 0
    
    for i in range(1, num_levels + 1):
        tag = f"mesh_file_{i}"
        
        if not dpg.does_item_exist(tag):
            continue
        
        file_path = dpg.get_value(tag)
        
        if not file_path:
            logger.warning(f"Mesh file {i} is not specified")
            all_valid = False
            continue
        
        is_valid, error_msg = validate_mesh_file(file_path)
        
        if is_valid:
            validated_count += 1
            logger.success(f"Mesh {i}: Valid")
        else:
            logger.error(f"Mesh {i}: {error_msg}")
            all_valid = False
    
    if all_valid and validated_count > 0:
        logger.success(f"All {validated_count} mesh files are valid")
    elif validated_count == 0:
        logger.warning("No mesh files to validate")
    else:
        logger.warning(f"Some mesh files are invalid ({validated_count}/{num_levels} valid)")


def clear_mesh_files_callback(sender, app_data, user_data):
    """
    Clear all mesh file paths
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Clearing all mesh file paths...")
    
    cleared_count = 0
    
    for i in range(1, MAX_MESH_LEVELS + 1):
        tag = f"mesh_file_{i}"
        
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, "")
            cleared_count += 1
    
    logger.success(f"Cleared {cleared_count} mesh file paths")


def browse_init_file_callback(sender, app_data, user_data):
    """
    Open file dialog for initialization file
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if dpg.does_item_exist("file_dialog_init"):
        dpg.configure_item("file_dialog_init", show=True)


def select_init_file_callback(sender, app_data, user_data):
    """
    Callback when initialization file is selected
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: Target input field tag
    """
    if 'file_path_name' not in app_data:
        return
    
    file_path = app_data['file_path_name']
    
    if dpg.does_item_exist("init_path"):
        dpg.set_value("init_path", file_path)
        app_state.init_file_path = file_path
        logger.info(f"Initialization file selected: {file_path}")


def copy_mesh_path_callback(sender, app_data, user_data):
    """
    Copy mesh file path from one level to another
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: Tuple of (source_level, target_level)
    """
    source_level, target_level = user_data
    
    source_tag = f"mesh_file_{source_level}"
    target_tag = f"mesh_file_{target_level}"
    
    if not dpg.does_item_exist(source_tag) or not dpg.does_item_exist(target_tag):
        logger.error("Source or target mesh field not found")
        return
    
    source_path = dpg.get_value(source_tag)
    
    if not source_path:
        logger.warning(f"Mesh file {source_level} is empty, nothing to copy")
        return
    
    dpg.set_value(target_tag, source_path)
    logger.info(f"Copied mesh path from level {source_level} to {target_level}")


def auto_fill_mesh_levels_callback(sender, app_data, user_data):
    """
    Auto-fill mesh levels by assuming sequential naming pattern
    (e.g., mesh_1.msh, mesh_2.msh, mesh_3.msh)
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # Get first mesh file
    if not dpg.does_item_exist("mesh_file_1"):
        return
    
    first_mesh = dpg.get_value("mesh_file_1")
    
    if not first_mesh:
        logger.warning("First mesh file must be specified to auto-fill")
        return
    
    # Try to extract pattern
    import re
    from pathlib import Path
    
    path = Path(first_mesh)
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    
    # Look for number pattern
    match = re.search(r'(\d+)', stem)
    
    if not match:
        logger.warning("Could not detect numbering pattern in first mesh file")
        return
    
    base_name = stem[:match.start()]
    start_num = int(match.group(1))
    end_name = stem[match.end():]
    
    num_levels = int(dpg.get_value("num_mesh_levels"))
    filled_count = 0
    
    for i in range(2, num_levels + 1):
        new_num = start_num + (i - 1)
        new_name = f"{base_name}{new_num}{end_name}{suffix}"
        new_path = parent / new_name
        
        tag = f"mesh_file_{i}"
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, str(new_path))
            filled_count += 1
    
    logger.success(f"Auto-filled {filled_count} mesh file paths")
    logger.info("Please verify the paths are correct")