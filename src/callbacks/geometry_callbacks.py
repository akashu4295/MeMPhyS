"""
Geometry and Gmsh-related callbacks
"""

import os
import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils.gmsh_bc_manager import (
    launch_gmsh,
    get_last_gmsh_file,
    set_last_gmsh_file,
    get_mesh_from_gmsh_file,
    read_physical_names_from_msh,
)


def launch_gmsh_callback(sender, app_data, user_data):
    """
    Launch Gmsh for geometry creation
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Launching Gmsh geometry editor...")
    
    # Get last Gmsh file if any
    last_file = get_last_gmsh_file()
    
    if last_file and os.path.exists(last_file):
        success = launch_gmsh(last_file)
    else:
        success = launch_gmsh()
    
    if success:
        logger.info("Gmsh launched successfully")
        logger.info("After creating geometry, save as .geo file and mesh it (2D/3D)")
        logger.info("Then click 'Set Mesh from Geometry' to load the mesh")
    else:
        logger.error("Failed to launch Gmsh")
        logger.info("Make sure Gmsh is installed and in your PATH")


def browse_geometry_file_callback(sender, app_data, user_data):
    """
    Browse for .geo file
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if dpg.does_item_exist("file_dialog_geometry"):
        dpg.configure_item("file_dialog_geometry", show=True)


def select_geometry_file_callback(sender, app_data, user_data):
    """
    Callback when .geo file is selected
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: User data (unused)
    """
    if 'file_path_name' not in app_data:
        return
    
    geo_file = app_data['file_path_name']
    
    # Save as last Gmsh file
    set_last_gmsh_file(geo_file)
    
    logger.info(f"Geometry file selected: {geo_file}")
    
    # Get corresponding mesh file
    mesh_file = get_mesh_from_gmsh_file(geo_file)
    
    if os.path.exists(mesh_file):
        # Auto-populate mesh file field
        if dpg.does_item_exist("mesh_file_1"):
            dpg.set_value("mesh_file_1", mesh_file)
            logger.success(f"Mesh file set: {mesh_file}")
            
            # Read physical names and refresh BC panel
            physical_names = read_physical_names_from_msh(mesh_file)
            
            if physical_names:
                from src.ui.bc_panel import refresh_bc_panel
                refresh_bc_panel()
        else:
            logger.warning("Could not set mesh file - widget not found")
    else:
        logger.warning(f"Mesh file not found: {mesh_file}")
        logger.info("Please generate mesh in Gmsh (Mesh → 2D/3D)")


def set_mesh_from_geometry_callback(sender, app_data, user_data):
    """
    Set mesh file from last opened geometry file
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    geo_file = get_last_gmsh_file()
    
    if not geo_file:
        logger.warning("No geometry file opened yet")
        logger.info("Please use 'Open Geometry' or 'New Geometry' first")
        return
    
    mesh_file = get_mesh_from_gmsh_file(geo_file)
    
    if not os.path.exists(mesh_file):
        logger.error(f"Mesh file not found: {mesh_file}")
        logger.info("Please generate mesh in Gmsh (Mesh → 2D/3D)")
        return
    
    # Set mesh file
    if dpg.does_item_exist("mesh_file_1"):
        dpg.set_value("mesh_file_1", mesh_file)
        logger.success(f"Mesh file set: {mesh_file}")
        
        # Read physical names and refresh BC panel
        physical_names = read_physical_names_from_msh(mesh_file)
        
        if physical_names:
            from src.ui.bc_panel import refresh_bc_panel
            refresh_bc_panel()
    else:
        logger.error("Could not set mesh file - widget not found")