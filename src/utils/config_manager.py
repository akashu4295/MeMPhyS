"""
Configuration management for MeMPhyS GUI

Handles saving and loading of:
- User preferences (fonts, UI settings)
- Application options
- Last session state (parameters, files)
"""

import os
import json
from pathlib import Path
from typing import Dict, Any, Optional
import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.config import (
    BASE_PARAMETERS,
    IMPLICIT_PARAMETERS,
    MULTIGRID_PARAMETERS,
    DEFAULT_OPTIONS,
    MAX_MESH_LEVELS,
)


# Configuration file paths
CONFIG_DIR = "./config"
USER_PREFERENCES_FILE = os.path.join(CONFIG_DIR, "user_preferences.json")
APP_OPTIONS_FILE = os.path.join(CONFIG_DIR, "app_options.json")
LAST_SESSION_FILE = os.path.join(CONFIG_DIR, "last_session.json")


def ensure_config_dir():
    """Ensure configuration directory exists"""
    Path(CONFIG_DIR).mkdir(parents=True, exist_ok=True)


# ==================== User Preferences ====================

def save_user_preferences() -> bool:
    """
    Save user preferences (fonts, UI settings)
    
    Returns:
        True if successful, False otherwise
    """
    try:
        ensure_config_dir()
        
        preferences = {
            "font_name": app_state.current_font_name,
            "font_size": app_state.current_font_size,
            "version": "2.2"
        }
        
        with open(USER_PREFERENCES_FILE, 'w') as f:
            json.dump(preferences, f, indent=4)
        
        logger.success(f"User preferences saved to {USER_PREFERENCES_FILE}")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error saving user preferences")
        return False


def load_user_preferences() -> Optional[Dict[str, Any]]:
    """
    Load user preferences
    
    Returns:
        Dictionary of preferences or None if failed
    """
    try:
        if not os.path.exists(USER_PREFERENCES_FILE):
            logger.info("No saved preferences found, using defaults")
            return None
        
        with open(USER_PREFERENCES_FILE, 'r') as f:
            preferences = json.load(f)
        
        logger.success("User preferences loaded")
        return preferences
    
    except Exception as e:
        logger.log_exception(e, "Error loading user preferences")
        return None


def apply_user_preferences(preferences: Dict[str, Any]) -> bool:
    """
    Apply loaded user preferences
    
    Args:
        preferences: Dictionary of preferences
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if not preferences:
            return False
        
        # Apply font settings
        font_name = preferences.get("font_name")
        font_size = preferences.get("font_size")
        
        if font_name and font_size:
            app_state.current_font_name = font_name
            app_state.current_font_size = font_size
            
            # Font will be applied during initialization
            logger.info(f"Will use saved font: {font_name} @ {font_size}px")
        
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error applying user preferences")
        return False


# ==================== Application Options ====================

def save_app_options() -> bool:
    """
    Save application options (all those checkboxes)
    
    Returns:
        True if successful, False otherwise
    """
    try:
        ensure_config_dir()
        
        options = app_state.get_all_options()
        
        with open(APP_OPTIONS_FILE, 'w') as f:
            json.dump(options, f, indent=4)
        
        logger.success(f"Application options saved to {APP_OPTIONS_FILE}")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error saving application options")
        return False


def load_app_options() -> Optional[Dict[str, bool]]:
    """
    Load application options
    
    Returns:
        Dictionary of options or None if failed
    """
    try:
        if not os.path.exists(APP_OPTIONS_FILE):
            logger.info("No saved options found, using defaults")
            return None
        
        with open(APP_OPTIONS_FILE, 'r') as f:
            options = json.load(f)
        
        logger.success("Application options loaded")
        return options
    
    except Exception as e:
        logger.log_exception(e, "Error loading application options")
        return None


def apply_app_options(options: Dict[str, bool]) -> bool:
    """
    Apply loaded application options
    
    Args:
        options: Dictionary of options
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if not options:
            return False
        
        # Apply options to app_state
        for key, value in options.items():
            app_state.set_option(key, value)
        
        logger.info("Application options applied")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error applying application options")
        return False


# ==================== Last Session State ====================

def save_session_state() -> bool:
    """
    Save current session state (parameters, files, etc.)
    
    Returns:
        True if successful, False otherwise
    """
    try:
        ensure_config_dir()
        
        # Check if GUI widgets still exist
        if not dpg.does_item_exist("solver_method"):
            logger.debug("GUI not available, skipping session save")
            return False
        
        session = {
            "version": "2.2",
            "solver_method": dpg.get_value("solver_method") if dpg.does_item_exist("solver_method") else "Fractional Step",
            "parameters": {},
            "multigrid_enabled": dpg.get_value("multigrid_toggle") if dpg.does_item_exist("multigrid_toggle") else False,
            "num_mesh_levels": dpg.get_value("num_mesh_levels") if dpg.does_item_exist("num_mesh_levels") else 1,
            "mesh_files": {},
            "init_file": dpg.get_value("init_path") if dpg.does_item_exist("init_path") else "",
            "vtk_file": dpg.get_value("contour_vtk_path") if dpg.does_item_exist("contour_vtk_path") else "",
            "plot_variable": dpg.get_value("contour_var") if dpg.does_item_exist("contour_var") else "velocity magnitude",
            "colormap": dpg.get_value("contour_cmap") if dpg.does_item_exist("contour_cmap") else "viridis",
        }
        
        # Save all parameters
        for pname in BASE_PARAMETERS.keys():
            tag = f"param_{pname}"
            if dpg.does_item_exist(tag):
                session["parameters"][pname] = dpg.get_value(tag)
        
        for pname in IMPLICIT_PARAMETERS.keys():
            tag = f"param_{pname}"
            if dpg.does_item_exist(tag):
                session["parameters"][pname] = dpg.get_value(tag)
        
        for pname in MULTIGRID_PARAMETERS.keys():
            tag = f"param_{pname}"
            if dpg.does_item_exist(tag):
                session["parameters"][pname] = dpg.get_value(tag)
        
        # Save mesh files
        for i in range(1, MAX_MESH_LEVELS + 1):
            tag = f"mesh_file_{i}"
            if dpg.does_item_exist(tag):
                path = dpg.get_value(tag)
                if path:
                    session["mesh_files"][str(i)] = path
        
        with open(LAST_SESSION_FILE, 'w') as f:
            json.dump(session, f, indent=4)
        
        logger.success(f"Session state saved to {LAST_SESSION_FILE}")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error saving session state")
        return False


def load_session_state() -> Optional[Dict[str, Any]]:
    """
    Load last session state
    
    Returns:
        Dictionary of session state or None if failed
    """
    try:
        if not os.path.exists(LAST_SESSION_FILE):
            logger.info("No saved session found")
            return None
        
        with open(LAST_SESSION_FILE, 'r') as f:
            session = json.load(f)
        
        logger.success("Last session state loaded")
        return session
    
    except Exception as e:
        logger.log_exception(e, "Error loading session state")
        return None


def restore_session_state(session: Dict[str, Any]) -> bool:
    """
    Restore session state to GUI
    
    Args:
        session: Dictionary of session state
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if not session:
            return False
        
        logger.info("Restoring previous session...")
        
        # Restore solver method
        if "solver_method" in session and dpg.does_item_exist("solver_method"):
            dpg.set_value("solver_method", session["solver_method"])
        
        # Restore parameters
        if "parameters" in session:
            for pname, value in session["parameters"].items():
                tag = f"param_{pname}"
                if dpg.does_item_exist(tag):
                    dpg.set_value(tag, value)
        
        # Restore multigrid settings
        if "multigrid_enabled" in session and dpg.does_item_exist("multigrid_toggle"):
            dpg.set_value("multigrid_toggle", session["multigrid_enabled"])
        
        if "num_mesh_levels" in session and dpg.does_item_exist("num_mesh_levels"):
            dpg.set_value("num_mesh_levels", session["num_mesh_levels"])
        
        # Restore mesh files
        if "mesh_files" in session:
            for level_str, path in session["mesh_files"].items():
                tag = f"mesh_file_{level_str}"
                if dpg.does_item_exist(tag):
                    dpg.set_value(tag, path)
        
        # Restore init file
        if "init_file" in session and dpg.does_item_exist("init_path"):
            dpg.set_value("init_path", session["init_file"])
        
        # Restore visualization settings
        if "vtk_file" in session and dpg.does_item_exist("contour_vtk_path"):
            dpg.set_value("contour_vtk_path", session["vtk_file"])
        
        if "plot_variable" in session and dpg.does_item_exist("contour_var"):
            dpg.set_value("contour_var", session["plot_variable"])
        
        if "colormap" in session and dpg.does_item_exist("contour_cmap"):
            dpg.set_value("contour_cmap", session["colormap"])
        
        logger.success("Previous session restored")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error restoring session state")
        return False


# ==================== Save/Load Configuration (Combined) ====================

def save_configuration(filename: Optional[str] = None) -> bool:
    """
    Save complete configuration to a file
    
    Args:
        filename: Optional custom filename
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if filename is None:
            ensure_config_dir()
            filename = os.path.join(CONFIG_DIR, "full_config.json")
        
        config = {
            "version": "2.2",
            "saved_at": str(Path(filename).name),
            "preferences": {
                "font_name": app_state.current_font_name,
                "font_size": app_state.current_font_size,
            },
            "options": app_state.get_all_options(),
            "session": {}
        }
        
        # Get session state
        if dpg.does_item_exist("solver_method"):
            session_data = {
                "solver_method": dpg.get_value("solver_method"),
                "parameters": {},
                "multigrid_enabled": dpg.get_value("multigrid_toggle") if dpg.does_item_exist("multigrid_toggle") else False,
                "num_mesh_levels": dpg.get_value("num_mesh_levels") if dpg.does_item_exist("num_mesh_levels") else 1,
                "mesh_files": {},
                "init_file": dpg.get_value("init_path") if dpg.does_item_exist("init_path") else "",
            }
            
            # Parameters
            for pname in list(BASE_PARAMETERS.keys()) + list(IMPLICIT_PARAMETERS.keys()) + list(MULTIGRID_PARAMETERS.keys()):
                tag = f"param_{pname}"
                if dpg.does_item_exist(tag):
                    session_data["parameters"][pname] = dpg.get_value(tag)
            
            # Mesh files
            for i in range(1, MAX_MESH_LEVELS + 1):
                tag = f"mesh_file_{i}"
                if dpg.does_item_exist(tag):
                    path = dpg.get_value(tag)
                    if path:
                        session_data["mesh_files"][str(i)] = path
            
            config["session"] = session_data
        
        with open(filename, 'w') as f:
            json.dump(config, f, indent=4)
        
        logger.success(f"Configuration saved to {filename}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error saving configuration to {filename}")
        return False


def load_configuration(filename: str) -> bool:
    """
    Load complete configuration from a file
    
    Args:
        filename: Configuration file path
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if not os.path.exists(filename):
            logger.error(f"Configuration file not found: {filename}")
            return False
        
        with open(filename, 'r') as f:
            config = json.load(f)
        
        logger.info(f"Loading configuration from {filename}")
        
        # Apply preferences
        if "preferences" in config:
            apply_user_preferences(config["preferences"])
        
        # Apply options
        if "options" in config:
            apply_app_options(config["options"])
        
        # Restore session
        if "session" in config and dpg.does_item_exist("solver_method"):
            restore_session_state(config["session"])
        
        logger.success("Configuration loaded successfully")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error loading configuration from {filename}")
        return False


# ==================== Auto-save on Exit ====================

def auto_save_on_exit():
    """
    Automatically save all configurations on exit
    """
    logger.info("Auto-saving configuration...")
    save_user_preferences()
    save_app_options()
    save_session_state()
    logger.info("Auto-save complete")