"""
File I/O utilities

Handles reading and writing configuration files, parameter files,
and grid files with proper validation and error handling.
"""

import os
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import dearpygui.dearpygui as dpg

from src.config import (
    BASE_PARAMETERS,
    IMPLICIT_PARAMETERS,
    MULTIGRID_PARAMETERS,
    FIXED_PARAMETERS,
    PARAMS_CSV,
    GRID_CSV,
)
from src.core import logger, app_state


def validate_file_path(filepath: str, must_exist: bool = True) -> bool:
    """
    Validate that a file path is valid
    
    Args:
        filepath: Path to validate
        must_exist: Whether the file must already exist
    
    Returns:
        True if valid, False otherwise
    """
    if not filepath or not isinstance(filepath, str):
        return False
    
    try:
        path = Path(filepath)
        
        if must_exist:
            return path.exists() and path.is_file()
        else:
            # Check if parent directory exists and is writable
            parent = path.parent
            return parent.exists() and os.access(str(parent), os.W_OK)
    
    except Exception as e:
        logger.debug(f"Path validation error: {e}")
        return False


def validate_mesh_file(filepath: str) -> Tuple[bool, str]:
    """
    Validate a mesh file
    
    Args:
        filepath: Path to the mesh file
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not filepath:
        return False, "No file path provided"
    
    if not validate_file_path(filepath, must_exist=True):
        return False, f"File not found: {filepath}"
    
    # Check file extension
    if not filepath.lower().endswith('.msh'):
        return False, f"Invalid file extension. Expected .msh, got {Path(filepath).suffix}"
    
    # Check if file is readable
    try:
        with open(filepath, 'r') as f:
            # Try to read first line to verify it's readable
            f.readline()
        return True, ""
    except Exception as e:
        return False, f"Cannot read file: {e}"


def read_parameters_from_gui() -> Dict[str, Any]:
    """
    Read all parameters from the GUI input fields
    
    Returns:
        Dictionary with all parameters
    """
    params = {}
    
    # Read base parameters
    for pname in BASE_PARAMETERS.keys():
        tag = f"param_{pname}"
        if dpg.does_item_exist(tag):
            value = dpg.get_value(tag)
            params[pname] = value
    
    # Read implicit parameters
    for pname in IMPLICIT_PARAMETERS.keys():
        tag = f"param_{pname}"
        if dpg.does_item_exist(tag):
            value = dpg.get_value(tag)
            params[pname] = value
    
    # Read multigrid parameters
    for pname in MULTIGRID_PARAMETERS.keys():
        tag = f"param_{pname}"
        if dpg.does_item_exist(tag):
            value = dpg.get_value(tag)
            params[pname] = value
    
    # Read solver method
    if dpg.does_item_exist("solver_method"):
        solver_method = dpg.get_value("solver_method")
        params["fractional_step"] = 1 if solver_method == "Fractional Step" else 0
    
    # Add fixed parameters
    params.update(FIXED_PARAMETERS)
    
    return params


def write_parameters_csv(filename: str = PARAMS_CSV) -> bool:
    """
    Write parameters to CSV file
    
    Args:
        filename: Output filename
    
    Returns:
        True if successful, False otherwise
    """
    try:
        logger.info(f"Writing parameters to {filename}")
        
        # Get parameters from GUI
        params = read_parameters_from_gui()
        
        # Create list of rows
        rows = [[param_name, value] for param_name, value in params.items()]
        
        # Write to CSV
        df = pd.DataFrame(rows, columns=["Parameter", "Value"])
        df.to_csv(filename, index=False, header=False)
        
        logger.log_file_operation("Wrote", filename, success=True)
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error writing {filename}")
        return False


def read_mesh_files_from_gui() -> List[str]:
    """
    Read mesh file paths from GUI
    
    Returns:
        List of mesh file paths
    """
    mesh_files = []
    num_levels = 1
    
    # Get number of levels if multigrid is enabled
    if dpg.does_item_exist("multigrid_toggle") and dpg.get_value("multigrid_toggle"):
        if dpg.does_item_exist("num_mesh_levels"):
            num_levels = int(dpg.get_value("num_mesh_levels"))
    
    # Read mesh file paths
    for i in range(1, num_levels + 1):
        tag = f"mesh_file_{i}"
        if dpg.does_item_exist(tag):
            path = dpg.get_value(tag)
            if path:
                mesh_files.append(path)
            else:
                logger.warning(f"Mesh file {i} is empty")
    
    return mesh_files


def write_grid_csv(filename: str = GRID_CSV) -> bool:
    """
    Write grid configuration to CSV file
    
    Args:
        filename: Output filename
    
    Returns:
        True if successful, False otherwise
    """
    try:
        logger.info(f"Writing grid configuration to {filename}")
        
        # Get mesh files from GUI
        mesh_files = read_mesh_files_from_gui()
        
        if not mesh_files:
            logger.error("No mesh files specified")
            return False
        
        # Validate all mesh files
        for i, mesh_file in enumerate(mesh_files, 1):
            is_valid, error_msg = validate_mesh_file(mesh_file)
            if not is_valid:
                logger.error(f"Mesh file {i} validation failed: {error_msg}")
                return False
        
        # Write to CSV
        num_levels = len(mesh_files)
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["num_levels", num_levels])
            for mesh_file in mesh_files:
                writer.writerow([mesh_file])
        
        logger.log_file_operation("Wrote", filename, success=True)
        logger.info(f"Configured {num_levels} mesh level(s)")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error writing {filename}")
        return False


def read_csv_file(filename: str) -> Optional[pd.DataFrame]:
    """
    Read a CSV file into a pandas DataFrame
    
    Args:
        filename: Path to CSV file
    
    Returns:
        DataFrame or None if reading fails
    """
    try:
        if not os.path.exists(filename):
            logger.warning(f"File not found: {filename}")
            return None
        
        df = pd.read_csv(filename)
        logger.debug(f"Read {len(df)} rows from {filename}")
        return df
    
    except Exception as e:
        logger.log_exception(e, f"Error reading {filename}")
        return None


def write_csv_file(df: pd.DataFrame, filename: str, **kwargs) -> bool:
    """
    Write a DataFrame to CSV file
    
    Args:
        df: DataFrame to write
        filename: Output filename
        **kwargs: Additional arguments passed to to_csv()
    
    Returns:
        True if successful, False otherwise
    """
    try:
        df.to_csv(filename, **kwargs)
        logger.log_file_operation("Wrote", filename, success=True)
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error writing {filename}")
        return False


def create_backup(filepath: str) -> bool:
    """
    Create a backup of a file
    
    Args:
        filepath: Path to file to backup
    
    Returns:
        True if successful, False otherwise
    """
    try:
        if not os.path.exists(filepath):
            logger.warning(f"Cannot backup non-existent file: {filepath}")
            return False
        
        # Create backup filename
        backup_path = f"{filepath}.bak"
        
        # Copy file
        import shutil
        shutil.copy2(filepath, backup_path)
        
        logger.info(f"Created backup: {backup_path}")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error creating backup")
        return False


def ensure_directory_exists(dirpath: str) -> bool:
    """
    Ensure a directory exists, creating it if necessary
    
    Args:
        dirpath: Path to directory
    
    Returns:
        True if directory exists or was created, False otherwise
    """
    try:
        Path(dirpath).mkdir(parents=True, exist_ok=True)
        return True
    except Exception as e:
        logger.log_exception(e, f"Error creating directory {dirpath}")
        return False


def get_file_info(filepath: str) -> Dict[str, Any]:
    """
    Get information about a file
    
    Args:
        filepath: Path to file
    
    Returns:
        Dictionary with file information
    """
    try:
        path = Path(filepath)
        
        if not path.exists():
            return {"exists": False}
        
        stat = path.stat()
        
        return {
            "exists": True,
            "name": path.name,
            "size": stat.st_size,
            "size_mb": stat.st_size / (1024 * 1024),
            "modified": stat.st_mtime,
            "is_file": path.is_file(),
            "is_dir": path.is_dir(),
            "extension": path.suffix,
        }
    
    except Exception as e:
        logger.log_exception(e, f"Error getting file info for {filepath}")
        return {"exists": False, "error": str(e)}


def load_configuration(config_file: str) -> Optional[Dict[str, Any]]:
    """
    Load configuration from a file (placeholder for future implementation)
    
    Args:
        config_file: Path to configuration file
    
    Returns:
        Configuration dictionary or None if loading fails
    """
    # TODO: Implement configuration loading
    logger.warning("Configuration loading not yet implemented")
    return None


def save_configuration(config_file: str, config: Dict[str, Any]) -> bool:
    """
    Save configuration to a file (placeholder for future implementation)
    
    Args:
        config_file: Path to save configuration
        config: Configuration dictionary
    
    Returns:
        True if successful, False otherwise
    """
    # TODO: Implement configuration saving
    logger.warning("Configuration saving not yet implemented")
    return False


def list_files_in_directory(dirpath: str, extension: Optional[str] = None) -> List[str]:
    """
    List all files in a directory
    
    Args:
        dirpath: Directory path
        extension: Optional file extension filter (e.g., ".msh")
    
    Returns:
        List of file paths
    """
    try:
        path = Path(dirpath)
        
        if not path.exists() or not path.is_dir():
            logger.warning(f"Directory not found: {dirpath}")
            return []
        
        if extension:
            files = list(path.glob(f"*{extension}"))
        else:
            files = [f for f in path.iterdir() if f.is_file()]
        
        return [str(f) for f in files]
    
    except Exception as e:
        logger.log_exception(e, f"Error listing files in {dirpath}")
        return []


def read_text_file(filepath: str) -> Optional[str]:
    """
    Read a text file
    
    Args:
        filepath: Path to text file
    
    Returns:
        File contents or None if reading fails
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        return content
    
    except Exception as e:
        logger.log_exception(e, f"Error reading {filepath}")
        return None


def write_text_file(filepath: str, content: str) -> bool:
    """
    Write content to a text file
    
    Args:
        filepath: Path to output file
        content: Content to write
    
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(filepath, 'w') as f:
            f.write(content)
        logger.log_file_operation("Wrote", filepath, success=True)
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error writing {filepath}")
        return False