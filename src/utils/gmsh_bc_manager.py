"""
Gmsh integration and boundary condition management

Handles:
- Launching Gmsh for geometry creation
- Reading physical names from mesh files
- Managing boundary conditions
- Writing bc.csv file
"""

import os
import subprocess
import csv
from typing import Dict, List, Tuple, Optional
from pathlib import Path

from src.core import logger, app_state
from src.config import GMSH_EXECUTABLE, GMSH_FILE_EXTENSION, BC_CSV_FILE, BC_TYPES, BC_VARIABLES


# Store boundary conditions in memory
_boundary_conditions: Dict[str, Dict] = {}
# Store physical names from mesh
_physical_names: List[str] = []
# Store last opened Gmsh file
_last_gmsh_file: str = ""


def launch_gmsh(geometry_file: Optional[str] = None) -> bool:
    """
    Launch Gmsh GUI
    
    Args:
        geometry_file: Optional .geo file to open
        
    Returns:
        True if launched successfully, False otherwise
    """
    try:
        from src.utils import check_command_exists
        
        # Check if Gmsh is installed
        if not check_command_exists(GMSH_EXECUTABLE):
            logger.error(f"{GMSH_EXECUTABLE} not found in PATH")
            logger.info("Please install Gmsh: https://gmsh.info/")
            return False
        
        # Build command
        cmd = [GMSH_EXECUTABLE]
        
        if geometry_file and os.path.exists(geometry_file):
            cmd.append(geometry_file)
            logger.info(f"Launching Gmsh with {geometry_file}")
        else:
            logger.info("Launching Gmsh")
        
        # Launch Gmsh (non-blocking)
        subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        logger.success("Gmsh launched")
        return True
    
    except Exception as e:
        logger.log_exception(e, "Error launching Gmsh")
        return False


def read_physical_names_from_msh(mesh_file: str) -> List[str]:
    """
    Read physical entity names from a Gmsh .msh file
    
    Args:
        mesh_file: Path to .msh file
        
    Returns:
        List of physical names
    """
    global _physical_names
    
    physical_names = []
    
    try:
        if not os.path.exists(mesh_file):
            logger.error(f"Mesh file not found: {mesh_file}")
            return []
        
        with open(mesh_file, 'r') as f:
            content = f.read()
        
        # Look for $PhysicalNames section
        if '$PhysicalNames' not in content:
            logger.warning(f"No physical names found in {mesh_file}")
            return []
        
        # Parse physical names
        lines = content.split('\n')
        in_physical_names = False
        
        for line in lines:
            line = line.strip()
            
            if line == '$PhysicalNames':
                in_physical_names = True
                continue
            
            if line == '$EndPhysicalNames':
                break
            
            if in_physical_names and line and not line.isdigit():
                # Format: dimension tag "name"
                parts = line.split('"')
                if len(parts) >= 2:
                    name = parts[1]
                    physical_names.append(name)
        
        _physical_names = physical_names
        logger.success(f"Found {len(physical_names)} physical names: {', '.join(physical_names)}")
        return physical_names
    
    except Exception as e:
        logger.log_exception(e, f"Error reading physical names from {mesh_file}")
        return []


def get_physical_names() -> List[str]:
    """Get cached physical names"""
    return _physical_names.copy()


def set_boundary_condition(physical_name: str, bc_type: str, variables: Dict[str, float]):
    """
    Set boundary condition for a physical name
    
    Args:
        physical_name: Name of physical entity
        bc_type: Type of boundary condition
        variables: Dictionary of variable: value pairs
    """
    global _boundary_conditions
    
    _boundary_conditions[physical_name] = {
        "type": bc_type,
        "variables": variables.copy()
    }
    
    logger.info(f"Set BC for '{physical_name}': {bc_type} with {variables}")


def get_boundary_condition(physical_name: str) -> Optional[Dict]:
    """
    Get boundary condition for a physical name
    
    Args:
        physical_name: Name of physical entity
        
    Returns:
        Dictionary with 'type' and 'variables', or None
    """
    return _boundary_conditions.get(physical_name)


def get_all_boundary_conditions() -> Dict[str, Dict]:
    """Get all boundary conditions"""
    return _boundary_conditions.copy()


def clear_boundary_conditions():
    """Clear all boundary conditions"""
    global _boundary_conditions
    _boundary_conditions = {}
    logger.info("Cleared all boundary conditions")


def write_bc_csv(filename: str = BC_CSV_FILE) -> bool:
    """
    Write boundary conditions to CSV file
    
    Format: physical_name,bc_type,var1=val1,var2=val2,...
    
    Args:
        filename: Output CSV filename
        
    Returns:
        True if successful, False otherwise
    """
    try:
        if not _boundary_conditions:
            logger.warning("No boundary conditions to write")
            return False
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            
            for physical_name, bc_data in _boundary_conditions.items():
                bc_type = bc_data['type']
                variables = bc_data['variables']
                
                # Build row: name, type, var1=val1, var2=val2, ...
                row = [physical_name, bc_type]
                
                for var_name, var_value in variables.items():
                    row.append(f"{var_name}={var_value}")
                
                writer.writerow(row)
        
        logger.success(f"Boundary conditions written to {filename}")
        logger.info(f"Wrote {len(_boundary_conditions)} boundary conditions")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error writing {filename}")
        return False


def read_bc_csv(filename: str = BC_CSV_FILE) -> bool:
    """
    Read boundary conditions from CSV file
    
    Args:
        filename: CSV filename to read
        
    Returns:
        True if successful, False otherwise
    """
    global _boundary_conditions
    
    try:
        if not os.path.exists(filename):
            logger.info(f"No existing {filename} found")
            return False
        
        _boundary_conditions = {}
        
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            
            for row in reader:
                if len(row) < 2:
                    continue
                
                physical_name = row[0]
                bc_type = row[1]
                variables = {}
                
                # Parse variable=value pairs
                for item in row[2:]:
                    if '=' in item:
                        var_name, var_value = item.split('=', 1)
                        try:
                            variables[var_name] = float(var_value)
                        except ValueError:
                            variables[var_name] = 0.0
                
                _boundary_conditions[physical_name] = {
                    "type": bc_type,
                    "variables": variables
                }
        
        logger.success(f"Loaded {len(_boundary_conditions)} boundary conditions from {filename}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error reading {filename}")
        return False


def get_last_gmsh_file() -> str:
    """Get the last opened/saved Gmsh file"""
    return _last_gmsh_file


def set_last_gmsh_file(filepath: str):
    """Set the last opened/saved Gmsh file"""
    global _last_gmsh_file
    _last_gmsh_file = filepath


def validate_bc_assignment() -> Tuple[bool, List[str]]:
    """
    Validate that all physical names have boundary conditions assigned
    
    Returns:
        Tuple of (all_assigned, list_of_missing_names)
    """
    missing = []
    
    for name in _physical_names:
        if name not in _boundary_conditions:
            missing.append(name)
    
    all_assigned = len(missing) == 0
    
    if all_assigned:
        logger.success("All physical names have boundary conditions assigned")
    else:
        logger.warning(f"{len(missing)} physical names missing BCs: {', '.join(missing)}")
    
    return all_assigned, missing