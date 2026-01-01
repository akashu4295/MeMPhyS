"""
Output file management for MeMPhyS GUI

Handles automatic organization of solver output files into dated folders.
"""

import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Optional

from src.core import logger, app_state
from src.config import OUTPUT_DIR, SOLUTION_CSV, DEFAULT_VTK, CONVERGENCE_CSV


def get_output_folder_path() -> str:
    """
    Get the output folder path for today
    
    Returns:
        Path to today's output folder
    """
    today_str = datetime.today().strftime("%Y-%m-%d")
    output_path = os.path.join(OUTPUT_DIR, today_str)
    return output_path


def ensure_output_folder_exists() -> str:
    """
    Ensure the output folder for today exists
    
    Returns:
        Path to the output folder
    """
    output_path = get_output_folder_path()
    
    try:
        Path(output_path).mkdir(parents=True, exist_ok=True)
        return output_path
    except Exception as e:
        logger.log_exception(e, f"Error creating output folder {output_path}")
        return OUTPUT_DIR


def get_unique_filename(directory: str, base_name: str) -> str:
    """
    Get a unique filename by appending numbers if file exists
    
    Args:
        directory: Target directory
        base_name: Base filename (e.g., "Solution.csv")
    
    Returns:
        Unique filename with number appended if necessary
    """
    file_path = os.path.join(directory, base_name)
    
    # If file doesn't exist, use base name
    if not os.path.exists(file_path):
        return base_name
    
    # Split name and extension
    name_parts = base_name.rsplit('.', 1)
    if len(name_parts) == 2:
        name, ext = name_parts
    else:
        name = base_name
        ext = ""
    
    # Find next available number
    counter = 1
    while True:
        if ext:
            new_name = f"{name}_{counter}.{ext}"
        else:
            new_name = f"{name}_{counter}"
        
        new_path = os.path.join(directory, new_name)
        
        if not os.path.exists(new_path):
            return new_name
        
        counter += 1
        
        # Safety limit
        if counter > 1000:
            logger.error("Too many files with same name, using timestamp")
            timestamp = datetime.now().strftime("%H%M%S")
            if ext:
                return f"{name}_{timestamp}.{ext}"
            else:
                return f"{name}_{timestamp}"


def move_file_to_output(source_file: str, destination_folder: str, 
                        use_unique_name: bool = True) -> Tuple[bool, str]:
    """
    Move a file to the output folder
    
    Args:
        source_file: Source file path
        destination_folder: Destination folder path
        use_unique_name: Whether to generate unique name if file exists
    
    Returns:
        Tuple of (success, destination_path)
    """
    if not os.path.exists(source_file):
        logger.warning(f"Source file not found: {source_file}")
        return False, ""
    
    try:
        base_name = os.path.basename(source_file)
        
        if use_unique_name:
            dest_name = get_unique_filename(destination_folder, base_name)
        else:
            dest_name = base_name
        
        dest_path = os.path.join(destination_folder, dest_name)
        
        # Move file
        shutil.move(source_file, dest_path)
        
        logger.success(f"Moved {base_name} → {destination_folder}/{dest_name}")
        return True, dest_path
    
    except Exception as e:
        logger.log_exception(e, f"Error moving {source_file}")
        return False, ""


def organize_solver_outputs() -> List[str]:
    """
    Organize all solver output files into dated output folder
    
    Returns:
        List of moved file paths
    """
    # Check if auto-move is enabled
    if not app_state.get_option("auto_move_outputs", True):
        logger.info("Auto-move outputs is disabled")
        return []
    
    logger.info("Organizing solver output files...")
    
    # Create output folder
    if app_state.get_option("create_dated_folders", True):
        output_folder = ensure_output_folder_exists()
    else:
        output_folder = OUTPUT_DIR
        Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Files to move
    output_files = [
        SOLUTION_CSV,
        DEFAULT_VTK,
        CONVERGENCE_CSV,
    ]
    
    moved_files = []
    use_unique = app_state.get_option("append_number_if_exists", True)
    
    for file_path in output_files:
        if os.path.exists(file_path):
            success, dest_path = move_file_to_output(
                file_path, 
                output_folder, 
                use_unique_name=use_unique
            )
            
            if success:
                moved_files.append(dest_path)
        else:
            logger.debug(f"Output file not found: {file_path}")
    
    if moved_files:
        logger.success(f"Organized {len(moved_files)} output files to {output_folder}")
    else:
        logger.info("No output files to organize")
    
    return moved_files


def list_output_runs(date_str: Optional[str] = None) -> List[str]:
    """
    List all output runs for a given date
    
    Args:
        date_str: Date string (YYYY-MM-DD), None for today
    
    Returns:
        List of output folders
    """
    if date_str is None:
        date_str = datetime.today().strftime("%Y-%m-%d")
    
    output_path = os.path.join(OUTPUT_DIR, date_str)
    
    if not os.path.exists(output_path):
        return []
    
    try:
        files = os.listdir(output_path)
        return sorted(files)
    except Exception as e:
        logger.log_exception(e, f"Error listing outputs for {date_str}")
        return []


def get_latest_output_folder() -> Optional[str]:
    """
    Get the path to the latest output folder
    
    Returns:
        Path to latest output folder or None
    """
    if not os.path.exists(OUTPUT_DIR):
        return None
    
    try:
        # Get all date folders
        folders = [f for f in os.listdir(OUTPUT_DIR) if os.path.isdir(os.path.join(OUTPUT_DIR, f))]
        
        if not folders:
            return None
        
        # Sort and get latest
        folders.sort(reverse=True)
        latest = os.path.join(OUTPUT_DIR, folders[0])
        
        return latest
    
    except Exception as e:
        logger.log_exception(e, "Error getting latest output folder")
        return None


def clean_old_outputs(days_to_keep: int = 30) -> int:
    """
    Clean output folders older than specified days
    
    Args:
        days_to_keep: Number of days to keep
    
    Returns:
        Number of folders deleted
    """
    if not os.path.exists(OUTPUT_DIR):
        return 0
    
    from datetime import timedelta
    
    try:
        cutoff_date = datetime.now() - timedelta(days=days_to_keep)
        deleted_count = 0
        
        for folder_name in os.listdir(OUTPUT_DIR):
            folder_path = os.path.join(OUTPUT_DIR, folder_name)
            
            if not os.path.isdir(folder_path):
                continue
            
            try:
                # Parse date from folder name (YYYY-MM-DD)
                folder_date = datetime.strptime(folder_name, "%Y-%m-%d")
                
                if folder_date < cutoff_date:
                    shutil.rmtree(folder_path)
                    logger.info(f"Deleted old output folder: {folder_name}")
                    deleted_count += 1
            
            except ValueError:
                # Not a date folder, skip
                continue
        
        if deleted_count > 0:
            logger.success(f"Cleaned {deleted_count} old output folders")
        
        return deleted_count
    
    except Exception as e:
        logger.log_exception(e, "Error cleaning old outputs")
        return 0