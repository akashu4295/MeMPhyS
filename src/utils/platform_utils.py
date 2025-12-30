"""
Platform-specific utilities for MeMPhyS GUI

Handles OS-specific operations like opening folders, detecting OS,
and platform-specific file paths.
"""

import os
import sys
import platform
import subprocess
import webbrowser
from pathlib import Path
from typing import Optional, List, Tuple

from src.core import logger


def get_platform() -> str:
    """
    Get the current platform
    
    Returns:
        Platform name: "Windows", "Darwin" (macOS), or "Linux"
    """
    return platform.system()


def is_windows() -> bool:
    """Check if running on Windows"""
    return get_platform() == "Windows"


def is_macos() -> bool:
    """Check if running on macOS"""
    return get_platform() == "Darwin"


def is_linux() -> bool:
    """Check if running on Linux"""
    return get_platform() == "Linux"


def get_platform_info() -> dict:
    """
    Get detailed platform information
    
    Returns:
        Dictionary with platform details
    """
    return {
        "system": platform.system(),
        "release": platform.release(),
        "version": platform.version(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python_version": platform.python_version(),
    }


def open_folder(path: str) -> bool:
    """
    Open a folder in the system file explorer
    
    Args:
        path: Path to folder to open
    
    Returns:
        True if successful, False otherwise
    """
    try:
        path = os.path.abspath(path)
        
        if not os.path.exists(path):
            logger.error(f"Folder does not exist: {path}")
            return False
        
        system = get_platform()
        
        if system == "Windows":
            os.startfile(path)
        elif system == "Darwin":
            subprocess.Popen(["open", path])
        else:  # Linux and others
            subprocess.Popen(["xdg-open", path])
        
        logger.info(f"Opened folder: {path}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error opening folder {path}")
        return False


def open_file(filepath: str) -> bool:
    """
    Open a file with the system default application
    
    Args:
        filepath: Path to file to open
    
    Returns:
        True if successful, False otherwise
    """
    try:
        filepath = os.path.abspath(filepath)
        
        if not os.path.exists(filepath):
            logger.error(f"File does not exist: {filepath}")
            return False
        
        system = get_platform()
        
        if system == "Windows":
            os.startfile(filepath)
        elif system == "Darwin":
            subprocess.Popen(["open", filepath])
        else:  # Linux
            subprocess.Popen(["xdg-open", filepath])
        
        logger.info(f"Opened file: {filepath}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error opening file {filepath}")
        return False


def open_url(url: str) -> bool:
    """
    Open a URL in the default web browser
    
    Args:
        url: URL to open
    
    Returns:
        True if successful, False otherwise
    """
    try:
        webbrowser.open(url)
        logger.info(f"Opened URL: {url}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error opening URL {url}")
        return False


def get_executable_name(base_name: str) -> str:
    """
    Get the appropriate executable name for the current platform
    
    Args:
        base_name: Base name of the executable (without extension)
    
    Returns:
        Full executable name with platform-specific extension
    """
    if is_windows():
        return f"{base_name}.exe"
    else:
        return base_name


def get_compiler_command() -> str:
    """
    Get the appropriate C compiler command for the platform
    
    Returns:
        Compiler command (typically "gcc")
    """
    return "gcc"


def build_compile_command(source_files: List[str], output_name: str, 
                         flags: Optional[List[str]] = None) -> List[str]:
    """
    Build a compilation command for the current platform
    
    Args:
        source_files: List of source file paths
        output_name: Name of output executable
        flags: Additional compiler flags
    
    Returns:
        List of command arguments
    """
    compiler = get_compiler_command()
    output_exe = get_executable_name(output_name)
    
    cmd = [compiler] + source_files
    
    if flags:
        cmd.extend(flags)
    
    cmd.extend(["-o", output_exe])
    
    return cmd


def get_path_separator() -> str:
    """
    Get the path separator for the current platform
    
    Returns:
        Path separator character
    """
    return os.sep


def normalize_path(path: str) -> str:
    """
    Normalize a path for the current platform
    
    Args:
        path: Path to normalize
    
    Returns:
        Normalized path
    """
    return os.path.normpath(path)


def get_temp_directory() -> str:
    """
    Get the system temporary directory
    
    Returns:
        Path to temp directory
    """
    import tempfile
    return tempfile.gettempdir()


def get_home_directory() -> str:
    """
    Get the user's home directory
    
    Returns:
        Path to home directory
    """
    return str(Path.home())


def check_command_exists(command: str) -> bool:
    """
    Check if a command exists in the system PATH
    
    Args:
        command: Command to check
    
    Returns:
        True if command exists, False otherwise
    """
    try:
        if is_windows():
            result = subprocess.run(["where", command], 
                                  capture_output=True, 
                                  text=True)
        else:
            result = subprocess.run(["which", command], 
                                  capture_output=True, 
                                  text=True)
        
        return result.returncode == 0
    
    except Exception:
        return False


def get_system_info_summary() -> str:
    """
    Get a formatted string with system information
    
    Returns:
        Formatted system information
    """
    info = get_platform_info()
    
    lines = [
        "=== System Information ===",
        f"OS: {info['system']} {info['release']}",
        f"Architecture: {info['machine']}",
        f"Python: {info['python_version']}",
        f"Processor: {info['processor']}",
    ]
    
    return "\n".join(lines)


def check_disk_space(path: str) -> Optional[Tuple[int, int, int]]:
    """
    Check available disk space at a given path
    
    Args:
        path: Path to check
    
    Returns:
        Tuple of (total, used, free) in bytes, or None if check fails
    """
    try:
        import shutil
        total, used, free = shutil.disk_usage(path)
        return (total, used, free)
    
    except Exception as e:
        logger.log_exception(e, f"Error checking disk space for {path}")
        return None


def format_bytes(bytes_value: int) -> str:
    """
    Format bytes into human-readable string
    
    Args:
        bytes_value: Number of bytes
    
    Returns:
        Formatted string (e.g., "1.5 GB")
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_value < 1024.0:
            return f"{bytes_value:.2f} {unit}"
        bytes_value /= 1024.0
    return f"{bytes_value:.2f} PB"


def get_environment_variable(var_name: str, default: Optional[str] = None) -> Optional[str]:
    """
    Get an environment variable
    
    Args:
        var_name: Variable name
        default: Default value if not found
    
    Returns:
        Variable value or default
    """
    return os.environ.get(var_name, default)


def set_environment_variable(var_name: str, value: str) -> bool:
    """
    Set an environment variable (for current process only)
    
    Args:
        var_name: Variable name
        value: Variable value
    
    Returns:
        True if successful, False otherwise
    """
    try:
        os.environ[var_name] = value
        return True
    except Exception as e:
        logger.log_exception(e, f"Error setting environment variable {var_name}")
        return False


def run_command(command: List[str], cwd: Optional[str] = None,
                capture_output: bool = True) -> Tuple[int, str, str]:
    """
    Run a shell command
    
    Args:
        command: Command and arguments as list
        cwd: Working directory
        capture_output: Whether to capture stdout/stderr
    
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            capture_output=capture_output,
            text=True
        )
        
        return (result.returncode, result.stdout, result.stderr)
    
    except Exception as e:
        logger.log_exception(e, f"Error running command: {' '.join(command)}")
        return (-1, "", str(e))


def is_path_writable(path: str) -> bool:
    """
    Check if a path is writable
    
    Args:
        path: Path to check
    
    Returns:
        True if writable, False otherwise
    """
    try:
        return os.access(path, os.W_OK)
    except Exception:
        return False


def is_path_readable(path: str) -> bool:
    """
    Check if a path is readable
    
    Args:
        path: Path to check
    
    Returns:
        True if readable, False otherwise
    """
    try:
        return os.access(path, os.R_OK)
    except Exception:
        return False


def kill_process_by_name(process_name: str) -> bool:
    """
    Kill a process by name (use with caution)
    
    Args:
        process_name: Name of process to kill
    
    Returns:
        True if successful, False otherwise
    """
    try:
        system = get_platform()
        
        if system == "Windows":
            subprocess.run(["taskkill", "/F", "/IM", process_name], 
                         capture_output=True)
        else:
            subprocess.run(["pkill", process_name], 
                         capture_output=True)
        
        logger.warning(f"Attempted to kill process: {process_name}")
        return True
    
    except Exception as e:
        logger.log_exception(e, f"Error killing process {process_name}")
        return False