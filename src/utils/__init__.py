"""
Shared utilities for MeMPhyS

This package contains file, platform, output, and Gmsh utilities.

Example:
    from src.utils import initialize_fonts, write_parameters_csv
    from src.utils.platform_utils import open_folder, get_platform
"""

# File I/O utilities
from .file_io import (
    validate_file_path,
    validate_mesh_file,
    read_parameters_from_gui,
    write_parameters_csv,
    read_mesh_files_from_gui,
    write_grid_csv,
    read_csv_file,
    write_csv_file,
    create_backup,
    ensure_directory_exists,
    get_file_info,
    load_configuration,
    save_configuration,
    list_files_in_directory,
    read_text_file,
    write_text_file,
)

# Platform utilities
from .platform_utils import (
    get_platform,
    is_windows,
    is_macos,
    is_linux,
    get_platform_info,
    open_folder,
    open_file,
    open_url,
    get_executable_name,
    get_compiler_command,
    build_compile_command,
    get_path_separator,
    normalize_path,
    get_temp_directory,
    get_home_directory,
    check_command_exists,
    get_system_info_summary,
    check_disk_space,
    format_bytes,
    get_environment_variable,
    set_environment_variable,
    run_command,
    is_path_writable,
    is_path_readable,
    kill_process_by_name,
)

# Output management
from .output_manager import (
    get_output_folder_path,
    ensure_output_folder_exists,
    get_unique_filename,
    move_file_to_output,
    organize_solver_outputs,
    list_output_runs,
    get_latest_output_folder,
    clean_old_outputs,
)

# Gmsh and BC management
from .gmsh_bc_manager import (
    launch_gmsh,
    read_physical_names_from_msh,
    get_physical_names,
    set_boundary_condition,
    get_boundary_condition,
    get_all_boundary_conditions,
    clear_boundary_conditions,
    write_bc_csv,
    read_bc_csv,
    get_last_gmsh_file,
    set_last_gmsh_file,
    # get_mesh_from_gmsh_file,
    validate_bc_assignment,
)

__all__ = [
    # File I/O utilities
    'validate_file_path',
    'validate_mesh_file',
    'read_parameters_from_gui',
    'write_parameters_csv',
    'read_mesh_files_from_gui',
    'write_grid_csv',
    'read_csv_file',
    'write_csv_file',
    'create_backup',
    'ensure_directory_exists',
    'get_file_info',
    'load_configuration',
    'save_configuration',
    'list_files_in_directory',
    'read_text_file',
    'write_text_file',
    
    # Platform utilities
    'get_platform',
    'is_windows',
    'is_macos',
    'is_linux',
    'get_platform_info',
    'open_folder',
    'open_file',
    'open_url',
    'get_executable_name',
    'get_compiler_command',
    'build_compile_command',
    'get_path_separator',
    'normalize_path',
    'get_temp_directory',
    'get_home_directory',
    'check_command_exists',
    'get_system_info_summary',
    'check_disk_space',
    'format_bytes',
    'get_environment_variable',
    'set_environment_variable',
    'run_command',
    'is_path_writable',
    'is_path_readable',
    'kill_process_by_name',
    
    # Output management
    'get_output_folder_path',
    'ensure_output_folder_exists',
    'get_unique_filename',
    'move_file_to_output',
    'organize_solver_outputs',
    'list_output_runs',
    'get_latest_output_folder',
    'clean_old_outputs',
]