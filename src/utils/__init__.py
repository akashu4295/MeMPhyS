"""
Utilities package

This package contains utility modules for:
- Font management (fonts.py)
- File I/O operations (file_io.py)
- Platform-specific utilities (platform_utils.py)

"""

# Font utilities
from .fonts import (
    find_system_font,
    load_font,
    initialize_fonts,
    set_default_font,
    change_font,
    get_available_fonts,
    get_available_sizes,
    get_current_font,
    list_system_fonts,
    validate_font_file,
    get_font_info,
    create_font_with_fallback,
)

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

__all__ = [
    # Font utilities
    'find_system_font',
    'load_font',
    'initialize_fonts',
    'set_default_font',
    'change_font',
    'get_available_fonts',
    'get_available_sizes',
    'get_current_font',
    'list_system_fonts',
    'validate_font_file',
    'get_font_info',
    'create_font_with_fallback',
    
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
]