"""
UI package for MeMPhyS GUI

This package contains all UI construction modules:
- main_window: Main application window assembly
- menu_bar: Menu bar creation
- parameters_panel: Left panel with parameters
- visualization_panel: Right panel with plots and logs
- dialogs: Modal dialogs (About, Preferences, etc.)

Example:
    from src.ui import create_main_window, initialize_ui
    
    # Initialize and create UI
    themes = initialize_ui()
    main_window = create_main_window(themes)
"""

from .main_window import create_main_window
from .menu_bar import create_menu_bar
from .parameters_panel import create_parameters_panel
from .visualization_panel import create_visualization_panel
from .dialogs import (
    create_about_dialog,
    create_preferences_dialog,
    create_options_dialog,
    show_error_dialog,
    show_info_dialog,
    show_confirmation_dialog,
    apply_options_callback,
    reset_options_to_defaults,
    # create_bc_dialog,
)

__all__ = [
    'create_main_window',
    'create_menu_bar',
    'create_parameters_panel',
    'create_visualization_panel',
    'create_about_dialog',
    'create_preferences_dialog',
    'create_options_dialog',
    'show_error_dialog',
    'show_info_dialog',
    'show_confirmation_dialog',
    'apply_options_callback',
    'reset_options_to_defaults',
    # 'create_bc_dialog',
]