"""
Menu bar construction for MeMPhyS GUI

Creates the application menu bar with File, Edit, and Help menus.
"""

import dearpygui.dearpygui as dpg

from src.callbacks import (
    open_logs_callback,
    open_help_callback,
    show_about_callback,
    show_preferences_callback,
    exit_application_callback,
    open_config_callback,
    save_config_callback,
    show_options_callback,
    load_config_file_callback,
    save_config_file_callback,
    toggle_theme_callback,
)


def create_menu_bar(themes: dict) -> int:
    """
    Create the main menu bar
    
    Args:
        themes: Dictionary of theme IDs
    
    Returns:
        Menu bar tag/ID
    """
    # Create file dialogs first (before menu bar)
    _create_config_file_dialogs()
    
    with dpg.menu_bar(tag="main_menu_bar") as menu_bar:
        # File Menu
        with dpg.menu(label="File"):
            dpg.add_menu_item(
                label="Open Configuration",
                callback=open_config_callback
            )
            
            dpg.add_menu_item(
                label="Open Logs",
                callback=open_logs_callback
            )
            
            dpg.add_separator()
            
            dpg.add_menu_item(
                label="Save Configuration",
                callback=save_config_callback
            )
            
            dpg.add_separator()
            
            dpg.add_menu_item(
                label="Exit",
                callback=exit_application_callback
            )
        
        # Spacer to push remaining items
        dpg.add_spacer(width=20)
        
        # Edit Menu
        with dpg.menu(label="Edit"):
            dpg.add_menu_item(
                label="Preferences",
                callback=show_preferences_callback
            )

        # Spacer to push remaining items
        dpg.add_spacer(width=20)

        # Options Menu
        with dpg.menu(label="Options"):
            dpg.add_menu_item(
                label="Application Options",
                callback=show_options_callback
            )
        
        # Spacer to push remaining items
        dpg.add_spacer(width=20)
        
        # View Menu
        with dpg.menu(label="View"):
            dpg.add_menu_item(
                label="Toggle Dark/Light Theme",
                callback=toggle_theme_callback
            )
        
        # Spacer to push Help to the right
        dpg.add_spacer(width=20)
        
        # Help Menu
        with dpg.menu(label="Help"):
            dpg.add_menu_item(
                label="Documentation",
                callback=open_help_callback
            )
            
            dpg.add_separator()
            
            dpg.add_menu_item(
                label="About",
                callback=show_about_callback
            )
    
    # Apply menu bar theme
    if "menu_bar" in themes:
        dpg.bind_item_theme(menu_bar, themes["menu_bar"])
    
    return menu_bar


def _create_config_file_dialogs():
    """Create file dialogs for configuration load/save"""
    from src.config import FILE_DIALOG_WIDTH, FILE_DIALOG_HEIGHT
    
    print("Creating config file dialogs...")
    
    # Load configuration dialog
    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_config_load",
        callback=load_config_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT,
        default_path="./config"
    )
    dpg.add_file_extension(".json", parent="file_dialog_config_load", color=(255, 255, 0, 255))
    print("Created config load dialog")
    
    # Save configuration dialog
    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_config_save",
        callback=save_config_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT,
        default_path="./config",
        default_filename="my_config.json"
    )
    dpg.add_file_extension(".json", parent="file_dialog_config_save", color=(255, 255, 0, 255))
    print("Created config save dialog")