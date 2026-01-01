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
)


def create_menu_bar(themes: dict) -> int:
    """
    Create the main menu bar
    
    Args:
        themes: Dictionary of theme IDs
    
    Returns:
        Menu bar tag/ID
    """
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
        
        # Options Menu
        with dpg.menu(label="Options"):
            dpg.add_menu_item(
                label="Application Options",
                callback=show_options_callback
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