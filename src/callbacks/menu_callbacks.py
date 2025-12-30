"""
Menu bar callbacks for MeMPhyS GUI

Handles callbacks for:
- File menu (open, save, exit)
- Edit menu (preferences)
- Help menu (help, about)
"""

import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.utils import open_folder, open_url, change_font, get_platform
from src.config import (
    LOG_DIR,
    HELP_URL,
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    FONT_SIZES,
    FONT_PREFERENCES,
)


def open_logs_callback(sender, app_data, user_data):
    """
    Open the logs directory
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    success = open_folder(LOG_DIR)
    
    if not success:
        logger.error(f"Failed to open logs directory: {LOG_DIR}")


def open_help_callback(sender, app_data, user_data):
    """
    Open the help documentation in web browser
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    success = open_url(HELP_URL)
    
    if success:
        logger.info("Opened help documentation in browser")
    else:
        logger.error("Failed to open help URL")


def show_about_callback(sender, app_data, user_data):
    """
    Show the About dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # If window already exists, just show it
    if dpg.does_item_exist("about_window"):
        dpg.configure_item("about_window", show=True)
        dpg.focus_item("about_window")
        return
    
    # Create About window
    with dpg.window(
        label="About",
        tag="about_window",
        modal=True,
        width=620,
        height=400,
        no_resize=True,
        pos=(330, 200)
    ):
        dpg.add_text(APP_FULL_NAME, color=(200, 220, 255))
        dpg.add_text(APP_SUBTITLE, color=(200, 220, 255))
        dpg.add_separator()
        
        dpg.add_text("GUI Development:", color=(255, 220, 160))
        for dev in DEVELOPERS["gui"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        dpg.add_text("Backend / Simulation:", color=(255, 220, 160))
        for dev in DEVELOPERS["backend"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        dpg.add_text("Affiliation:", color=(255, 220, 160))
        for affil in DEVELOPERS["affiliation"]:
            dpg.add_text(f"\t{affil}")
        
        dpg.add_spacer(height=12)
        
        close_btn = dpg.add_button(
            label="Close",
            width=80,
            callback=lambda: dpg.configure_item("about_window", show=False)
        )


def show_preferences_callback(sender, app_data, user_data):
    """
    Show the Preferences dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # If window already exists, just show it
    if dpg.does_item_exist("preferences_window"):
        dpg.configure_item("preferences_window", show=True)
        dpg.focus_item("preferences_window")
        return
    
    # Get current font settings
    current_font_name = app_state.current_font_name
    current_font_size = app_state.current_font_size
    
    # Get available fonts for current platform
    platform = get_platform()
    available_fonts = FONT_PREFERENCES.get(platform, [])
    
    # Create Preferences window
    with dpg.window(
        label="Preferences",
        tag="preferences_window",
        modal=True,
        width=450,
        height=250,
        no_resize=True,
        pos=(415, 275)
    ):
        dpg.add_text("Preferences", color=(200, 220, 255))
        dpg.add_separator()
        
        dpg.add_combo(
            label="Font family",
            items=available_fonts,
            default_value=current_font_name,
            tag="pref_font_family",
            width=200
        )
        
        dpg.add_combo(
            label="Font size",
            items=FONT_SIZES,
            default_value=current_font_size,
            tag="pref_font_size",
            width=200
        )
        
        dpg.add_spacer(height=12)
        
        dpg.add_button(
            label="Apply",
            width=80,
            callback=apply_preferences_callback
        )
        
        dpg.add_same_line()
        
        dpg.add_button(
            label="Close",
            width=80,
            callback=lambda: dpg.configure_item("preferences_window", show=False)
        )


def apply_preferences_callback(sender, app_data, user_data):
    """
    Apply preferences changes
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if not dpg.does_item_exist("pref_font_family") or not dpg.does_item_exist("pref_font_size"):
        logger.error("Preference widgets not found")
        return
    
    font_name = dpg.get_value("pref_font_family")
    font_size = dpg.get_value("pref_font_size")
    
    logger.info(f"Attempting to change font to {font_name} @ {font_size}px")
    
    # Try to change font
    success = change_font(font_name, font_size)
    
    if success:
        logger.success(f"Font changed to {font_name} @ {font_size}px")
        logger.info("Note: Some UI elements may require application restart to update")
    else:
        logger.error("Failed to change font")
        logger.info("Font may not be available. Try another font or size.")


def exit_application_callback(sender, app_data, user_data):
    """
    Exit the application with cleanup
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Exiting application...")
    
    # Perform cleanup
    app_state.cleanup()
    
    # Stop DearPyGUI
    dpg.stop_dearpygui()


def open_config_callback(sender, app_data, user_data):
    """
    Open a configuration file (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Configuration loading not yet implemented")
    
    if dpg.does_item_exist("file_dialog_config"):
        dpg.configure_item("file_dialog_config", show=True)


def save_config_callback(sender, app_data, user_data):
    """
    Save current configuration (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Configuration saving not yet implemented")
    
    if dpg.does_item_exist("file_dialog_config_save"):
        dpg.configure_item("file_dialog_config_save", show=True)


def show_system_info_callback(sender, app_data, user_data):
    """
    Show system information dialog
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    from src.utils import get_system_info_summary
    
    info = get_system_info_summary()
    
    logger.separator()
    logger.info("System Information:")
    for line in info.split('\n'):
        if line.strip():
            logger.info(line)
    logger.separator()


def show_application_state_callback(sender, app_data, user_data):
    """
    Show current application state (for debugging)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    state_summary = app_state.get_state_summary()
    
    logger.separator()
    logger.info("Application State:")
    for key, value in state_summary.items():
        logger.info(f"  {key}: {value}")
    logger.separator()


def check_for_updates_callback(sender, app_data, user_data):
    """
    Check for application updates (placeholder)
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    logger.info("Update checking not yet implemented")
    logger.info(f"Current version: {APP_FULL_NAME}")


def report_bug_callback(sender, app_data, user_data):
    """
    Open bug report page
    
    Args:
        sender: Menu item tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    bug_url = f"{HELP_URL}/issues"
    success = open_url(bug_url)
    
    if success:
        logger.info("Opened bug report page in browser")
    else:
        logger.error("Failed to open bug report URL")
        logger.info(f"Please visit: {bug_url}")