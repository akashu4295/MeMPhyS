"""
Dialog windows for MeMPhyS GUI

Creates modal dialogs for About and Preferences.
Note: These are created on-demand by callbacks, not at startup.
"""

import dearpygui.dearpygui as dpg

from src.config import (
    APP_FULL_NAME,
    APP_SUBTITLE,
    DEVELOPERS,
    FONT_SIZES,
    FONT_PREFERENCES,
    COLORS,
)
from src.utils import get_platform
from src.core import app_state
from src.callbacks import apply_preferences_callback


def create_about_dialog():
    """
    Create About dialog window
    
    This function is called by the menu callback when needed.
    """
    # Check if already exists
    if dpg.does_item_exist("about_window"):
        dpg.configure_item("about_window", show=True)
        dpg.focus_item("about_window")
        return
    
    with dpg.window(
        label="About",
        tag="about_window",
        modal=True,
        width=620,
        height=400,
        no_resize=True,
        pos=(330, 200)
    ):
        # Application title
        dpg.add_text(APP_FULL_NAME, color=COLORS["header"])
        dpg.add_text(APP_SUBTITLE, color=COLORS["header"])
        dpg.add_separator()
        
        # GUI Development section
        dpg.add_text("GUI Development:", color=COLORS["subheader"])
        for dev in DEVELOPERS["gui"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        # Backend/Simulation section
        dpg.add_text("Backend / Simulation:", color=COLORS["subheader"])
        for dev in DEVELOPERS["backend"]:
            dpg.add_text(f"\t{dev}")
        
        dpg.add_spacer(height=6)
        
        # Affiliation section
        dpg.add_text("Affiliation:", color=COLORS["subheader"])
        for affil in DEVELOPERS["affiliation"]:
            dpg.add_text(f"\t{affil}")
        
        dpg.add_spacer(height=12)
        
        # Close button
        dpg.add_button(
            label="Close",
            width=80,
            callback=lambda: dpg.configure_item("about_window", show=False)
        )


def create_preferences_dialog():
    """
    Create Preferences dialog window
    
    This function is called by the menu callback when needed.
    """
    # Check if already exists
    if dpg.does_item_exist("preferences_window"):
        dpg.configure_item("preferences_window", show=True)
        dpg.focus_item("preferences_window")
        return
    
    # Get current settings
    current_font_name = app_state.current_font_name
    current_font_size = app_state.current_font_size
    
    # Get available fonts for platform
    platform = get_platform()
    available_fonts = FONT_PREFERENCES.get(platform, [])
    
    with dpg.window(
        label="Preferences",
        tag="preferences_window",
        modal=True,
        width=450,
        height=250,
        no_resize=True,
        pos=(415, 275)
    ):
        dpg.add_text("Preferences", color=COLORS["header"])
        dpg.add_separator()
        
        # Font family selector
        dpg.add_combo(
            label="Font family",
            items=available_fonts,
            default_value=current_font_name,
            tag="pref_font_family",
            width=200
        )
        
        # Font size selector
        dpg.add_combo(
            label="Font size",
            items=FONT_SIZES,
            default_value=current_font_size,
            tag="pref_font_size",
            width=200
        )
        
        dpg.add_spacer(height=12)
        
        # Buttons
        with dpg.group(horizontal=True):
            dpg.add_button(
                label="Apply",
                width=80,
                callback=apply_preferences_callback
            )
            
            dpg.add_button(
                label="Close",
                width=80,
                callback=lambda: dpg.configure_item("preferences_window", show=False)
            )


def show_error_dialog(title: str, message: str):
    """
    Show an error dialog
    
    Args:
        title: Dialog title
        message: Error message to display
    """
    dialog_tag = f"error_dialog_{id(message)}"
    
    with dpg.window(
        label=title,
        tag=dialog_tag,
        modal=True,
        width=400,
        height=200,
        no_resize=True,
        pos=(440, 300)
    ):
        dpg.add_text("Error:", color=COLORS["error"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        dpg.add_button(
            label="OK",
            width=80,
            callback=lambda: dpg.delete_item(dialog_tag)
        )


def show_info_dialog(title: str, message: str):
    """
    Show an information dialog
    
    Args:
        title: Dialog title
        message: Information message to display
    """
    dialog_tag = f"info_dialog_{id(message)}"
    
    with dpg.window(
        label=title,
        tag=dialog_tag,
        modal=True,
        width=400,
        height=200,
        no_resize=True,
        pos=(440, 300)
    ):
        dpg.add_text("Information:", color=COLORS["info"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        dpg.add_button(
            label="OK",
            width=80,
            callback=lambda: dpg.delete_item(dialog_tag)
        )


def show_confirmation_dialog(title: str, message: str, on_confirm, on_cancel=None):
    """
    Show a confirmation dialog with Yes/No buttons
    
    Args:
        title: Dialog title
        message: Confirmation message
        on_confirm: Callback function to call on Yes
        on_cancel: Optional callback function to call on No
    """
    dialog_tag = f"confirm_dialog_{id(message)}"
    
    def confirm_callback():
        dpg.delete_item(dialog_tag)
        if on_confirm:
            on_confirm()
    
    def cancel_callback():
        dpg.delete_item(dialog_tag)
        if on_cancel:
            on_cancel()
    
    with dpg.window(
        label=title,
        tag=dialog_tag,
        modal=True,
        width=400,
        height=200,
        no_resize=True,
        pos=(440, 300)
    ):
        dpg.add_text("Confirm:", color=COLORS["warning"])
        dpg.add_separator()
        dpg.add_text(message, wrap=380)
        dpg.add_spacer(height=12)
        
        with dpg.group(horizontal=True):
            dpg.add_button(
                label="Yes",
                width=80,
                callback=confirm_callback
            )
            
            dpg.add_button(
                label="No",
                width=80,
                callback=cancel_callback
            )