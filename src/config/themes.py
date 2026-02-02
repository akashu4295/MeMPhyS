"""
DearPyGUI theme definitions for MeMPhyS GUI
Contains all visual styling and theme configurations
"""

import dearpygui.dearpygui as dpg


def create_button_theme():
    """
    Create primary button theme with rounded corners and hover effects
    Used for main action buttons like "Compile and Run Solver"
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            # Colors
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))          # Normal state
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))   # Hover state
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))    # Pressed state
            
            # Styling
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)               # Rounded corners
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)             # Inner padding (x, y)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)             # Spacing between items
    
    return theme


def create_secondary_button_theme():
    """
    Create secondary button theme for less prominent actions
    Used for utility buttons like "Browse", "Save", "Clear Logs"
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            # Colors - slightly different from primary
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))
            
            # Styling
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)
    
    return theme


def create_menu_bar_theme():
    """
    Create menu bar theme with dark background
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvMenuBar):
            # Dark menu bar background
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg, (28, 28, 32))
            
            # Light text color
            dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 220, 220))
            
            # Padding
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)
    
    return theme


def create_input_theme():
    """
    Create theme for input fields
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvInputText):
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 8, 4)
    
    return theme


def create_window_theme():
    """
    Create theme for child windows and panels
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvChildWindow):
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (30, 30, 35))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (60, 60, 70))
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 10, 10)
            dpg.add_theme_style(dpg.mvStyleVar_ChildRounding, 6)
    
    return theme


def create_plot_theme():
    """
    Create theme for plots
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvPlot):
            dpg.add_theme_color(dpg.mvPlotCol_FrameBg, (25, 25, 30))
            dpg.add_theme_color(dpg.mvPlotCol_PlotBg, (20, 20, 25))
            dpg.add_theme_color(dpg.mvPlotCol_PlotBorder, (60, 60, 70))
            dpg.add_theme_style(dpg.mvPlotStyleVar_PlotPadding, 12, 12)
    
    return theme


def create_combo_theme():
    """
    Create theme for combo boxes (dropdowns)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvCombo):
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 40, 45))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (50, 50, 55))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 8, 4)
    
    return theme


def create_checkbox_theme():
    """
    Create theme for checkboxes
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvCheckbox):
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark, (40, 120, 200))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (50, 50, 55))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (70, 70, 75))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 3)
    
    return theme


def create_separator_theme():
    """
    Create theme for separators
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvSeparator):
            dpg.add_theme_color(dpg.mvThemeCol_Separator, (80, 80, 90))
    
    return theme


def create_modal_theme():
    """
    Create theme for modal windows (About, Preferences)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvWindowAppItem):
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (35, 35, 40))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (80, 80, 90))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (28, 28, 32))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (40, 120, 200))
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 8)
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 15, 15)
    
    return theme


def create_disabled_theme():
    """
    Create theme for disabled items
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton, enabled_state=False):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (60, 60, 65))
            dpg.add_theme_color(dpg.mvThemeCol_Text, (120, 120, 125))
    
    return theme


def create_success_button_theme():
    """
    Create theme for success state buttons (green)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 180, 80))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (50, 200, 100))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 160, 70))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
    
    return theme


def create_error_button_theme():
    """
    Create theme for error state buttons (red)
    """
    with dpg.theme() as theme:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button, (200, 60, 60))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (220, 80, 80))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (180, 50, 50))
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)
    
    return theme


def initialize_all_themes():
    """
    Initialize all themes and return them in a dictionary for easy access
    
    Returns:
        dict: Dictionary mapping theme names to theme IDs
    """
    themes = {
        "button": create_button_theme(),
        "button_secondary": create_secondary_button_theme(),
        "button_success": create_success_button_theme(),
        "button_error": create_error_button_theme(),
        "menu_bar": create_menu_bar_theme(),
        "input": create_input_theme(),
        "window": create_window_theme(),
        "plot": create_plot_theme(),
        "combo": create_combo_theme(),
        "checkbox": create_checkbox_theme(),
        "separator": create_separator_theme(),
        "modal": create_modal_theme(),
        "disabled": create_disabled_theme(),
    }
    
    return themes


def apply_global_theme():
    """
    Apply global theme settings that affect all widgets
    """
    with dpg.theme() as global_theme:
        with dpg.theme_component(dpg.mvAll):
            # Global colors - Dark theme
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (25, 25, 30))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (60, 60, 70))
            dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 220, 220))
            
            # Global styles
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 10, 10)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 8, 4)
            dpg.add_theme_style(dpg.mvStyleVar_ItemInnerSpacing, 4, 4)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarSize, 14)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarRounding, 8)
            dpg.add_theme_style(dpg.mvStyleVar_GrabMinSize, 10)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 4)
    
    dpg.bind_theme(global_theme)
    return global_theme


def apply_light_theme():
    """
    Apply light theme settings
    """
    with dpg.theme() as light_theme:
        with dpg.theme_component(dpg.mvAll):
            # Global colors - Light theme
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (240, 240, 245))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (180, 180, 190))
            dpg.add_theme_color(dpg.mvThemeCol_Text, (20, 20, 20))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (255, 255, 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (245, 245, 250))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (235, 235, 245))
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (250, 250, 252))
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg, (230, 230, 235))
            
            # Global styles
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 10, 10)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 8, 4)
            dpg.add_theme_style(dpg.mvStyleVar_ItemInnerSpacing, 4, 4)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarSize, 14)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarRounding, 8)
            dpg.add_theme_style(dpg.mvStyleVar_GrabMinSize, 10)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 4)
    
    dpg.bind_theme(light_theme)
    return light_theme


def toggle_theme(current_dark: bool):
    """
    Toggle between dark and light themes
    
    Args:
        current_dark: True if currently dark theme
        
    Returns:
        New theme state (True = dark, False = light)
    """
    if current_dark:
        apply_light_theme()
        return False
    else:
        apply_global_theme()
        return True