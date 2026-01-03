"""
Visualization panel construction for MeMPhyS GUI

Creates the right panel with convergence plot, contour plotting controls,
and log window.
"""

import dearpygui.dearpygui as dpg

from src.config import (
    CONVERGENCE_PLOT_HEIGHT,
    LOG_WINDOW_HEIGHT,
    DEFAULT_VTK,
    VTK_PATH_INPUT_WIDTH,
    COMBO_WIDTH,
    SAVE_PATH_INPUT_WIDTH,
    DEFAULT_SAVE_PATH,
    COLORMAPS,
    DEFAULT_COLORMAP,
    PLOT_VARIABLES,
    DEFAULT_PLOT_VARIABLE,
    COLORS,
    FILE_EXTENSIONS,
    FILE_DIALOG_WIDTH,
    FILE_DIALOG_HEIGHT,
)
from src.callbacks import (
    update_plot_callback,
    save_plot_image_callback,
    open_file_dialog_callback,
    select_vtk_file_callback,
    change_colormap_callback,
    change_variable_callback,
)
# from src.core import clear_logs


def create_visualization_panel(themes: dict) -> int:
    """
    Create the right panel with visualization controls
    
    Args:
        themes: Dictionary of theme IDs
    
    Returns:
        Child window tag/ID
    """
    with dpg.child_window(
        width=-1,
        border=True,
        tag="visualization_panel"
    ) as panel:
        
        # Convergence Plot Section
        _create_convergence_plot_section(themes)
        
        dpg.add_separator()
        dpg.add_spacer(height=6)
        
        # Plotting Options Section
        _create_plotting_options_section(themes)
        
        dpg.add_separator()
        dpg.add_spacer(height=6)
        
        # Logs Section
        _create_logs_section(themes)
    
    return panel


def _create_convergence_plot_section(themes: dict):
    """Create convergence plot with controls"""
    dpg.add_text("Convergence Plot", color=COLORS["success"])
    
    # Create plot
    with dpg.plot(
        label="Residual",
        height=CONVERGENCE_PLOT_HEIGHT,
        width=-1,
        tag="convergence_plot"
    ):
        # X-axis
        dpg.add_plot_axis(
            dpg.mvXAxis,
            label="Timestep",
            tag="x_axis_conv"
        )
        
        # Y-axis (log scale)
        dpg.add_plot_axis(
            dpg.mvYAxis,
            label="Total Steady-State Error",
            tag="y_axis_conv",
            log_scale=True
        )
        
        # Line series
        dpg.add_line_series(
            [],
            [],
            parent="y_axis_conv",
            tag="conv_series",
            label="Convergence"
        )
    
    # Apply plot theme
    if "plot" in themes:
        dpg.bind_item_theme("convergence_plot", themes["plot"])


def _create_plotting_options_section(themes: dict):
    """Create contour plotting controls"""
    dpg.add_text("Plotting options", color=COLORS["success"])
    
    # First row: VTK file, variable, colormap
    with dpg.group(horizontal=True):
        dpg.add_text("VTK File:")
        
        dpg.add_input_text(
            hint="VTK File",
            default_value=DEFAULT_VTK,
            tag="contour_vtk_path",
            width=VTK_PATH_INPUT_WIDTH,
            show=True
        )
        
        dpg.add_button(
            label="Browse",
            tag="vtk_browse",
            callback=open_file_dialog_callback,
            user_data="vtk",
            show=True
        )
        
        dpg.add_text("Variable:")
        
        dpg.add_combo(
            PLOT_VARIABLES,
            label="##Variable",
            default_value=DEFAULT_PLOT_VARIABLE,
            tag="contour_var",
            width=COMBO_WIDTH,
            callback=change_variable_callback
        )
        
        dpg.add_text("Colormap:")
        
        dpg.add_combo(
            COLORMAPS,
            label="##Colormap",
            default_value=DEFAULT_COLORMAP,
            tag="contour_cmap",
            width=COMBO_WIDTH,
            callback=change_colormap_callback
        )
        
        plot_btn = dpg.add_button(
            label="Plot",
            callback=update_plot_callback,
            tag="plot_button"
        )
        
        # Apply theme to plot button
        if "button_secondary" in themes:
            dpg.bind_item_theme(plot_btn, themes["button_secondary"])
    
    # # Second row: Save image controls
    # with dpg.group(horizontal=True):
    #     dpg.add_text("Save Image as:")
        
    #     dpg.add_input_text(
    #         hint="Save Path",
    #         default_value=DEFAULT_SAVE_PATH,
    #         tag="contour_save_path",
    #         width=SAVE_PATH_INPUT_WIDTH
    #     )
        
    #     save_btn = dpg.add_button(
    #         label="Save",
    #         callback=save_plot_image_callback,
    #         tag="save_plot_button"
    #     )
        
    #     # Apply theme to save button
    #     if "button_secondary" in themes:
    #         dpg.bind_item_theme(save_btn, themes["button_secondary"])
    
    # File dialog for VTK file
    dpg.add_file_dialog(
        directory_selector=False,
        tag="file_dialog_vtk",
        user_data="contour_vtk_path",
        callback=select_vtk_file_callback,
        show=False,
        width=FILE_DIALOG_WIDTH,
        height=FILE_DIALOG_HEIGHT
    )
    
    # Add .vtk file extension filter
    if ".vtk" in FILE_EXTENSIONS:
        dpg.add_file_extension(
            ".vtk",
            parent="file_dialog_vtk",
            color=FILE_EXTENSIONS[".vtk"]
        )


def _create_logs_section(themes: dict):
    """Create logs display section"""
    dpg.add_text("Logs", color=COLORS["success"])
    
    # Scrollable child window for logs
    with dpg.child_window(
        height=LOG_WINDOW_HEIGHT,
        border=True,
        tag="log_child"
    ):
        dpg.add_input_text(
            tag="log_window",
            multiline=True,
            readonly=True,
            width=-1,
            height=-1,
            default_value="Ready.\n"
        )
    
    # Clear logs button
    clear_log_btn = dpg.add_button(
        label="Clear Logs",
        callback=lambda: clear_logs(),
        tag="clear_logs_button"
    )
    
    # Apply theme to clear button
    if "button_secondary" in themes:
        dpg.bind_item_theme(clear_log_btn, themes["button_secondary"])