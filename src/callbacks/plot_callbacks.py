"""
Visualization and plotting callbacks for MeMPhyS GUI

Handles callbacks for:
- Contour plotting
- Plot updates
- Image saving
- Convergence plot control
"""

import os
import sys
import subprocess
import dearpygui.dearpygui as dpg

from src.core import logger, app_state
from src.solver import convergence_monitor
from src.config import DEFAULT_VTK, DEFAULT_SAVE_PATH
from src.utils import validate_file_path


def update_plot_callback(sender, app_data, user_data):
    """
    Update the contour plot with current settings
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # Get plot settings
    vtk_path = dpg.get_value("contour_vtk_path") or DEFAULT_VTK
    var_choice = dpg.get_value("contour_var") or "velocity magnitude"
    cmap = dpg.get_value("contour_cmap") or "viridis"
    
    # Get dimension from parameters
    if dpg.does_item_exist("param_domain_dimensions"):
        dimension = dpg.get_value("param_domain_dimensions")
    else:
        dimension = "2"
    
    # Validate VTK file exists
    if not validate_file_path(vtk_path, must_exist=True):
        logger.error(f"VTK file not found: {vtk_path}")
        logger.info("Please ensure the solver has run and generated output")
        return
    
    logger.info(f"Plotting {var_choice} from {vtk_path}")
    logger.info(f"Colormap: {cmap}, Dimension: {dimension}")
    
    # Get path to plotter script
    script_path = os.path.join("src", "plotter.py")
    
    if not os.path.exists(script_path):
        logger.error(f"Plotter script not found: {script_path}")
        return
    
    # Terminate existing plotter process if running
    if app_state.plotter_process is not None:
        try:
            if app_state.plotter_process.poll() is None:
                app_state.plotter_process.terminate()
                app_state.plotter_process.wait(timeout=2)
                logger.info("Previous plotter window closed")
        except Exception as e:
            logger.warning(f"Error closing previous plotter: {e}")
    
    # Launch new plotter process
    try:
        proc = subprocess.Popen(
            [
                sys.executable,
                script_path,
                vtk_path,
                var_choice,
                cmap,
                str(dimension),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        
        app_state.plotter_process = proc
        logger.success("Plotter window opened")
    
    except Exception as e:
        logger.log_exception(e, "Error launching plotter")


def save_plot_image_callback(sender, app_data, user_data):
    """
    Save plot image (currently shows message about using external plotter)
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    save_path = dpg.get_value("contour_save_path") or DEFAULT_SAVE_PATH
    
    logger.info("Screenshot feature requires the external plotting window")
    logger.info("Please use the plot window's built-in save/screenshot functionality")
    logger.info(f"Suggested filename: {save_path}")


def reset_convergence_plot_callback(sender, app_data, user_data):
    """
    Reset the convergence plot
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    convergence_monitor.reset()
    logger.info("Convergence plot reset")


def force_convergence_update_callback(sender, app_data, user_data):
    """
    Force an immediate update of the convergence plot
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    convergence_monitor.force_update()
    logger.info("Convergence plot updated")


def toggle_convergence_monitor_callback(sender, app_data, user_data):
    """
    Toggle convergence monitoring on/off
    
    Args:
        sender: Checkbox tag
        app_data: Checkbox state
        user_data: User data (unused)
    """
    enabled = app_data
    
    if enabled:
        if not convergence_monitor.is_running():
            convergence_monitor.start()
            logger.info("Convergence monitoring enabled")
    else:
        if convergence_monitor.is_running():
            convergence_monitor.stop()
            logger.info("Convergence monitoring disabled")


def show_convergence_status_callback(sender, app_data, user_data):
    """
    Display detailed convergence status
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    status = convergence_monitor.get_convergence_status()
    
    logger.separator()
    logger.info("Convergence Status:")
    
    if status["has_data"]:
        logger.info(f"  Iterations: {status['iterations']}")
        
        if status['current_error'] is not None:
            logger.info(f"  Current Error: {status['current_error']:.6e}")
        
        if status['min_error'] is not None:
            logger.info(f"  Minimum Error: {status['min_error']:.6e}")
        
        converging_str = "Yes" if status.get('is_converging', False) else "No"
        logger.info(f"  Converging: {converging_str}")
    else:
        logger.info("  No convergence data available")
    
    logger.separator()


def export_convergence_data_callback(sender, app_data, user_data):
    """
    Export convergence data to a file
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    data = convergence_monitor.get_latest_data()
    
    if data is None:
        logger.warning("No convergence data available to export")
        return
    
    x, y = data
    
    # Create export filename with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"convergence_export_{timestamp}.csv"
    
    try:
        import pandas as pd
        df = pd.DataFrame({"Iteration": x, "Error": y})
        df.to_csv(filename, index=False)
        
        logger.success(f"Convergence data exported to {filename}")
        logger.info(f"Exported {len(x)} data points")
    
    except Exception as e:
        logger.log_exception(e, "Error exporting convergence data")


def browse_vtk_file_callback(sender, app_data, user_data):
    """
    Open file dialog for VTK file selection
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    if dpg.does_item_exist("file_dialog_vtk"):
        dpg.configure_item("file_dialog_vtk", show=True)


def select_vtk_file_callback(sender, app_data, user_data):
    """
    Callback when VTK file is selected
    
    Args:
        sender: File dialog tag
        app_data: Dictionary with file info
        user_data: Target input field tag
    """
    if 'file_path_name' not in app_data:
        return
    
    file_path = app_data['file_path_name']
    
    if dpg.does_item_exist("contour_vtk_path"):
        dpg.set_value("contour_vtk_path", file_path)
        app_state.current_vtk_file = file_path
        logger.info(f"VTK file selected: {file_path}")


def change_colormap_callback(sender, app_data, user_data):
    """
    Callback when colormap is changed
    
    Args:
        sender: Combo box tag
        app_data: Selected colormap
        user_data: User data (unused)
    """
    app_state.current_colormap = app_data
    logger.debug(f"Colormap changed to: {app_data}")


def change_variable_callback(sender, app_data, user_data):
    """
    Callback when plot variable is changed
    
    Args:
        sender: Combo box tag
        app_data: Selected variable
        user_data: User data (unused)
    """
    app_state.current_variable = app_data
    logger.debug(f"Plot variable changed to: {app_data}")


def open_plot_settings_callback(sender, app_data, user_data):
    """
    Open a dialog with advanced plot settings
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    # TODO: Implement advanced plot settings dialog
    logger.info("Advanced plot settings not yet implemented")


def close_all_plot_windows_callback(sender, app_data, user_data):
    """
    Close all external plot windows
    
    Args:
        sender: Button tag
        app_data: Application data (unused)
        user_data: User data (unused)
    """
    closed = False
    
    if app_state.plotter_process is not None:
        try:
            if app_state.plotter_process.poll() is None:
                app_state.plotter_process.terminate()
                app_state.plotter_process.wait(timeout=2)
                logger.info("Plotter window closed")
                closed = True
        except Exception as e:
            logger.warning(f"Error closing plotter: {e}")
        finally:
            app_state.plotter_process = None
    
    if not closed:
        logger.info("No plot windows to close")