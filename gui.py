import os
import time
import threading
import numpy as np
import pandas as pd
import dearpygui.dearpygui as dpg
import json
import socket
import csv
import subprocess
import sys
import platform
import glob
import pyvista as pv
from pyvistaqt import BackgroundPlotter
import matplotlib.pyplot as plt
import io
import functools
import webbrowser
from datetime import datetime

# ============================================================
# Default parameters for the GUI

STATE = {"plotter": None}
LOG_DIR = "./logs"
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)

today_str = datetime.today().strftime("%Y-%m-%d")
log_file = f"logs/log_{today_str}.txt"
if not os.path.exists(log_file):
    file = open(log_file, "w")
else:
    file = open(log_file, "a")
file.write("DATE: " + today_str + "\n\n")

BASE_PARAMETERS = {
    "domain_dimensions": 2,
    "poly_deg": 3,
    "phs_deg": 3,
    "cloud_size_multiplier": 2,
    "test_derivative": 0,
    "courant_number": 0.3,
    "steady_tolerance": 1e-8,
    "poisson_solver_tolerance": 1e-8,
    "sor_parameter": 1.6,
    "time_step": 0.1,
    "num_time_steps": 100000,
    "write_interval": 50,
    "Re": 10,
}

MULTIGRID_PARAMETERS = {
    "num_vcycles": 10,
    "num_relax": 50,
}

IMPLICIT_PARAMETERS = {
    "iter_momentum": 5,
    "iter_timple": 1,
}

DEFAULT_VTK = "Solution.vtk"
CMAPS = sorted([m for m in plt.colormaps()]) # Gather available matplotlib colormaps for dropdown


# ============================================================
# Helper functions for processing, and CSV Writers

def write_params_csv(filename="flow_parameters.csv"):
    rows = []
    for pname, default in BASE_PARAMETERS.items():
        val = dpg.get_value(f"param_{pname}")
        rows.append([pname, val])
    for pname, default in IMPLICIT_PARAMETERS.items():
        val = dpg.get_value(f"param_{pname}")
        rows.append([pname, val])
    for pname, default in MULTIGRID_PARAMETERS.items():
        val = dpg.get_value(f"param_{pname}")
        rows.append([pname, val])
    rows.append(["neumann_flag_boundary", "1"])
    rows.append(["facRe", "1"])
    rows.append(["facdt", "1"])
    rows.append(["fractional_step", 1 if dpg.get_value("solver_method") == "Fractional Step" else 0])
    pd.DataFrame(rows, columns=["Parameter", "Value"]).to_csv(filename, index=False, header=False)
    append_log("Parameters written to flow_parameters.csv\n")
    append_log_file(file, "Parameters written to flow_parameters.csv\n")


def write_grid_csv():
    """Writes mesh file paths to grid_files.csv"""
    mg_enabled = dpg.get_value("multigrid_toggle")
    num_levels = 1
    mesh_paths = []
    num_levels = int(dpg.get_value("num_mesh_levels"))
    for i in range(num_levels):
        path = dpg.get_value(f"mesh_file_{i+1}")
        if path:
            mesh_paths.append(path)
    with open("grid_filenames.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["num_levels", num_levels])
        for path in mesh_paths:
            writer.writerow([path])
    append_log("Grid details written to grid_filenames.csv\n")
    append_log_file(file, "Grid details written to grid_filenames.csv\n")

# =============================================================
# Function to stream logs to GUI

def append_log_file(file, msg):
    file.write(msg + "\n")
    file.flush()

def append_log(msg):
    current = dpg.get_value("log_window") or ""
    dpg.set_value("log_window", current + msg + "\n")
    dpg.split_frame()
    try:
        # Scroll parent child window (replace 'log_child' with your actual ID)
        max_scroll = dpg.get_y_scroll_max("log_child")
        current_scroll = dpg.get_y_scroll("log_child")
        if abs(max_scroll - current_scroll) < 5:
            dpg.set_y_scroll("log_child", max_scroll)
    except Exception:
        pass

def clear_logs():
    dpg.set_value("log_window", "")
    dpg.set_y_scroll("log_child", 0)

# ============================================================
# Callback function to compile and run the C solver

def run_solver(sender):
    """Compile and run the C solver, streaming logs into the GUI."""
    # disable button while running
    write_grid_csv()
    write_params_csv()
    dpg.disable_item(sender)
    dpg.configure_item(sender, label="Compiling and Running...")
    
    dpg.set_value("log_window", "Starting new run...\n")
    append_log("Saving Meshfile details and parameters...\n")
    append_log_file(file, "Saving Meshfile details and parameters...\n")
    write_grid_csv()
    # Detect OS
    system_type = platform.system()
    compiler_cmd = []
    run_cmd = []

    init_path = dpg.get_value("init_path")    # Path to the init .c file
    header_dir = "header_files"
    header_c_files = glob.glob(os.path.join(header_dir, "*.c")) # Collect all .c files in the 'header_files' directory

    # Build the compile command
    if system_type == "Windows":
        compiler_cmd = ["gcc"] + header_c_files + [init_path, "mg_NS_solver.c", "-lm", "-o", "solver.exe"]
        run_cmd = ["solver.exe"]
    else:
        compiler_cmd = ["gcc"] + header_c_files + [init_path, "mg_NS_solver.c", "-lm", "-o", "solver"]
        run_cmd = ["./solver"]

    append_log("Compiling solver...")
    append_log_file(file, "Compiling solver...")

    # Compile first
    compile_process = subprocess.Popen(
        compiler_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    out, err = compile_process.communicate()

    if compile_process.returncode != 0:
        append_log("Compilation failed:\n" + err)
        append_log_file(file, "Compilation failed:\n" + err)
        dpg.configure_item(sender, label="Compilation Failed! Try Again")
        dpg.enable_item(sender)
        return
    else:
        append_log("Compilation successful. Running solver...\n")
        append_log_file(file, "Compilation successful. Running solver...\n")

    # Run solver asynchronously (so GUI stays responsive)
    def solver_thread():
        process = subprocess.Popen(
            run_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        for line in process.stdout:
            append_log(line.strip())
        for line in process.stderr:
            append_log(line.strip())
        process.wait()

        if process.returncode == 0:
            append_log("Solver completed successfully.\n Solution file saved as Solution.csv")
            append_log_file(file, "Solver completed successfully.\n Solution file saved as Solution.csv")
        else:
            append_log(f"Solver exited with code {process.returncode}.")
            append_log_file(file, f"Solver exited with code {process.returncode}.")

        dpg.configure_item(sender, label="Done! Run Again")
        dpg.enable_item(sender)

    threading.Thread(target=solver_thread, daemon=True).start()

def on_write_callback(sender, app_data, user_data):
    write_params_csv()
    write_grid_csv()
    dpg.set_value("log_window", "Files written!\n")
    time.sleep(0.5)

def show_multigrid_callback():
    multigrid_on = dpg.get_value("multigrid_toggle")

    # show/hide entire multigrid parameter group
    if dpg.does_item_exist("multigrid_parameters_section"):
        dpg.configure_item("multigrid_parameters_section", show=multigrid_on)

    # show/hide mesh groups: if multigrid_on -> show up to num_mesh_levels, else show only first
    if dpg.does_item_exist("num_mesh_levels"):
        n = int(dpg.get_value("num_mesh_levels"))
    else:
        n = 1

    if multigrid_on:
        # show first n mesh groups and browse buttons/dialogs
        for i in range(1, 11):
            show_i = (i <= n)
            if dpg.does_item_exist(f"mesh_group_{i}"):
                dpg.configure_item(f"mesh_group_{i}", show=show_i)
            if dpg.does_item_exist(f"mesh_file_{i}"):
                dpg.configure_item(f"mesh_file_{i}", show=show_i)
            if dpg.does_item_exist(f"browse_{i}"):
                dpg.configure_item(f"browse_{i}", show=show_i)
            if dpg.does_item_exist(f"file_dialog_{i}"):
                dpg.configure_item(f"file_dialog_{i}", show=False)  # keep dialog hidden until browse pressed
    else:
        # hide all except the first
        for i in range(1, 11):
            show_i = (i == 1)
            if dpg.does_item_exist(f"mesh_group_{i}"):
                dpg.configure_item(f"mesh_group_{i}", show=show_i)
            if dpg.does_item_exist(f"mesh_file_{i}"):
                dpg.configure_item(f"mesh_file_{i}", show=show_i)
            if dpg.does_item_exist(f"browse_{i}"):
                dpg.configure_item(f"browse_{i}", show=show_i)
            if dpg.does_item_exist(f"file_dialog_{i}"):
                dpg.configure_item(f"file_dialog_{i}", show=False)

    if multigrid_on:
        update_mesh_inputs()
    
def select_mesh_file(sender, app_data, user_data):
    # user_data contains the tag of the input_text to update
    file_path = app_data['file_path_name']
    dpg.set_value(user_data, file_path)

def update_mesh_inputs():
    """Show the correct number of mesh input groups when num_mesh_levels changes."""
    num_levels = int(dpg.get_value("num_mesh_levels"))
    multigrid_on = dpg.get_value("multigrid_toggle")

    for i in range(1, 11):
        show_i = multigrid_on and (i <= num_levels)
        if dpg.does_item_exist(f"mesh_group_{i}"):
            dpg.configure_item(f"mesh_group_{i}", show=show_i)
        if dpg.does_item_exist(f"mesh_file_{i}"):
            dpg.configure_item(f"mesh_file_{i}", show=show_i)
        if dpg.does_item_exist(f"browse_{i}"):
            dpg.configure_item(f"browse_{i}", show=show_i)
        if dpg.does_item_exist(f"file_dialog_{i}"):
            dpg.configure_item(f"file_dialog_{i}", show=False)  # dialogs hidden until needed

# ============================================================
# Menu bar callbacks

def open_folder(path):
    path = os.path.abspath(path)
    if sys.platform.startswith("win"):
        os.startfile(path)
    elif sys.platform.startswith("darwin"):
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])

def open_help():
    webbrowser.open(
        "https://github.com/akashu4295/Meshless_methods"
    )

def show_about():
    if dpg.does_item_exist("about_window"):
        dpg.configure_item("about_window", show=True)
        dpg.focus_item("about_window")
        return
    with dpg.window(label="About", tag="about_window", modal=True, width=620, height=400, no_resize=True):
        dpg.add_text("MeMPhyS v2.2", color=(200, 220, 255))
        dpg.add_text("Meshless Multi-Physics Solver", color=(200, 220, 255))
        dpg.add_separator()
        dpg.add_text("GUI Development:", color=(255, 220, 160))
        dpg.add_text("\tAkash Unnikrishnan")
        dpg.add_spacer(height=6)
        dpg.add_text("Backend / Simulation:", color=(255, 220, 160))
        dpg.add_text("\tAkash Unnikrishnan")
        dpg.add_text("\tProf. Surya Pratap Vanka")
        dpg.add_text("\tProf. Vinod Narayanan")
        dpg.add_spacer(height=6)
        dpg.add_text("Affiliation:", color=(255, 220, 160))
        dpg.add_text("\tIndian Institute of Technology Gandhinagar")
        dpg.add_text("\tFaculty of Physics, University of Warsaw")
        dpg.add_text("\tUniversity of Illinois Urbana-Champaign")
        dpg.add_spacer(height=12)
        dpg.add_button(label="Close", width=80, callback=lambda: dpg.configure_item("about_window", show=False))

def open_logs():
    if not os.path.exists(LOG_DIR):
        os.makedirs(LOG_DIR)
    open_folder(LOG_DIR)

def save_logs():
    if not os.path.exists(LOG_DIR):
        os.makedirs(LOG_DIR)
    open_folder(LOG_DIR)

def config_dir_selected(sender, app_data):
    selected_dir = app_data["file_path_name"]
    print("Selected config directory:", selected_dir)
    # store this in STATE, load config, etc.

def open_file_dialog(sender, app_data, user_data):
    """Callback for each 'Browse' button."""
    dialog_tag = f"file_dialog_{user_data}"
    if dpg.does_item_exist(dialog_tag):
        dpg.configure_item(dialog_tag, show=True)
    else:
        append_log(f"File dialog {dialog_tag} not found!")
        append_log_file(file, f"File dialog {dialog_tag} not found!")

def show_implicit_callback():
    solver = dpg.get_value("solver_method")
    show_implicit = solver == "Time Implicit"
    for pname in IMPLICIT_PARAMETERS.keys():
        dpg.configure_item(f"param_{pname}", show=show_implicit)

# ============================================================
# Utility / plotting helpers

def update_plot():
    vtk_path = dpg.get_value("contour_vtk_path") or DEFAULT_VTK
    var_choice = dpg.get_value("contour_var")
    cmap = dpg.get_value("contour_cmap") or "viridis"
    dimension = dpg.get_value("param_domain_dimensions")

    script_path = os.path.join("src", "plotter.py")
    proc = STATE.get("plotter_process")
    if proc is not None and proc.poll() is None:
        proc.terminate()
        STATE["plotter_process"] = None

    try:
        proc = subprocess.Popen(
            [
                sys.executable,
                script_path,
                vtk_path,
                var_choice,
                cmap,
                dimension,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        STATE["plotter_process"] = proc
    except Exception as e:
        append_log(f"Plotting Error: {e}")
        append_log_file(file, f"Plotting Error: {e}")

def do_save_image():
    plotter = STATE["plotter"]
    if plotter is None:
        append_log("Save failed: No active plot window found.")
        append_log_file(file, "Save failed: No active plot window found.")
        return
    filename = dpg.get_value("contour_save_path") or "contour.png"
    if not any(filename.lower().endswith(ext) for ext in [".png", ".jpg", ".tif", ".pdf"]):
        filename = filename + ".png"
    try:
        plotter.screenshot(filename, scale=3) 
        append_log(f"Saved high-res image: {filename}")
        append_log_file(file, f"Saved high-res image: {filename}")
    except Exception as e:
        append_log(f"Save failed: {e}")
        append_log_file(file, f"Save failed: {e}")

def update_convergence_plot():
    while dpg.is_dearpygui_running():
        if os.path.exists("Convergence.csv"):
            try:
                df = pd.read_csv("Convergence.csv")
                if df.shape[1] >= 2:
                    x = df.iloc[:, 0].values
                    y = df.iloc[:, 1].values
                    dpg.set_value("conv_series", [x.tolist(), y.tolist()])
                    if len(x) > 0 and len(y) > 0:
                        # Get current y limits to adjust smoothly
                        y_min = np.nanmin(y[y > 0]) if np.any(y > 0) else 1e-8
                        y_max = np.nanmax(y)
                        # Add one order of magnitude padding in log scale
                        y_min = 10 ** (np.floor(np.log10(y_min)) - 1)
                        y_max = 10 ** (np.ceil(np.log10(y_max)) + 0.2)
                        # X-axis: limit to last 500 timesteps or full range
                        x_min = 0
                        x_max = np.max(x) + 10
                        dpg.set_axis_limits("x_axis_conv", x_min, x_max)
                        dpg.set_axis_limits("y_axis_conv", y_min, y_max)
            except Exception as e:
                print("Error updating plot:", e)
        time.sleep(2)

# ============================================================
# GUI Construction
dpg.create_context()

# Define a custom button theme
with dpg.theme() as button_theme:
    with dpg.theme_component(dpg.mvButton):
        dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))          # normal
        dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))   # hover
        dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))    # pressed
        dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)               # rounded corners
        dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)             # inner padding
        dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)

# Define a custom button theme
with dpg.theme() as button_theme2:
    with dpg.theme_component(dpg.mvButton):
        dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))          # normal
        dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))   # hover
        dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))    # pressed
        dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)               # rounded corners
        dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)             # inner padding
        dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)

with dpg.theme() as menu_theme:
    with dpg.theme_component(dpg.mvMenuBar):
        dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg, (28, 28, 32))
        dpg.add_theme_color(dpg.mvThemeCol_Text, (220, 220, 220))
        dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 10, 6)


# Main Window
with dpg.window(label="MeMPhyS GUI", tag="MainWindow", no_close=True):

    with dpg.group(horizontal=True):

        # Left Column
        with dpg.child_window(width=320, border=True, pos=(0, 30)):
            dpg.add_text("Input & Parameters", color=(200, 220, 255))
            dpg.add_spacer(height=5)

            dpg.add_text("Solver Method")
            dpg.add_combo(items=["Fractional Step", "Time Implicit"],
                          default_value="Fractional Step",
                          tag="solver_method", width=200,
                          callback=lambda s, a, u: show_implicit_callback())

            dpg.add_separator()
            dpg.add_text("Flow Parameters", color=(255, 220, 160))
            for pname, pval in BASE_PARAMETERS.items():
                tag = f"param_{pname}"
                dpg.add_input_text(label=pname, tag=tag, default_value=str(pval), width=100)

            # Implicit params (hidden unless selected)
            for pname, pval in IMPLICIT_PARAMETERS.items():
                dpg.add_input_text(label=pname, tag=f"param_{pname}", default_value=str(pval), show=False, width=100)

            # Multigrid section
            dpg.add_separator()
            dpg.add_checkbox(label="Enable Multigrid", tag="multigrid_toggle", callback=lambda s, a, u: show_multigrid_callback())
            
            with dpg.group(horizontal=False, tag="multigrid_parameters_section", show=False):
                dpg.add_text("Multigrid Parameters", color=(200, 255, 200))
                # multigrid parameter inputs
                for pname, pval in MULTIGRID_PARAMETERS.items():
                    dpg.add_input_text(label=pname, tag=f"param_{pname}", default_value=str(pval), width=100)
                # number of mesh levels (create once)
                dpg.add_input_int(label="Number of mesh levels", default_value=1, min_value=1, max_value=10, width=100,
                  tag="num_mesh_levels", callback=lambda s, a, u: update_mesh_inputs())

            # Mesh inputs         
            with dpg.group(tag="multigrid_mesh_section"):
                dpg.add_text("Choose Mesh Files", color=(200, 255, 200))
                dpg.add_text("(finest to coarsest if multigrid)", color=(200, 255, 200))

                for i in range(10):
                    idx = i + 1
                    show_initial = (i == 0)

                    with dpg.group(horizontal=True, show=show_initial, tag=f"mesh_group_{idx}"):
                        dpg.add_input_text(tag=f"mesh_file_{idx}", hint=f"Mesh file {idx} (choose or type path)",width=200,show=show_initial)

                        # Pass index as user_data
                        dpg.add_button(label="Browse", tag=f"browse_{idx}", callback=open_file_dialog, user_data=idx, show=show_initial)

                # File dialogs (defined *after* so they exist when we call show_item)
                for i in range(10):
                    dialog_tag = f"file_dialog_{i+1}"
                    input_tag = f"mesh_file_{i+1}"

                    dpg.add_file_dialog(directory_selector=False, tag=dialog_tag, user_data=input_tag, callback=select_mesh_file, show=False, width=600, height=400)
                    dpg.add_file_extension(".msh", parent=dialog_tag)
            
            dpg.add_separator()               
            dpg.add_text("Choose the initialisation file", color=(255, 220, 160))     
            with dpg.group(horizontal=True, show=True, tag="init_group"):
                dpg.add_input_text(hint="Path to Initialisation file", tag="init_path",width = 200, show = True)
                dpg.add_button(label="Browse", tag="init_browse", callback=open_file_dialog, user_data="init", show=True)
                dpg.add_file_dialog(directory_selector=False, tag="file_dialog_init", user_data="init_path", callback=select_mesh_file, show=False, width=600, height=400)
                dpg.add_file_extension(".c", parent="file_dialog_init")
                        
            # Write button and spacer
            dpg.add_spacer(height=10)
            compile_btn = dpg.add_button(label="Compile and Run Solver", width=250, height=45, callback=run_solver)
            with dpg.tooltip(compile_btn):
                dpg.add_text("Click to compile and execute the solver", color=(220, 230, 255))
                dpg.add_separator()
                # dpg.add_text("This recompiles all .c files and runs the solver executable.")
            dpg.bind_item_theme(compile_btn, button_theme)

            
        # Right Column
        with dpg.child_window(width=-1, border=True):
            dpg.add_text("Convergence Plot", color=(200, 255, 200))
            with dpg.plot(label="Residual", height=300, width=-1):
                dpg.add_plot_axis(dpg.mvXAxis, label="Timestep", tag="x_axis_conv")
                dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", tag="y_axis_conv", log_scale=True)
                dpg.add_line_series([], [], parent="y_axis_conv", tag="conv_series", label="Convergence")
            dpg.add_separator()
            dpg.add_spacer(height=6)

            dpg.add_text("Plotting options", color=(200, 255, 200))
            with dpg.group(horizontal=True):
                dpg.add_text("VTK File:")
                dpg.add_input_text(hint="VTK File", default_value=DEFAULT_VTK, tag="contour_vtk_path", width=150, show = True)
                dpg.add_button(label="Browse", tag="vtk_browse", callback=open_file_dialog, user_data="vtk", show=True)
                dpg.add_file_dialog(directory_selector=False, tag="file_dialog_vtk", user_data="contour_vtk_path", callback=select_mesh_file, show=False, width=600, height=400)
                dpg.add_file_extension(".vtk", parent="file_dialog_vtk")
                dpg.add_text("Variable:")
                dpg.add_combo(("u","v","w","velocity magnitude","p"), label="##Variable", default_value="velocity magnitude", tag="contour_var", width=150)
                dpg.add_text("Colormap:")
                dpg.add_combo(CMAPS, label="##Colormap", default_value="viridis", tag="contour_cmap", width=150)
                plot_btn = dpg.add_button(label="Plot", callback=lambda s,a,u: update_plot())
                dpg.bind_item_theme(plot_btn, button_theme2)
                
            with dpg.group(horizontal=True):
                dpg.add_text("Save Image as:")
                dpg.add_input_text(hint="Save Path", default_value="contour.png", tag="contour_save_path", width=200)
                save_btn = dpg.add_button(label="Save", callback=lambda s,a,u: do_save_image())
                dpg.bind_item_theme(save_btn, button_theme2)

            dpg.add_separator()
            dpg.add_spacer(height=6)
            dpg.add_text("Logs", color=(200, 255, 200))
            
            # with dpg.child_window(tag="log_child", autosize_x=True, height=275, horizontal_scrollbar=True):
            dpg.add_input_text(tag="log_window", multiline=True, readonly=True, width=-1, height = 250,
                 default_value="Ready.\n")
    
            clear_log_btn = dpg.add_button(label="Clear Logs", callback=clear_logs)
            dpg.bind_item_theme(clear_log_btn, button_theme2)

# Viewport setup
dpg.create_viewport(title="MeMPhyS GUI", width=1280, height=800, resizable=True)

with dpg.viewport_menu_bar() as menu_bar:
    with dpg.menu(label="File"):
        dpg.add_menu_item(label="Open Configuration", callback=open_file_dialog, user_data="config")
        dpg.add_file_dialog(directory_selector=False, tag="file_dialog_config", user_data="config_path", callback=select_mesh_file, show=False, width=600, height=400)
        dpg.add_file_extension(".c", parent="file_dialog_config", color=(255, 255, 255, 255))
        dpg.add_menu_item(label="Open Log", callback=open_logs)
        dpg.add_separator()
        dpg.add_menu_item(label="Save Configuration", callback=open_file_dialog, user_data="config_save")
        dpg.add_file_dialog(directory_selector=False, tag="file_dialog_config_save", user_data="config_save_path", callback=select_mesh_file, show=False, width=600, height=400)
        dpg.add_file_extension(".c", parent="file_dialog_config_save", color=(255, 255, 255, 255))
        dpg.add_separator()
        dpg.add_menu_item(label="Exit", callback=dpg.stop_dearpygui)
    # Spacer pushes Help menu to the right
    dpg.add_spacer()
    with dpg.menu(label="Help"):
        dpg.add_menu_item(label="Help", callback=open_help)
        dpg.add_menu_item(label="About", callback=show_about)

dpg.bind_item_theme(menu_bar, menu_theme)

dpg.set_primary_window("MainWindow", True)
dpg.setup_dearpygui()
dpg.show_viewport()
threading.Thread(target=update_convergence_plot, daemon=True).start()
dpg.start_dearpygui()
dpg.destroy_context()