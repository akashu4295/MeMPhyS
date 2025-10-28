import os
import time
import threading
import numpy as np
import pandas as pd
import dearpygui.dearpygui as dpg
import json
import csv
import sys
import subprocess
import platform
# import glob

# ============================================================
# Default parameters for the GUI

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

PARAM_TAGS = []

# ---------------------- Load parameter info ---------------------- #
PARAM_ORDER_FILE = "param_order.json"
PARAM_DEFAULTS_FILE = "param_defaults.json"

# Load order and defaults
if os.path.exists(PARAM_ORDER_FILE):
    with open(PARAM_ORDER_FILE, "r") as f:
        PARAM_ORDER = json.load(f)
else:
    PARAM_ORDER = []

if os.path.exists(PARAM_DEFAULTS_FILE):
    with open(PARAM_DEFAULTS_FILE, "r") as f:
        PARAM_DEFAULTS = json.load(f)
else:
    PARAM_DEFAULTS = {}
    
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
    dpg.set_value("log_window", f"Parameters written to {filename}\n")

def write_grid_csv():
    """Writes mesh file paths to test_grid_files.csv"""
    mg_enabled = dpg.get_value("multigrid_toggle")
    num_levels = 1
    mesh_paths = []

    # if mg_enabled:
    num_levels = int(dpg.get_value("num_mesh_levels"))
    for i in range(num_levels):
        path = dpg.get_value(f"mesh_file_{i+1}")
        if path:
            mesh_paths.append(path)
    # else:
    #     path = dpg.get_value("single_mesh_file")
    #     if path:
    #         mesh_paths.append(path)

    with open("test_grid_files.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["num_levels", num_levels])
        for path in mesh_paths:
            writer.writerow([path])

    dpg.set_value("log_window", "test_grid_files.csv written successfully.\n")

# =============================================================
# Function to compile and run the C solver, streaming logs to GUI

def append_log(msg):
    current = dpg.get_value("log_window") or ""
    dpg.set_value("log_window", current + msg + "\n")

    # Scroll the parent child_window instead
    try:
        dpg.set_y_scroll("log_child", 999999)
    except Exception:
        pass  # just in case it's not created yet

def clear_logs():
    dpg.set_value("log_window", "")
    dpg.set_y_scroll("log_child", 0)
        
def run_solver():
    """Compile and run the C solver, streaming logs into the GUI."""
    dpg.set_value("log_window", "Starting new run...\n")
    
    # Detect OS
    system_type = platform.system()
    compiler_cmd = []
    run_cmd = []

    # --- Build compile command based on OS ---
    if system_type == "Windows":
        compiler_cmd = ["gcc", "mg_NS_solver.c", "-lm", "-o", "solver.exe"]
        run_cmd = ["solver.exe"]
    else:
        compiler_cmd = ["gcc", "mg_NS_solver.c", "-lm", "-o", "solver"]
        run_cmd = ["./solver"]

    append_log("Compiling solver...")

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
        return
    else:
        append_log("Compilation successful. Running solver...\n")

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
            append_log("Solver completed successfully.")
        else:
            append_log(f"Solver exited with code {process.returncode}.")

    threading.Thread(target=solver_thread, daemon=True).start()


# ============================================================
# Callbacks : functions triggered by GUI events

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

def browse_directory_callback():
    num_levels = int(dpg.get_value("num_mesh_levels"))
    for i in range(10):
        if dpg.does_item_exist(f"file_dialog_{i+1}"):
            dpg.configure_item(f"file_dialog_{i+1}", show=(i+1 <= num_levels))

def show_implicit_callback():
    solver = dpg.get_value("solver_method")
    show_implicit = solver == "Implicit"
    for pname in IMPLICIT_PARAMETERS.keys():
        dpg.configure_item(f"param_{pname}", show=show_implicit)

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
        
def open_file_dialog(sender, app_data, user_data):
    """Callback for each 'Browse' button."""
    dialog_tag = f"file_dialog_{user_data}"
    if dpg.does_item_exist(dialog_tag):
        dpg.configure_item(dialog_tag, show=True)
    else:
        append_log(f"File dialog {dialog_tag} not found!")

# ============================================================
# GUI Construction

dpg.create_context()

# ---------- Main Window ----------
with dpg.window(label="MeMPhyS GUI", tag="MainWindow", no_close=True):
    with dpg.group(horizontal=True):

        # ===== Left Column =====
        with dpg.child_window(width=420, border=True):
            dpg.add_text("Input & Parameters", color=(200, 220, 255))
            dpg.add_spacer(height=5)

            dpg.add_text("Solver Method")
            dpg.add_combo(items=["Fractional Step", "Implicit"],
                          default_value="Fractional Step",
                          tag="solver_method",
                          callback=lambda s, a, u: show_implicit_callback())

            dpg.add_separator()
            dpg.add_text("Flow Parameters", color=(255, 220, 160))
            for pname, pval in BASE_PARAMETERS.items():
                tag = f"param_{pname}"
                PARAM_TAGS.append(tag)
                dpg.add_input_text(label=pname, tag=tag, default_value=str(pval))

            # Implicit params (hidden unless selected)
            for pname, pval in IMPLICIT_PARAMETERS.items():
                dpg.add_input_text(label=pname, tag=f"param_{pname}", default_value=str(pval), show=False)

            # Multigrid section
            dpg.add_separator()
            dpg.add_checkbox(label="Enable Multigrid", tag="multigrid_toggle", callback=lambda s, a, u: show_multigrid_callback())
            
            with dpg.group(horizontal=False, tag="multigrid_parameters_section", show=False):
                dpg.add_text("Multigrid Parameters", color=(200, 255, 200))
                # multigrid parameter inputs
                for pname, pval in MULTIGRID_PARAMETERS.items():
                    dpg.add_input_text(label=pname, tag=f"param_{pname}", default_value=str(pval))
                # number of mesh levels (create once)
                dpg.add_input_int(label="Number of mesh levels", default_value=1, min_value=1, max_value=10,
                  tag="num_mesh_levels", callback=lambda s, a, u: update_mesh_inputs())

            # ---- Mesh inputs ---- #            
            with dpg.group(tag="multigrid_mesh_section"):
                dpg.add_text("Choose Mesh Files (finest to coarsest if multigrid):", color=(200, 255, 255))

                for i in range(10):
                    idx = i + 1
                    show_initial = (i == 0)

                    with dpg.group(horizontal=True, show=show_initial, tag=f"mesh_group_{idx}"):
                        dpg.add_input_text(
                            tag=f"mesh_file_{idx}",
                            hint=f"Mesh file {idx} (choose or type path)",
                            width=300,
                            show=show_initial
                        )

                        # Pass index as user_data
                        dpg.add_button(
                            label="Browse",
                            tag=f"browse_{idx}",
                            callback=open_file_dialog,
                            user_data=idx,
                            show=show_initial
                        )

                # --- File dialogs (defined *after* so they exist when we call show_item) ---
                for i in range(10):
                    dpg.add_file_dialog(
                        directory_selector=False,
                        tag=f"file_dialog_{i+1}",
                        callback=select_mesh_file,
                        user_data=f"mesh_file_{i+1}",
                        show=False,
                        width=600,
                        height=400
                    )

            # Write button and spacer
            dpg.add_spacer(height=30)
            dpg.add_button(label="Save Parameters and mesh details", callback=on_write_callback)
            dpg.add_spacer(height=30)
            dpg.add_button(label="Compile and Run Solver", callback=run_solver)

        # ===== Right Column =====
        with dpg.child_window(width=-1, border=True):
            dpg.add_text("Convergence Plot", color=(200, 255, 200))
            with dpg.plot(label="Residual", height=300, width=-1):
                dpg.add_plot_axis(dpg.mvXAxis, label="Timestep", tag="x_axis_conv")
                dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", tag="y_axis_conv", log_scale=True)
                # dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", tag="y_axis_conv", scale=dpg.mvScale_Log10)
                dpg.add_line_series([], [], parent="y_axis_conv", tag="conv_series", label="Convergence")

            # dpg.add_text("Convergence Plot", color=(200, 255, 200))
            # with dpg.plot(label="Residual", height=300, width=-1):
            #     dpg.add_plot_axis(dpg.mvXAxis, label="Timestep")
            #     y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", log_scale=True)
            #     dpg.add_line_series([], [], parent=y_axis, tag="conv_series", label="Convergence")

            dpg.add_separator()
            dpg.add_text("Flow Contour (|U| with streamlines)", color=(200, 255, 255))
            # Placeholder image (light gray texture)
            w0, h0 = 200, 150
            empty = np.zeros((h0, w0, 4), dtype=np.float32) + 0.8  # light gray RGBA
            empty[:, :, 3] = 1.0

            # Register the texture before using it
            with dpg.texture_registry(show=False):
                dpg.add_static_texture(w0, h0, empty.flatten().tolist(), tag="contour_placeholder")

            # Now add the image referring to that texture
            dpg.add_image("contour_placeholder", tag="contour_image")

            dpg.add_separator()
            dpg.add_text("Logs")
            # with dpg.child_window(tag="log_child", autosize_x=True, height=150, horizontal_scrollbar=True):
            #     dpg.add_input_text(label="Logs", tag="log_window", multiline=True, readonly=True, height=150, width=-1, default_value="Ready.\n")

            with dpg.group(horizontal=True):  # Put button next to label if you prefer
                dpg.add_button(label="Clear Logs", callback=clear_logs)

            with dpg.child_window(tag="log_child", autosize_x=True, height=150, horizontal_scrollbar=True):
                dpg.add_input_text(
                    tag="log_window",
                    multiline=True,
                    readonly=True,
                    width=-1,
                    default_value="Ready.\n"
                )

# ============================================================
# Viewport setup
dpg.create_viewport(title="MeMPhyS GUI", width=1280, height=800, resizable=True)
dpg.set_primary_window("MainWindow", True)
dpg.setup_dearpygui()
dpg.show_viewport()

# Background thread for live plotting
threading.Thread(target=update_convergence_plot, daemon=True).start()

# Start DearPyGui
dpg.start_dearpygui()
dpg.destroy_context()