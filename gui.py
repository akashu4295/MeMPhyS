import os
import time
import threading
import numpy as np
import pandas as pd
import dearpygui.dearpygui as dpg
import json
import csv
import subprocess
import platform
import glob
import pyvista as pv
import matplotlib.pyplot as plt
from PIL import Image
import io

# ---------------------------
# Utility / plotting helpers
# ---------------------------
DEFAULT_VTK = "Solution.vtk"
DEFAULT_WIDTH = 900
DEFAULT_HEIGHT = 600
# Gather available matplotlib colormaps for dropdown
CMAPS = sorted([m for m in plt.colormaps()])

# ---------------------------------------------------------
# ---  Window and UI construction
# ---------------------------------------------------------
def build_contour_ui():
    # file selection row
    # with dpg.group(horizontal=True):
    #     dpg.add_input_text(label="VTK File", default_value=DEFAULT_VTK, tag="contour_vtk_path", width=600)
    #     def browse_cb(sender, app_data):
    #         # open file dialog
    #         dpg.show_item("file_dialog_id")
    #     dpg.add_button(label="Browse...", callback=browse_cb)
    # dpg.add_spacer(height=6)
    with dpg.group(horizontal=True):
        dpg.add_input_text(label="VTK File", default_value=DEFAULT_VTK, tag="contour_vtk_path", width=600, show = True)
        dpg.add_button(label="Browse", tag="vtk_browse", callback=open_file_dialog, user_data="vtk", show=True)
        dpg.add_file_dialog(directory_selector=False, tag="file_dialog_vtk", user_data="contour_vtk_path", callback=select_mesh_file, show=False, width=600, height=400)
        dpg.add_file_extension(".vtk", parent="file_dialog_vtk")
    dpg.add_spacer(height=6)
    
    # options: variable, colormap, levels, filled
    with dpg.group(horizontal=True):
        dpg.add_combo(("u","v","w","velocity magnitude","p"), label="Variable", default_value="velocity magnitude", tag="contour_var", width=200)
        dpg.add_combo(CMAPS, label="Colormap", default_value="viridis", tag="contour_cmap", width=200)
        dpg.add_input_int(label="Contour Levels", default_value=10, tag="contour_levels", width=120)
        dpg.add_checkbox(label="Filled contours", default_value=True, tag="contour_filled")

    # custom image size inputs (so user can control resolution)
    with dpg.group(horizontal=True):
        dpg.add_input_int(label="Image Width (px)", default_value=DEFAULT_WIDTH, tag="contour_img_w", width=120)
        dpg.add_input_int(label="Image Height (px)", default_value=DEFAULT_HEIGHT, tag="contour_img_h", width=120)
        dpg.add_button(label="Plot", callback=lambda s,a,u: do_plot_async())
        dpg.add_button(label="Save", callback=lambda s,a,u: do_save_image())
        dpg.add_input_text(label="Save Path", default_value="contour.png", tag="contour_save_path", width=300)

    dpg.add_separator()
    dpg.add_text("Preview (rendered):")
    # image holder (we will add image dynamically)
    with dpg.child_window(tag="contour_image_holder", width=-1, height=400, autosize_x=True, autosize_y=False):
        # placeholder text
        dpg.add_text("No plot yet. Click Plot to render.", tag="contour_image_placeholder")

    dpg.add_separator()
    dpg.add_text("Status log:")
    with dpg.child_window(tag="contour_status_child", width=-1, height=120, autosize_x=True, autosize_y=False):
        dpg.add_input_text(multiline=True, readonly=True, tag="contour_status", default_value="", height=110)

def open_contour_window():
    """Open or show contour plotting window."""
    if dpg.does_item_exist("contour_window"):
        dpg.show_item("contour_window")
        return

    with dpg.window(label="Contour Plot", tag="contour_window",
                    width=dpg.get_viewport_client_width() - 100,
                    height=dpg.get_viewport_client_height() - 100,
                    pos=(50, 50),
                    no_close=False):
        build_contour_ui()


def read_dataset(vtk_path):
    """Read a vtk file using pyvista. Returns a pyvista object."""
    if not os.path.isfile(vtk_path):
        raise FileNotFoundError(f"VTK file not found: {vtk_path}")
    mesh = pv.read(vtk_path)
    return mesh

def get_scalar_field(mesh, varname):
    """Return scalar/vector array and an appropriate name/key for plotting.
       varname: 'u','v','w','velocity','p'."""
    # common keys: try exact matches first, then lower-case matches
    keys = mesh.point_data.keys()
    # normalize keys for lookup
    lower_map = {k.lower(): k for k in keys}
    var = varname.lower()
    if var in ['u','v','w']:
        # attempt to find velocity vector first
        # if vector exists (e.g., 'velocity' or 'vtkVec'), extract the component
        vec_key = None
        # common vector keys:
        for cand in ['velocity', 'vel', 'u_velocity', 'velocity_field']:
            if cand in lower_map:
                vec_key = lower_map[cand]; break
        # else take first vector-like (ndim==3)
        if vec_key is None:
            for k in keys:
                arr = mesh.point_data[k]
                if arr.ndim == 2 and arr.shape[1] in (2,3):
                    vec_key = k; break
        if vec_key is None:
            raise KeyError("No vector velocity field found in VTK file to extract u/v/w.")
        vec = mesh.point_data[vec_key]
        comp_index = {'u':0,'v':1,'w':2}[var]
        if vec.shape[1] <= comp_index:
            # if w not present, pad with zeros
            if comp_index == 2:
                return np.zeros(mesh.n_points), f"{vec_key}_w"
            else:
                raise KeyError(f"Requested component {var} not in vector field shape {vec.shape}.")
        return vec[:, comp_index], f"{vec_key}_{var}"
    elif var == 'velocity' or var == 'magnitude':
        # compute magnitude from vector
        vec_key = None
        for cand in ['velocity', 'vel', 'u_velocity', 'velocity_field']:
            if cand in lower_map:
                vec_key = lower_map[cand]; break
        if vec_key is None:
            for k in keys:
                arr = mesh.point_data[k]
                if arr.ndim == 2 and arr.shape[1] in (2,3):
                    vec_key = k; break
        if vec_key is None:
            raise KeyError("No vector velocity field found in VTK file to compute magnitude.")
        vec = mesh.point_data[vec_key]
        # if 2D vector, append zero z
        if vec.shape[1] == 2:
            vec = np.hstack([vec, np.zeros((vec.shape[0],1))])
        mag = np.linalg.norm(vec, axis=1)
        return mag, f"{vec_key}_mag"
    else:
        # pressure / scalar fields - try to find matching name
        if var in lower_map:
            key = lower_map[var]
            return mesh.point_data[key], key
        # fallback: if 'p' maybe map to 'pressure'
        if var == 'p' and 'pressure' in lower_map:
            key = lower_map['pressure']; return mesh.point_data[key], key
        # else try to find any scalar present
        for k in keys:
            arr = mesh.point_data[k]
            if arr.ndim == 1:
                # return first scalar if nothing else
                return arr, k
        raise KeyError(f"Scalar field '{varname}' not found in VTK point data.")

def pyvista_render_to_image(mesh, scalar_name, cmap, n_levels, filled=True, view='xy', img_size=(DEFAULT_WIDTH, DEFAULT_HEIGHT)):
    """
    Render mesh with pyvista off-screen and return an RGBA numpy array (H x W x 4 uint8).
    scalar_name: string field name in mesh.point_data or a computed array attached to mesh.
    """
    # clone to avoid modifying original
    grid = mesh.copy(deep=True)

    # attach scalar to mesh if scalar_not already present (scalar_name could be key or temporary)
    if isinstance(scalar_name, tuple):
        # (array, name)
        array, name = scalar_name
        grid.point_data[name] = array
        scalar_key = name
    else:
        scalar_key = scalar_name
    
    # Handle scalar range safely
    arr = grid.point_data[scalar_key]
    if np.allclose(np.nanmax(arr), np.nanmin(arr)):
        arr = arr + 1e-6  # prevent uniform color
        grid.point_data[scalar_key] = arr

    clim = (float(np.nanmin(arr)), float(np.nanmax(arr)))
    
    # create plotter off-screen
    plotter = pv.Plotter(off_screen=True, window_size=img_size)
    plotter.set_background("white")
    
    # Add the mesh
    plotter.add_mesh(grid,scalars=scalar_key, cmap=cmap,clim=clim,show_scalar_bar=True,
        smooth_shading=False,lighting=False)
    
    # choose camera orientation for 2D-like datasets if requested
    if view == 'xy':
        plotter.camera_position = 'xy'
    elif view == 'xz':
        plotter.camera_position = 'xz'
    elif view == 'yz':
        plotter.camera_position = 'yz'
    # Add mesh; for filled contours we show scalar shading, for unfilled we add contour lines
    # For filled: use add_mesh with smooth_shading False to show scalar heatmap
    clim = (np.nanmin(grid.point_data[scalar_key]), np.nanmax(grid.point_data[scalar_key]))
    # add geometry
    plotter.add_mesh(grid, scalars=scalar_key, cmap=cmap, clim=clim, show_scalar_bar=True)
    if not filled:
        # create contour lines and overlay
        try:
            iso = grid.contour(n_levels, scalars=scalar_key)
            plotter.add_mesh(iso, color="k", line_width=1.0, opacity=1.0)
        except Exception:
            # fallback: use contour on range manually
            pass

    # screenshot to array (RGB uint8)
    img = plotter.screenshot(return_img=True)
    plotter.close()
    # img is H x W x 3 (RGB); convert to RGBA by adding alpha channel = 255
    if img.ndim == 3 and img.shape[2] == 3:
        a = np.ones((img.shape[0], img.shape[1], 1), dtype=np.uint8) * 255
        img = np.concatenate([img, a], axis=2)
    return img


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
    append_log("Parameters written to flow_parameters.csv\n")

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

    with open("grid_filenames.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["num_levels", num_levels])
        for path in mesh_paths:
            writer.writerow([path])

    append_log("Grid details written to grid_filesnames.csv\n")

# =============================================================
# Function to compile and run the C solver, streaming logs to GUI

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
        
def run_solver(sender):
    """Compile and run the C solver, streaming logs into the GUI."""
    # disable button while running
    write_grid_csv()
    write_params_csv()
    dpg.disable_item(sender)
    dpg.configure_item(sender, label="Compiling and Running...")
    
    dpg.set_value("log_window", "Starting new run...\n")
    append_log("Saving Meshfile details and parameters...\n")
    write_grid_csv()
    # Detect OS
    system_type = platform.system()
    compiler_cmd = []
    run_cmd = []

    # Path to the init .c file
    init_path = dpg.get_value("init_path")    #"init/init_TC.c"

    # Collect all .c files in the 'header_files' directory
    header_dir = "header_files"
    header_c_files = glob.glob(os.path.join(header_dir, "*.c"))

    # Build the compile command
    if system_type == "Windows":
        compiler_cmd = ["gcc"] + header_c_files + [init_path, "mg_NS_solver.c", "-lm", "-o", "solver.exe"]
        run_cmd = ["solver.exe"]
    else:
        compiler_cmd = ["gcc"] + header_c_files + [init_path, "mg_NS_solver.c", "-lm", "-o", "solver"]
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
            append_log("Solver completed successfully.\n Solution file saved as Solution.csv")
        else:
            append_log(f"Solver exited with code {process.returncode}.")
        
        dpg.configure_item(sender, label="Done! Run Again")
        dpg.enable_item(sender)

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

##### Plotting callback functions
# Texture registry and global state for current image
CURRENT_IMAGE_TAG = "contour_img_texture"
CURRENT_IMAGE_WIDTH_TAG = "contour_img_w"
CURRENT_IMAGE_HEIGHT_TAG = "contour_img_h"
STATE = {"last_image_array": None}

def numpy_to_dpg_texture(arr, tag="CURRENT_IMAGE_TAG"):
    """Convert an HxWx4 uint8 array to float 0-1 list and register as DPG static texture."""
    h, w, c = arr.shape
    assert c == 4, "Array must be RGBA (HxWx4)"
    
    data = (arr.astype(np.float32) / 255.0).flatten().tolist()

    # Ensure a texture registry exists
    if not dpg.does_item_exist("texture_registry"):
        with dpg.texture_registry(tag="texture_registry"):
            pass

    # Delete existing texture if present
    if dpg.does_item_exist(tag):
        dpg.delete_item(tag)

    # Always create textures inside the registry
    with dpg.texture_registry():
        dpg.add_static_texture(w, h, data, tag=tag)

    # # Optionally store metadata
    dpg.set_value("contour_img_w", w)
    dpg.set_value("contour_img_h", h)


def append_status(msg):
    """Append to status log (simple multiline value)."""
    cur = dpg.get_value("contour_status") or ""
    dpg.set_value("contour_status", cur + msg + "\n")
    dpg.split_frame()
    # autoscroll the child
    try:
        dpg.set_y_scroll("contour_status_child", dpg.get_y_scroll_max("contour_status_child"))
    except Exception:
        pass

def do_plot_async(sender=None):
    """Collect UI options and run plotting in a background thread."""
    def worker():
        try:
            vtk_path = dpg.get_value("contour_vtk_path") or DEFAULT_VTK
            append_status(f"Reading: {vtk_path}")
            mesh = read_dataset(vtk_path)

            var_choice = dpg.get_value("contour_var")
            if var_choice == "velocity magnitude":
                arr, name = get_scalar_field(mesh, "velocity")
                scalar = (arr, name)
            else:
                arr, name = get_scalar_field(mesh, var_choice)
                scalar = (arr, name)

            cmap = dpg.get_value("contour_cmap") or "viridis"
            n_levels = int(dpg.get_value("contour_levels") or 10)
            filled = dpg.get_value("contour_filled")

            append_status(f"Preparing plot: var={var_choice}, cmap={cmap}, levels={n_levels}, filled={filled}")
            img_size = (dpg.get_value("contour_img_w") or DEFAULT_WIDTH, dpg.get_value("contour_img_h") or DEFAULT_HEIGHT)

            # Render to image
            img = pyvista_render_to_image(mesh, scalar, cmap=cmap, n_levels=n_levels, filled=filled, img_size=img_size)
            STATE["last_image_array"] = img

            # register texture (must be done on main thread)
            def finish():
                numpy_to_dpg_texture(img, tag=CURRENT_IMAGE_TAG)
                # create or update image widget
                if not dpg.does_item_exist("contour_image"):
                    dpg.add_image(CURRENT_IMAGE_TAG, tag="contour_image", parent="contour_image_holder")
                else:
                    dpg.configure_item("contour_image", texture_tag=CURRENT_IMAGE_TAG)
                append_status("Plot ready.")
            finish()

        except Exception as e:
            lambda: append_status(f"Error: {e}")

    threading.Thread(target=worker, daemon=True).start()

def pyvista_render_to_image2(mesh, scalar_name, cmap, n_levels, filled=True,
                            view='xy', img_size=(800, 600)):
    """Render mesh with pyvista off-screen and return RGBA image (HxWx4)."""
    grid = mesh.copy(deep=True)

    # Handle scalar_name being tuple
    if isinstance(scalar_name, tuple):
        array, name = scalar_name
        grid.point_data[name] = array
        scalar_key = name
    else:
        scalar_key = scalar_name

    # Handle scalar range safely
    arr = grid.point_data[scalar_key]
    if np.allclose(np.nanmax(arr), np.nanmin(arr)):
        arr = arr + 1e-6  # prevent uniform color
        grid.point_data[scalar_key] = arr

    clim = (float(np.nanmin(arr)), float(np.nanmax(arr)))

    plotter = pv.Plotter(off_screen=True, window_size=img_size)
    plotter.set_background("white")

    # Add the mesh
    plotter.add_mesh(
        grid,
        scalars=scalar_key,
        cmap=cmap,
        clim=clim,
        show_scalar_bar=True,
        smooth_shading=False,
        lighting=False,
    )

    # Overlay contour lines if requested
    if not filled:
        try:
            iso = grid.contour(n_levels, scalars=scalar_key)
            plotter.add_mesh(iso, color="black", line_width=1.0)
        except Exception as e:
            print("[Warning] Contour lines failed:", e)

    # Reset camera to fit mesh in view
    plotter.camera_position = view
    # plotter.show_bounds(False)
    plotter.view_xy()
    plotter.reset_camera()

    # Capture image
    img = plotter.screenshot(return_img=True)
    plotter.close()

    # Convert to RGBA
    if img.ndim == 3 and img.shape[2] == 3:
        alpha = np.full((*img.shape[:2], 1), 255, dtype=np.uint8)
        img = np.concatenate([img, alpha], axis=2)
    return img


def do_save_image(sender):
    """Save current image using PIL with high DPI metadata."""
    if STATE.get("last_image_array") is None:
        append_status("No image to save.")
        return
    # ask user for filename via file dialog
    filename = dpg.get_value("contour_save_path") or "contour.png"
    # ensure extension
    if not any(filename.lower().endswith(ext) for ext in [".png", ".jpg", ".tif", ".tiff"]):
        filename = filename + ".png"
    try:
        arr = STATE["last_image_array"]
        img = Image.fromarray(arr)
        # Save with high DPI metadata (600)
        img.save(filename, dpi=(600,600))
        append_status(f"Saved image: {filename} (dpi metadata=600)")
    except Exception as e:
        append_status(f"Save failed: {e}")


# ============================================================
# GUI Construction

dpg.create_context()

# ----------------------------
# Define a custom button theme
# ----------------------------
with dpg.theme() as button_theme:
    with dpg.theme_component(dpg.mvButton):
        dpg.add_theme_color(dpg.mvThemeCol_Button, (40, 120, 200))          # normal
        dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (60, 140, 230))   # hover
        dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (30, 100, 180))    # pressed
        dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)               # rounded corners
        dpg.add_theme_style(dpg.mvStyleVar_FramePadding, 12, 6)             # inner padding
        dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing, 10, 10)

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
                        dpg.add_input_text(tag=f"mesh_file_{idx}", hint=f"Mesh file {idx} (choose or type path)",width=300,show=show_initial)

                        # Pass index as user_data
                        dpg.add_button(label="Browse", tag=f"browse_{idx}", callback=open_file_dialog, user_data=idx, show=show_initial)

                # --- File dialogs (defined *after* so they exist when we call show_item) ---
                for i in range(10):
                    dialog_tag = f"file_dialog_{i+1}"
                    input_tag = f"mesh_file_{i+1}"

                    dpg.add_file_dialog(directory_selector=False, tag=dialog_tag, user_data=input_tag, callback=select_mesh_file, show=False, width=600, height=400)
                    dpg.add_file_extension(".msh", parent=dialog_tag)
            
            dpg.add_separator()               
            dpg.add_text("Choose the initialisation file", color=(255, 220, 160))     
            with dpg.group(horizontal=True, show=True, tag="init_group"):
                dpg.add_input_text(hint="Path to Initialisation file", tag="init_path",width = 300, show = True)
                dpg.add_button(label="Browse", tag="init_browse", callback=open_file_dialog, user_data="init", show=True)
                dpg.add_file_dialog(directory_selector=False, tag="file_dialog_init", user_data="init_path", callback=select_mesh_file, show=False, width=600, height=400)
                dpg.add_file_extension(".c", parent="file_dialog_init")
                        
            # Write button and spacer
            # dpg.add_spacer(height=10)
            # dpg.add_button(label="Save Parameters and mesh details", callback=on_write_callback)
            dpg.add_spacer(height=10)
            # dpg.add_button(label="Compile and Run Solver", callback=run_solver)
            compile_btn = dpg.add_button(label="Compile and Run Solver", width=250, height=45, callback=run_solver)
            with dpg.tooltip(compile_btn):
                dpg.add_text("Click to compile and execute the solver", color=(220, 230, 255))
                dpg.add_separator()
                dpg.add_text("This recompiles all .c files and runs the solver executable.")
            dpg.bind_item_theme(compile_btn, button_theme)
            
            dpg.add_spacer(height=10)
            plot_btn = dpg.add_button(label="Plot Contours", width=250, height=45, callback=open_contour_window)
            with dpg.tooltip(plot_btn):
                dpg.add_text("Click to plot contours of the solution data", color=(220, 230, 255))
                dpg.add_separator()
                dpg.add_text("This opens a new window for all the plotting options")
            dpg.bind_item_theme(plot_btn, button_theme)
            

        # ===== Right Column =====
        with dpg.child_window(width=-1, border=True):
            dpg.add_text("Convergence Plot", color=(200, 255, 200))
            with dpg.plot(label="Residual", height=300, width=-1):
                dpg.add_plot_axis(dpg.mvXAxis, label="Timestep", tag="x_axis_conv")
                dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", tag="y_axis_conv", log_scale=True)
                # dpg.add_plot_axis(dpg.mvYAxis, label="Total Steady-State Error", tag="y_axis_conv", scale=dpg.mvScale_Log10)
                dpg.add_line_series([], [], parent="y_axis_conv", tag="conv_series", label="Convergence")

            # dpg.add_separator()
            # dpg.add_text("Flow Contour (|U| with streamlines)", color=(200, 255, 255))
            # # Placeholder image (light gray texture)
            # w0, h0 = 200, 150
            # empty = np.zeros((h0, w0, 4), dtype=np.float32) + 0.8  # light gray RGBA
            # empty[:, :, 3] = 1.0

            # # Register the texture before using it
            # with dpg.texture_registry(show=False):
            #     dpg.add_static_texture(w0, h0, empty.flatten().tolist(), tag="contour_placeholder")

            # # Now add the image referring to that texture
            # dpg.add_image("contour_placeholder", tag="contour_image")

            dpg.add_separator()
            dpg.add_text("Logs")

            with dpg.group(horizontal=True):  # Put button next to label if you prefer
                dpg.add_button(label="Clear Logs", callback=clear_logs)

            with dpg.child_window(tag="log_child", autosize_x=True, height=350, horizontal_scrollbar=True):
                dpg.add_input_text(
                    tag="log_window",
                    multiline=True,
                    readonly=True,
                    width=-1, height = 330,
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