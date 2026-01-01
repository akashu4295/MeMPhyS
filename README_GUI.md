# MeMPhyS GUI Guide

# [TODO]MeMPhyS GUI Restructuring 

### High Priority
- [o] Test full solver workflow (compile → run → plot)
- [o] Verify mesh file validation
- [ ] Test all callbacks thoroughly
- [ ] Remove debug print statements from main.py

### Medium Priority
- [ ] Enable custom fonts (set `ENABLE_CUSTOM_FONTS = True`)
- [ ] Implement config save/load functionality
- [ ] Add unit tests for modules
- [ ] Add more validation rules

### Low Priority
- [ ] Advanced plot settings dialog
- [ ] Keyboard shortcuts
- [ ] Drag-and-drop mesh files
- [ ] Recent files menu

## Known Limitations

1. **Font Changes**: Require application restart (DearPyGUI limitation)
2. **Screenshot Feature**: External plotter window only
3. **Config Save/Load**: Not yet implemented (placeholders exist)

## New Directory Structure

```
MeMPhyS/
├── main.py                      # New entry point (replaces old gui.py)
├── src/
│   ├── config/                  # Configuration and constants
│   │   ├── __init__.py
│   │   ├── constants.py         # All constants and parameters
│   │   └── themes.py            # DearPyGUI themes
│   │
│   ├── core/                    # Core functionality
│   │   ├── __init__.py
│   │   ├── state.py             # Application state management
│   │   └── logger.py            # Logging system
│   │
│   ├── utils/                   # Utility functions
│   │   ├── __init__.py
│   │   ├── fonts.py             # Font management
│   │   ├── file_io.py           # File operations
│   │   └── platform_utils.py   # OS-specific utilities
│   │
│   ├── solver/                  # Solver management
│   │   ├── __init__.py
│   │   ├── runner.py            # Compilation and execution
│   │   └── monitoring.py        # Convergence monitoring
│   │
│   ├── callbacks/               # Callback functions
│   │   ├── __init__.py
│   │   ├── solver_callbacks.py  # Solver operations
│   │   ├── mesh_callbacks.py    # Mesh management
│   │   ├── plot_callbacks.py    # Visualization
│   │   └── menu_callbacks.py    # Menu actions
│   │
│   ├── ui/                      # User interface
│   │   ├── __init__.py
│   │   ├── main_window.py       # Main window assembly
│   │   ├── menu_bar.py          # Menu bar
│   │   ├── parameters_panel.py  # Left panel
│   │   ├── visualization_panel.py # Right panel
│   │   └── dialogs.py           # Modal dialogs
│   │
│   └── plotting/
│       └── plotter.py           # External plotting (existing)
│
├── header_files/                # C header files (unchanged)
├── logs/                        # Log files
└── requirements.txt             # Python dependencies
```

## Quick Start

### Running the Application

```bash
python gui.py
```

## Module Descriptions

### 1. **src/config/** - Configuration

#### constants.py
Contains all application constants:
- Solver parameters (BASE_PARAMETERS, IMPLICIT_PARAMETERS, MULTIGRID_PARAMETERS)
- GUI dimensions and settings
- File paths and names
- Colors and styling
- Validation rules (PARAMETER_CONSTRAINTS)

#### themes.py
All DearPyGUI themes:
- Button themes (primary, secondary, success, error)
- Input themes
- Plot themes
- Modal themes
- `initialize_all_themes()` - Create all themes
- `apply_global_theme()` - Apply app-wide styling

### 2. **src/core/** - Core Functionality

#### state.py
Application state management via `AppState` class:
- Process management (plotter, solver)
- File handle management
- Font registry
- UI state (multigrid enabled, mesh levels, etc.)
- `app_state` - Global instance

**Key Methods:**
```python
app_state.solver_running = True
app_state.set_mesh_file(1, "mesh.msh")
app_state.cleanup()  # Call on exit
```

#### logger.py
Structured logging system via `Logger` class:
- Multiple outputs (file, GUI, console)
- Log levels (INFO, WARNING, ERROR, SUCCESS, DEBUG)
- Colored console output
- Auto-scrolling GUI log
- `logger` - Global instance

**Key Methods:**
```python
logger.info("Information message")
logger.error("Error message")
logger.success("Success message")
logger.log_exception(e, "Context")
logger.separator()
```

### 3. **src/utils/** - Utilities

#### fonts.py
Font management:
- `find_system_font()` - Cross-platform font discovery
- `initialize_fonts()` - Load all fonts at startup
- `change_font()` - Change current font
- Platform-specific font paths

#### file_io.py
File operations:
- `write_parameters_csv()` - Write parameters to CSV
- `write_grid_csv()` - Write grid configuration
- `validate_mesh_file()` - Validate mesh files
- `read_csv_file()` / `write_csv_file()` - CSV utilities

#### platform_utils.py
OS-specific utilities:
- `open_folder()` - Open in file explorer
- `open_url()` - Open in browser
- `get_platform()` - Get OS name
- `build_compile_command()` - Build compilation command

### 4. **src/solver/** - Solver Management

#### runner.py
Solver compilation and execution via `SolverRunner` class:
- `compile_and_run()` - Main entry point
- `stop_solver()` - Terminate running solver
- `check_dependencies()` - Verify compiler, files exist
- Asynchronous execution with live output
- `solver_runner` - Global instance

**Usage:**
```python
from src.solver import solver_runner

success = solver_runner.compile_and_run(
    init_file="init_cavity.c",
    button_tag="run_button",
    on_complete=my_callback
)
```

#### monitoring.py
Convergence monitoring via `ConvergenceMonitor` class:
- Real-time CSV file monitoring
- Automatic plot updates
- Smart axis scaling
- Convergence status tracking
- `convergence_monitor` - Global instance

**Usage:**
```python
from src.solver import convergence_monitor

convergence_monitor.start()
status = convergence_monitor.get_convergence_status()
convergence_monitor.stop()
```

### 5. **src/callbacks/** - Callbacks

All callback functions organized by function:

#### solver_callbacks.py
- `run_solver_callback()` - Run solver
- `validate_numeric_input()` - Parameter validation
- `show_implicit_callback()` - Show/hide implicit params
- `stop_solver_callback()` - Stop solver

#### mesh_callbacks.py
- `show_multigrid_callback()` - Enable/disable multigrid
- `update_mesh_inputs_callback()` - Update mesh input visibility
- `select_mesh_file_callback()` - File selection
- `validate_all_mesh_files_callback()` - Validate all meshes
- `auto_fill_mesh_levels_callback()` - Auto-generate mesh names

#### plot_callbacks.py
- `update_plot_callback()` - Launch plotter
- `reset_convergence_plot_callback()` - Reset plot
- `export_convergence_data_callback()` - Export data to CSV

#### menu_callbacks.py
- `open_logs_callback()` - Open logs folder
- `show_about_callback()` - Show About dialog
- `show_preferences_callback()` - Show Preferences
- `exit_application_callback()` - Clean exit

### 6. **src/ui/** - User Interface

#### main_window.py
- `create_main_window()` - Assemble main window

#### parameters_panel.py
- `create_parameters_panel()` - Left panel with parameters

#### visualization_panel.py
- `create_visualization_panel()` - Right panel with plots

#### dialogs.py
- `create_about_dialog()` - About dialog
- `create_preferences_dialog()` - Preferences dialog
- `show_error_dialog()` - Error message
- `show_confirmation_dialog()` - Yes/No confirmation


## Configuration

### Enable Custom Fonts
In `main.py`, line ~47:
```python
ENABLE_CUSTOM_FONTS = True  # Change to True
```

### Adjust Convergence Update Interval
In `src/config/constants.py`:
```python
CONVERGENCE_UPDATE_INTERVAL = 2.0  # seconds
```

### Change Log Level
In your code:
```python
logger.set_enable_console(True)   # Console output
logger.set_enable_file(True)      # File output
logger.set_enable_gui(True)       # GUI output
```

## Code Examples

### Using the Logger
```python
from src.core import logger

logger.info("Information message")
logger.success("Success message")
logger.error("Error message")
logger.warning("Warning message")
logger.separator()
```

### Using App State
```python
from src.core import app_state

app_state.solver_running = True
app_state.set_mesh_file(1, "mesh.msh")
app_state.cleanup()  # On exit
```

### Using Solver
```python
from src.solver import solver_runner

solver_runner.compile_and_run(
    init_file="init_cavity.c",
    button_tag="run_button"
)
```

## Troubleshooting

### Convergence plot not updating
- Check if `Convergence.csv` exists
- Verify convergence monitor started (check logs)
- Run solver to generate data

### Fonts not loading
- Check available fonts for your OS
- See `src/utils/fonts.py` for font paths
- Try different font from Preferences

### Solver won't compile
- Verify gcc is installed: `gcc --version`
- Check init file path is correct
- Check mesh files exist
- Review logs for compilation errors

## Getting Help

1. Check logs: `logs/log_YYYY-MM-DD.txt`
2. Check GitHub issues
3. Enable debug logging


---

**Version**: 2.2 Restructured  
**Date**: December 2025  
**Maintainer**: Akash Unnikrishnan
