# MeMPhyS GUI Guide

## Overview

The MeMPhyS GUI has been completely restructured from a single monolithic file into a well-organized, modular architecture. This guide explains the new structure, how to use it, and the migration process.

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
# Old way (deprecated)
python gui.py

# New way
python main.py
```

### Basic Usage Example

```python
from src.core import app_state, logger
from src.solver import solver_runner
from src.config import BASE_PARAMETERS

# Log something
logger.info("Starting simulation")

# Check solver status
if solver_runner.is_solver_running():
    logger.info("Solver is running")

# Access parameters
reynolds = BASE_PARAMETERS["Re"]
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

## Migration from Old Code

### Old Code (Single File)
```python
# Everything in one file
STATE = {"plotter": None}
BASE_PARAMETERS = {...}

def run_solver(sender):
    # Implementation

dpg.create_context()
# 700+ lines of UI construction
dpg.start_dearpygui()
```

### New Code (Modular)
```python
# main.py - Clean entry point
from src.config import BASE_PARAMETERS
from src.core import app_state, logger
from src.ui import create_main_window

def main():
    themes = initialize_application()
    create_gui(themes)
    dpg.start_dearpygui()
```

## Benefits of New Structure

1. **Maintainability**: Easy to find and modify code
2. **Testability**: Can unit test individual modules
3. **Scalability**: Easy to add features
4. **Readability**: Clear separation of concerns
5. **Reusability**: Utility functions can be reused
6. **Collaboration**: Multiple developers can work simultaneously
7. **Debugging**: Easier to track down issues
8. **Documentation**: Self-documenting structure

## Bugs Fixed in Restructuring

1. Convergence monitor stops cleanly on exit
2. Log file properly closed
3. Processes terminated safely
4. Thread-safe GUI updates
5. Input validation for parameters
6. Better error handling throughout
7. No more resource leaks

## Adding New Features

### Adding a New Parameter

1. Add to `src/config/constants.py`:
```python
BASE_PARAMETERS = {
    ...
    "new_parameter": 1.0,
}
```

2. Add constraints (optional):
```python
PARAMETER_CONSTRAINTS = {
    ...
    "new_parameter": (0.0, 10.0, float),
}
```

The parameter will automatically appear in the GUI!

### Adding a New Callback

1. Create callback in appropriate file (e.g., `src/callbacks/solver_callbacks.py`):
```python
def my_new_callback(sender, app_data, user_data):
    logger.info("New callback triggered")
    # Implementation
```

2. Export in `src/callbacks/__init__.py`:
```python
from .solver_callbacks import my_new_callback

__all__ = [..., 'my_new_callback']
```

3. Use in UI:
```python
dpg.add_button(label="New Feature", callback=my_new_callback)
```

### Adding a New UI Panel

1. Create file `src/ui/my_panel.py`:
```python
def create_my_panel(themes: dict) -> int:
    with dpg.child_window(...) as panel:
        # Panel contents
    return panel
```

2. Import and use in `main_window.py`

## Testing

### Testing Individual Modules

```python
# Test logger
from src.core import logger
logger.info("Test message")
logger.error("Test error")

# Test file operations
from src.utils import validate_mesh_file
is_valid, msg = validate_mesh_file("mesh.msh")

# Test solver
from src.solver import solver_runner
deps_ok, missing = solver_runner.check_dependencies()
```

## Import Reference

```python
# Configuration
from src.config import BASE_PARAMETERS, COLORS
from src.config.themes import initialize_all_themes

# Core
from src.core import app_state, logger

# Utils
from src.utils import (
    initialize_fonts,
    write_parameters_csv,
    open_folder
)

# Solver
from src.solver import solver_runner, convergence_monitor

# Callbacks
from src.callbacks import (
    run_solver_callback,
    show_multigrid_callback
)

# UI
from src.ui import create_main_window
```

## Troubleshooting

### Issue: Fonts not loading
- Check `src/config/constants.py` - `FONT_PREFERENCES`
- Verify fonts exist on your system
- Logs will show which fonts were found

### Issue: Solver won't compile
- Check `solver_runner.check_dependencies()`
- Verify gcc is in PATH
- Check init file path is correct

### Issue: Convergence plot not updating
- Check `Convergence.csv` exists
- Verify convergence monitor is running: `convergence_monitor.is_running()`
- Check logs for errors
 
## Next Steps

1. **Code is restructured**
2. **Test thoroughly** - Run through all features
3. **Add unit tests** - Test individual modules
4. **Update documentation** - User manual
5. **Add more features** - Config save/load, advanced validation

## Support

For issues or questions:
- Check logs in `logs/` directory
- Enable debug logging: `logger.set_enable_console(True)`
- Review this guide
- Check GitHub repository

---

**Version**: 2.2 Restructured  
**Date**: December 2024  
**Maintainer**: Development Team