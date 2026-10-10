"""PySide6 desktop interface for the MeMPhyS Navier-Stokes solver."""

import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QProcess, QProcessEnvironment, QSettings, Qt, QUrl
from PySide6.QtGui import QAction, QColor, QDesktopServices, QFont, QFontDatabase, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFontDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.config.constants import APP_FULL_NAME, HELP_URL
from src.core import app_state, logger
from src.qt_model import PARAMETERS, ParameterSpec, read_grid_csv, read_parameter_csv, write_grid_csv, write_parameter_csv
from src.utils.gmsh_bc_manager import get_all_boundary_conditions, read_physical_names_from_msh, set_boundary_condition, write_bc_csv


ROOT = Path(__file__).resolve().parent.parent
CONVERGENCE_PATTERN = re.compile(r"Time step:\s*(\d+),\s*Steady state error:\s*([0-9eE+\-.]+)")


LIGHT_STYLE = """
QWidget { color: #202927; background-color: #f3f5f4; font-size: 10pt; selection-background-color: #007d6c; selection-color: #ffffff; }
QMainWindow, QDialog, QMenuBar, QMenu { background-color: #f3f5f4; color: #202927; }
QTabWidget::pane { background-color: #ffffff; border: 1px solid #c8d2ce; }
QTabBar::tab { color: #202927; background-color: #e4eae7; border: 1px solid #c8d2ce; padding: 7px 12px; }
QTabBar::tab:selected { color: #005d51; background-color: #ffffff; }
QMenuBar::item { padding: 5px 9px; background: transparent; }
QMenuBar::item:selected, QMenu::item:selected { color: #004f44; background-color: #d8ebe5; }
QMenu { border: 1px solid #bdc9c4; }
QMenu::item { padding: 6px 28px 6px 22px; }
QGroupBox { border: 1px solid #c8d2ce; margin-top: 12px; padding: 10px; font-weight: 600; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; }
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit, QPlainTextEdit, QTableWidget { color: #202927; background-color: #ffffff; border: 1px solid #bdc9c4; padding: 5px; selection-background-color: #007d6c; selection-color: #ffffff; }
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled { color: #78837e; background-color: #e8ecea; }
QComboBox::drop-down { border: 0; width: 22px; }
QComboBox QAbstractItemView { color: #202927; background-color: #ffffff; border: 1px solid #bdc9c4; selection-background-color: #007d6c; selection-color: #ffffff; outline: 0; }
QCheckBox::indicator { width: 16px; height: 16px; border: 1px solid #8b9992; background-color: #ffffff; }
QCheckBox::indicator:checked { border-color: #007d6c; background-color: #007d6c; }
QCheckBox::indicator:disabled { border-color: #c8d2ce; background-color: #e8ecea; }
QPushButton, QToolButton { color: #202927; background-color: #e4eae7; border: 1px solid #bdc9c4; padding: 7px 12px; }
QPushButton:hover, QToolButton:hover { background-color: #d5e4df; }
QPushButton:disabled { color: #78837e; background-color: #e8ecea; }
QToolButton:disabled { color: #78837e; background-color: #e8ecea; }
QToolButton#accordionHeader { color: #202927; text-align: left; font-weight: 700; padding: 10px 14px; background-color: #e4eae7; border: 1px solid #bdc9c4; border-radius: 7px; }
QToolButton#accordionHeader:hover { background-color: #d5e4df; border-color: #91aaa0; }
QToolButton#accordionHeader:checked { color: #005d51; background-color: #d5e7e1; border-color: #91b5aa; border-bottom-left-radius: 0; border-bottom-right-radius: 0; }
QWidget#accordionContent { background-color: #ffffff; border: 1px solid #c8d2ce; border-top: 0; border-bottom-left-radius: 7px; border-bottom-right-radius: 7px; }
QPushButton#runButton { color: #ffffff; background-color: #007d6c; border: 0; font-weight: 700; padding: 10px 16px; }
QPushButton#runButton[runState="compiling"] { color: #202927; background-color: #e4ae43; }
QPushButton#runButton[runState="running"] { color: #ffffff; background-color: #16824b; }
QPushButton#runButton[runState="success"] { color: #ffffff; background-color: #16824b; }
QPushButton#runButton[runState="failed"] { color: #ffffff; background-color: #b33d32; }
QPushButton#runButton[runState="stopping"] { color: #202927; background-color: #c3cbc7; }
QPushButton#stopButton { color: #a33027; }
QWidget#runFooter { background-color: #e7ece9; border-top: 1px solid #c8d2ce; }
QHeaderView::section { color: #202927; background-color: #e4eae7; border: 1px solid #c8d2ce; padding: 5px; }
QToolTip { color: #202927; background-color: #ffffff; border: 1px solid #9eaaa4; padding: 4px; }
QScrollBar:vertical { background: #e8ecea; width: 12px; margin: 0; }
QScrollBar:horizontal { background: #e8ecea; height: 12px; margin: 0; }
QScrollBar::handle:vertical { background: #aebbb5; min-height: 24px; }
QScrollBar::handle:vertical:hover { background: #899991; }
QScrollBar::handle:horizontal { background: #aebbb5; min-width: 24px; }
QScrollBar::handle:horizontal:hover { background: #899991; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QStatusBar { color: #202927; background-color: #e7ece9; }
"""

DARK_STYLE = """
QWidget { color: #e1e9e5; background-color: #202725; font-size: 10pt; selection-background-color: #008773; selection-color: #ffffff; }
QMainWindow, QDialog, QMenuBar, QMenu { background-color: #202725; color: #e1e9e5; }
QTabWidget::pane { background-color: #29312e; border: 1px solid #475550; }
QTabBar::tab { color: #e1e9e5; background-color: #303a35; border: 1px solid #475550; padding: 7px 12px; }
QTabBar::tab:selected { color: #8ce0ca; background-color: #29312e; }
QMenuBar::item { padding: 5px 9px; background: transparent; }
QMenuBar::item:selected, QMenu::item:selected { color: #ffffff; background-color: #34564c; }
QMenu { border: 1px solid #53615b; }
QMenu::item { padding: 6px 28px 6px 22px; }
QGroupBox { color: #e1e9e5; border: 1px solid #475550; margin-top: 12px; padding: 10px; font-weight: 600; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; }
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit, QPlainTextEdit, QTableWidget { color: #e1e9e5; background-color: #252d2a; border: 1px solid #53615b; padding: 5px; selection-background-color: #008773; selection-color: #ffffff; }
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled { color: #899791; background-color: #29312e; }
QComboBox::drop-down { border: 0; width: 22px; }
QComboBox QAbstractItemView { color: #e1e9e5; background-color: #252d2a; border: 1px solid #53615b; selection-background-color: #008773; selection-color: #ffffff; outline: 0; }
QCheckBox::indicator { width: 16px; height: 16px; border: 1px solid #718078; background-color: #252d2a; }
QCheckBox::indicator:checked { border-color: #008773; background-color: #008773; }
QCheckBox::indicator:disabled { border-color: #475550; background-color: #29312e; }
QPushButton, QToolButton { color: #e1e9e5; background-color: #35403b; border: 1px solid #53615b; padding: 7px 12px; }
QPushButton:hover, QToolButton:hover { background-color: #3e5149; }
QPushButton:disabled { color: #899791; background-color: #29312e; }
QToolButton:disabled { color: #899791; background-color: #29312e; }
QToolButton#accordionHeader { color: #e1e9e5; text-align: left; font-weight: 700; padding: 10px 14px; background-color: #303a35; border: 1px solid #53615b; border-radius: 7px; }
QToolButton#accordionHeader:hover { background-color: #3e5149; border-color: #72857c; }
QToolButton#accordionHeader:checked { color: #8ce0ca; background-color: #314a41; border-color: #56766a; border-bottom-left-radius: 0; border-bottom-right-radius: 0; }
QWidget#accordionContent { background-color: #252d2a; border: 1px solid #475550; border-top: 0; border-bottom-left-radius: 7px; border-bottom-right-radius: 7px; }
QPushButton#runButton { color: #ffffff; background-color: #008773; border: 0; font-weight: 700; padding: 10px 16px; }
QPushButton#runButton[runState="compiling"] { color: #202725; background-color: #e4ae43; }
QPushButton#runButton[runState="running"] { color: #ffffff; background-color: #16824b; }
QPushButton#runButton[runState="success"] { color: #ffffff; background-color: #16824b; }
QPushButton#runButton[runState="failed"] { color: #ffffff; background-color: #b33d32; }
QPushButton#runButton[runState="stopping"] { color: #e1e9e5; background-color: #59635e; }
QPushButton#stopButton { color: #ff9b90; }
QWidget#runFooter { background-color: #29312e; border-top: 1px solid #475550; }
QHeaderView::section { color: #e1e9e5; background-color: #303a35; border: 1px solid #475550; padding: 5px; }
QToolTip { color: #e1e9e5; background-color: #303a35; border: 1px solid #62716a; padding: 4px; }
QScrollBar:vertical { background: #252d2a; width: 12px; margin: 0; }
QScrollBar:horizontal { background: #252d2a; height: 12px; margin: 0; }
QScrollBar::handle:vertical { background: #53615b; min-height: 24px; }
QScrollBar::handle:vertical:hover { background: #718078; }
QScrollBar::handle:horizontal { background: #53615b; min-width: 24px; }
QScrollBar::handle:horizontal:hover { background: #718078; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QStatusBar { color: #e1e9e5; background-color: #29312e; }
"""


def _make_input(spec: ParameterSpec, value):
    if spec.kind == "bool":
        widget = QCheckBox()
        widget.setChecked(bool(int(value)))
    elif spec.kind == "choice":
        widget = QComboBox()
        for title, code in spec.choices:
            widget.addItem(title, code)
        index = widget.findData(int(value))
        widget.setCurrentIndex(max(index, 0))
    elif spec.kind == "int":
        widget = QSpinBox()
        widget.setRange(-2_000_000_000, 2_000_000_000)
        widget.setValue(int(value))
    elif spec.kind == "float":
        widget = QDoubleSpinBox()
        widget.setDecimals(12)
        widget.setRange(-1e100, 1e100)
        widget.setSingleStep(max(abs(float(value)) * 0.1, 0.01))
        widget.setValue(float(value))
    else:
        widget = QLineEdit(str(value))
    if spec.tooltip:
        widget.setToolTip(spec.tooltip)
    return widget


def _widget_value(spec: ParameterSpec, widget):
    if spec.kind == "bool":
        return int(widget.isChecked())
    if spec.kind == "choice":
        return int(widget.currentData())
    if spec.kind == "int":
        return widget.value()
    if spec.kind == "float":
        return widget.value()
    return widget.text().strip()


class BoundaryDialog(QDialog):
    """Assign BC types and values to physical names from a selected mesh."""

    def __init__(self, mesh_path: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Boundary conditions")
        self.resize(720, 420)
        layout = QVBoxLayout(self)
        names = read_physical_names_from_msh(mesh_path)
        self.table = QTableWidget(len(names), 3)
        self.table.setHorizontalHeaderLabels(("Physical name", "Condition", "Values (u=0,v=0,p=0)"))
        self.table.horizontalHeader().setStretchLastSection(True)
        current = get_all_boundary_conditions()
        self._type_boxes = []
        for row, name in enumerate(names):
            self.table.setItem(row, 0, QTableWidgetItem(name))
            type_box = QComboBox()
            type_box.addItems(("velocity_inlet", "pressure_outlet", "wall", "symmetry", "periodic", "outflow", "interior"))
            existing = current.get(name, {})
            type_box.setCurrentText(existing.get("type", "wall"))
            self._type_boxes.append(type_box)
            self.table.setCellWidget(row, 1, type_box)
            values = existing.get("variables", {})
            self.table.setItem(row, 2, QTableWidgetItem(", ".join(f"{key}={value}" for key, value in values.items())))
        layout.addWidget(self.table)
        if not names:
            layout.addWidget(QLabel("No physical names were found in this mesh."))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def save_assignments(self):
        for row in range(self.table.rowCount()):
            name = self.table.item(row, 0).text()
            kind = self.table.cellWidget(row, 1).currentText()
            raw_values = self.table.item(row, 2)
            variables = {}
            if raw_values and raw_values.text().strip():
                for part in raw_values.text().split(","):
                    key, separator, value = part.partition("=")
                    if separator:
                        try:
                            variables[key.strip()] = float(value.strip())
                        except ValueError as error:
                            raise ValueError(f"Invalid value for {name}: {part.strip()}") from error
            set_boundary_condition(name, kind, variables)
        return write_bc_csv()


class CollapsibleSection(QWidget):
    def __init__(self, title: str, expanded: bool = False, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.header = QToolButton()
        self.header.setObjectName("accordionHeader")
        self.header.setText(title)
        self.header.setCheckable(True)
        self.header.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.header.setMinimumHeight(42)
        self.header.toggled.connect(self._set_expanded)
        layout.addWidget(self.header)
        self.content = QWidget()
        self.content.setObjectName("accordionContent")
        self.content_layout = QVBoxLayout(self.content)
        self.content_layout.setContentsMargins(12, 10, 12, 12)
        self.content_layout.setSpacing(10)
        layout.addWidget(self.content)
        self._set_expanded(expanded)

    def set_expanded(self, expanded: bool):
        self.header.setChecked(expanded)

    def _set_expanded(self, expanded: bool):
        self.header.setArrowType(Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow)
        self.content.setVisible(expanded)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_FULL_NAME} | Meshless Flow Solver")
        screen = QApplication.primaryScreen()
        available = screen.availableGeometry() if screen else None
        max_width = min(1440, int(available.width() * 0.96)) if available else 1440
        max_height = min(940, int(available.height() * 0.92)) if available else 860
        self.setMaximumSize(max_width, max_height)
        self.resize(max_width, max_height)
        self.settings = QSettings("MeMPhyS", "MeMPhySQt")
        self.parameters = read_parameter_csv(ROOT / "flow_parameters.csv")
        self.mesh_files = read_grid_csv(ROOT / "grid_filenames.csv")
        self.parameter_widgets = {}
        self.mesh_edits = []
        self.convergence_steps = []
        self.convergence_errors = []
        self.process = None
        self.process_lines = ""
        self.solver_executable = None
        self.init_path = ""
        self._build_ui()
        self._build_menus()
        self._restore_preferences()
        self._refresh_conditional_fields()
        self._drain_timer = self.startTimer(100)
        logger.set_enable_gui(True)
        logger.info("MeMPhyS Qt interface ready")

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer = QVBoxLayout(central)
        outer.setContentsMargins(8, 8, 8, 8)
        outer.setSpacing(8)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)

        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(6)

        self.left_scroll = QScrollArea()
        self.left_scroll.setWidgetResizable(True)
        self.left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        accordion = QWidget()
        self.accordion_layout = QVBoxLayout(accordion)
        self.accordion_layout.setContentsMargins(6, 6, 6, 6)
        self.accordion_layout.setSpacing(5)
        self.accordion_sections = []

        setup = CollapsibleSection("Run Setup")
        setup.content_layout.addWidget(self._build_setup_section())
        self._add_accordion_section(setup)

        sections = (
            "Meshless discretization",
            "Flow and time",
            "Time integration",
            "Pressure solver",
            "Advanced numerics",
            "Compressible properties",
        )
        for title in sections:
            section = CollapsibleSection(title)
            section.content_layout.addWidget(self._build_parameter_section(title))
            if title == "Compressible properties":
                self.compressible_section = section
            self._add_accordion_section(section)
        self.accordion_layout.addStretch(1)
        self.left_scroll.setWidget(accordion)
        left_layout.addWidget(self.left_scroll, 1)
        left_layout.addWidget(self._build_run_footer())
        left_panel.setMinimumWidth(350)
        left_panel.setMaximumWidth(620)
        splitter.addWidget(left_panel)
        splitter.addWidget(self._build_output_panel())
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        splitter.setSizes((520, 880))
        self._set_run_state("idle", "Ready")
        self.statusBar().showMessage("Ready")

    def _add_accordion_section(self, section):
        self.accordion_sections.append(section)
        section.header.toggled.connect(
            lambda expanded, active=section: self._accordion_toggled(active, expanded)
        )
        self.accordion_layout.addWidget(section)

    def _accordion_toggled(self, active, expanded):
        if expanded:
            for section in self.accordion_sections:
                if section is not active:
                    section.set_expanded(False)

    def _build_setup_section(self):
        page = QWidget()
        body_layout = QVBoxLayout(page)
        body_layout.setContentsMargins(0, 0, 0, 0)

        geometry = QGroupBox("Geometry and mesh hierarchy")
        geometry_layout = QVBoxLayout(geometry)
        row = QHBoxLayout()
        self.geometry_path = QLineEdit()
        self.geometry_path.setPlaceholderText("Optional .geo source")
        browse_geo = QPushButton("Browse .geo")
        browse_geo.clicked.connect(self._browse_geometry)
        launch_gmsh = QPushButton("Launch Gmsh")
        launch_gmsh.clicked.connect(self._launch_gmsh)
        row.addWidget(self.geometry_path, 1)
        row.addWidget(browse_geo)
        row.addWidget(launch_gmsh)
        geometry_layout.addLayout(row)
        geometry_layout.addWidget(QLabel("Mesh files are ordered finest to coarsest."))
        self.mesh_count = QSpinBox()
        self.mesh_count.setRange(1, 10)
        self.mesh_count.setValue(max(1, len(self.mesh_files)))
        self.mesh_count.valueChanged.connect(self._set_mesh_count)
        count_row = QHBoxLayout()
        count_row.addWidget(QLabel("Mesh levels"))
        count_row.addWidget(self.mesh_count)
        count_row.addStretch(1)
        geometry_layout.addLayout(count_row)
        self.mesh_list = QWidget()
        self.mesh_list_layout = QVBoxLayout(self.mesh_list)
        self.mesh_list_layout.setContentsMargins(0, 0, 0, 0)
        geometry_layout.addWidget(self.mesh_list)
        self._set_mesh_count(self.mesh_count.value())
        mesh_buttons = QHBoxLayout()
        add_mesh = QPushButton("Add mesh file")
        add_mesh.clicked.connect(self._add_mesh_file)
        self.bc_button = QPushButton("Boundary conditions")
        self.bc_button.clicked.connect(self._edit_boundary_conditions)
        mesh_buttons.addWidget(add_mesh)
        mesh_buttons.addWidget(self.bc_button)
        geometry_layout.addLayout(mesh_buttons)
        body_layout.addWidget(geometry)

        initial = QGroupBox("Initial conditions and restart")
        initial_layout = QFormLayout(initial)
        self.use_init = QCheckBox("Compile a custom initialization .c file")
        self.init_edit = QLineEdit()
        init_row = self._file_row(self.init_edit, "C source (*.c)")
        initial_layout.addRow(self.use_init)
        initial_layout.addRow("Initialization file", init_row)
        restart_spec = next(spec for spec in PARAMETERS if spec.key == "restart_filename")
        self.parameter_widgets[restart_spec.key] = _make_input(restart_spec, self.parameters[restart_spec.key])
        initial_layout.addRow("Restart VTK", self._file_row(self.parameter_widgets[restart_spec.key], "VTK files (*.vtk)"))
        initial_layout.addRow(self._widget_label("Restart from previous run"), self.parameter_widgets["restart"] if "restart" in self.parameter_widgets else self._add_restart_checkbox(initial_layout))
        body_layout.addWidget(initial)

        compiler_group = QGroupBox("Compiler")
        compiler_layout = QVBoxLayout(compiler_group)
        self.gpu_check = QCheckBox("Use NVIDIA HPC SDK / OpenACC GPU")
        compiler_layout.addWidget(self.gpu_check)
        body_layout.addWidget(compiler_group)
        return page

    def _build_run_footer(self):
        footer = QWidget()
        footer.setObjectName("runFooter")
        layout = QHBoxLayout(footer)
        layout.setContentsMargins(12, 8, 12, 8)
        self.run_button = QPushButton("Compile and run solver")
        self.run_button.setObjectName("runButton")
        self.run_button.setMinimumSize(220, 40)
        self.run_button.clicked.connect(self._start_solver)
        layout.addWidget(self.run_button)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setObjectName("stopButton")
        self.stop_button.setMinimumSize(90, 40)
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop_solver)
        layout.addWidget(self.stop_button)
        layout.addStretch(1)
        return footer

    def _add_restart_checkbox(self, form):
        spec = next(spec for spec in PARAMETERS if spec.key == "restart")
        checkbox = _make_input(spec, self.parameters[spec.key])
        self.parameter_widgets[spec.key] = checkbox
        checkbox.toggled.connect(self._refresh_conditional_fields)
        return checkbox

    @staticmethod
    def _widget_label(text):
        return QLabel(text)

    def _file_row(self, edit, file_filter):
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        browse = QPushButton("Browse")
        browse.clicked.connect(lambda: self._browse_into(edit, file_filter))
        layout.addWidget(edit, 1)
        layout.addWidget(browse)
        return row

    def _build_parameter_section(self, section):
        content = QWidget()
        form = QFormLayout(content)
        form.setContentsMargins(0, 0, 0, 0)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        for spec in PARAMETERS:
            if spec.section != section or spec.key in {"restart", "restart_filename"}:
                continue
            widget = _make_input(spec, self.parameters[spec.key])
            self.parameter_widgets[spec.key] = widget
            if spec.key in {"use_timple", "use_compressible_flow"}:
                if spec.kind == "choice":
                    widget.currentIndexChanged.connect(self._refresh_conditional_fields)
                else:
                    widget.toggled.connect(self._refresh_conditional_fields)
            form.addRow(spec.label, widget)
        return content

    def _build_output_panel(self):
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.output_splitter = QSplitter(Qt.Orientation.Vertical)
        self.output_splitter.setChildrenCollapsible(False)
        self.output_splitter.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        convergence = QWidget()
        convergence_layout = QVBoxLayout(convergence)
        controls = QHBoxLayout()
        self.vtk_edit = QLineEdit("Solution.vtk")
        controls.addWidget(QLabel("VTK"))
        controls.addWidget(self.vtk_edit, 1)
        browse_vtk = QPushButton("Browse")
        browse_vtk.clicked.connect(lambda: self._browse_into(self.vtk_edit, "VTK files (*.vtk)"))
        controls.addWidget(browse_vtk)
        self.plot_variable = QComboBox()
        self.plot_variable.addItems(("velocity magnitude", "u", "v", "w", "p"))
        controls.addWidget(self.plot_variable)
        self.colormap = QComboBox()
        self.colormap.addItems(("viridis", "plasma", "inferno", "magma", "cividis", "coolwarm", "turbo"))
        controls.addWidget(self.colormap)
        plot_button = QPushButton("Open plot")
        plot_button.clicked.connect(self._open_plot)
        controls.addWidget(plot_button)
        paraview_button = QPushButton("ParaView")
        paraview_button.clicked.connect(self._open_paraview)
        controls.addWidget(paraview_button)
        convergence_layout.addLayout(controls)
        try:
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
            from matplotlib.figure import Figure

            self.figure = Figure(figsize=(7, 4), tight_layout=True)
            self.axes = self.figure.add_subplot(111)
            self.axes.set_xlabel("Time step")
            self.axes.set_ylabel("Steady-state error")
            self.axes.set_yscale("log")
            self.axes.grid(True, which="both", alpha=0.22)
            self.canvas = FigureCanvasQTAgg(self.figure)
            convergence_layout.addWidget(self.canvas, 1)
        except ImportError:
            self.canvas = None
            self.axes = None
            convergence_layout.addWidget(QLabel("Install matplotlib to show the convergence plot."))
        self.output_splitter.addWidget(convergence)

        logs = QWidget()
        logs_layout = QVBoxLayout(logs)
        logs_layout.setContentsMargins(0, 0, 0, 0)
        logs_layout.setSpacing(6)
        log_header = QHBoxLayout()
        log_header.addWidget(QLabel("Solver log"))
        log_header.addStretch(1)
        clear_logs = QPushButton("Clear log")
        clear_logs.clicked.connect(self._clear_logs)
        log_header.addWidget(clear_logs)
        logs_layout.addLayout(log_header)
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.log_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.log_view.setMinimumHeight(120)
        logs_layout.addWidget(self.log_view)
        self.output_splitter.addWidget(logs)
        self.output_splitter.setStretchFactor(0, 3)
        self.output_splitter.setStretchFactor(1, 2)
        self.output_splitter.setSizes((560, 300))
        layout.addWidget(self.output_splitter)
        return panel

    def _clear_logs(self):
        self.log_view.clear()
        logger.clear()

    def _build_menus(self):
        file_menu = self.menuBar().addMenu("File")
        self._action(file_menu, "Open configuration...", self._open_configuration)
        self._action(file_menu, "Save configuration...", self._save_configuration)
        file_menu.addSeparator()
        self._action(file_menu, "Exit", self.close)
        settings = self.menuBar().addMenu("Settings")
        self._action(settings, "Choose font...", self._choose_font)
        self._action(settings, "Reset font", self._reset_font)
        settings.addSeparator()
        self._action(settings, "Toggle light/dark theme", self._toggle_theme)
        help_menu = self.menuBar().addMenu("Help")
        self._action(help_menu, "Documentation", lambda: QDesktopServices.openUrl(QUrl(HELP_URL)))
        self._action(help_menu, "About", self._about)

    @staticmethod
    def _action(menu, label, callback):
        action = QAction(label, menu)
        action.triggered.connect(callback)
        menu.addAction(action)
        return action

    def _set_mesh_count(self, count):
        current = [edit.text() for edit in self.mesh_edits]
        while self.mesh_list_layout.count():
            item = self.mesh_list_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.mesh_edits = []
        for index in range(count):
            row = QWidget()
            row_layout = QHBoxLayout(row)
            row_layout.setContentsMargins(0, 2, 0, 2)
            label = QLabel(f"Level {index + 1}")
            edit = QLineEdit(current[index] if index < len(current) else (self.mesh_files[index] if index < len(self.mesh_files) else ""))
            edit.setPlaceholderText("Select Gmsh .msh file")
            browse = QPushButton("Browse")
            browse.clicked.connect(lambda _checked=False, target=edit: self._browse_into(target, "Gmsh meshes (*.msh)"))
            row_layout.addWidget(label)
            row_layout.addWidget(edit, 1)
            row_layout.addWidget(browse)
            self.mesh_list_layout.addWidget(row)
            self.mesh_edits.append(edit)

    def _add_mesh_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select mesh", str(ROOT), "Gmsh meshes (*.msh)")
        if not path:
            return
        if self.mesh_edits[-1].text():
            if self.mesh_count.value() == self.mesh_count.maximum():
                QMessageBox.information(self, "Mesh level limit", "The maximum number of mesh levels is 10.")
                return
            self.mesh_count.setValue(self.mesh_count.value() + 1)
        self.mesh_edits[-1].setText(path)

    def _browse_geometry(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Gmsh geometry", str(ROOT), "Gmsh geometry (*.geo *.geo_unrolled)")
        if path:
            self.geometry_path.setText(path)

    def _launch_gmsh(self):
        executable = shutil.which("gmsh")
        if not executable:
            QMessageBox.warning(self, "Gmsh not found", "Install Gmsh or add it to PATH.")
            return
        args = [self.geometry_path.text()] if self.geometry_path.text() else []
        try:
            subprocess.Popen([executable, *args], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self._append_log("Launched Gmsh")
        except OSError as error:
            QMessageBox.critical(self, "Could not launch Gmsh", str(error))

    def _edit_boundary_conditions(self):
        mesh_path = next((edit.text() for edit in self.mesh_edits if edit.text()), "")
        if not mesh_path:
            QMessageBox.information(self, "Select a mesh", "Choose a mesh file before assigning boundary conditions.")
            return
        dialog = BoundaryDialog(mesh_path, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            try:
                if dialog.save_assignments():
                    self._append_log("Boundary conditions saved to bc.csv")
            except ValueError as error:
                QMessageBox.warning(self, "Invalid boundary value", str(error))

    def _browse_into(self, edit, file_filter):
        path, _ = QFileDialog.getOpenFileName(self, "Select file", str(ROOT), file_filter)
        if path:
            edit.setText(path)

    def _refresh_conditional_fields(self, *_):
        timple = self.parameter_widgets.get("use_timple")
        implicit = timple.isChecked() if timple else bool(self.parameters.get("use_timple", 0))
        for key in ("iter_momentum", "iter_timple"):
            widget = self.parameter_widgets.get(key)
            if widget:
                widget.setEnabled(implicit)
        restart = self.parameter_widgets.get("restart")
        if restart:
            self.parameter_widgets["restart_filename"].parentWidget().setVisible(bool(restart.isChecked()))
        compressible = self.parameter_widgets.get("use_compressible_flow")
        if compressible:
            self.compressible_section.setEnabled(compressible.isChecked())

    def _collect_parameters(self):
        values = {}
        for spec in PARAMETERS:
            widget = self.parameter_widgets.get(spec.key)
            values[spec.key] = _widget_value(spec, widget) if widget else self.parameters[spec.key]
        return values

    def _start_solver(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            return
        meshes = [edit.text().strip() for edit in self.mesh_edits]
        if not meshes or any(not path for path in meshes):
            QMessageBox.warning(self, "Mesh files required", "Select a valid mesh file for each configured level.")
            return
        for mesh in meshes:
            if not Path(mesh).is_file() or Path(mesh).suffix.lower() != ".msh":
                QMessageBox.warning(self, "Invalid mesh file", f"Mesh file does not exist or is not .msh:\n{mesh}")
                return
        init_file = self.init_edit.text().strip()
        if self.use_init.isChecked() and (not init_file or not Path(init_file).is_file()):
            QMessageBox.warning(self, "Initialization file required", "Select an existing custom initialization .c file.")
            return
        try:
            values = self._collect_parameters()
            if values["restart"] and not Path(values["restart_filename"]).is_file():
                QMessageBox.warning(self, "Restart file required", "Select an existing restart VTK file.")
                return
            write_parameter_csv(values, ROOT / "flow_parameters.csv")
            write_grid_csv(meshes, ROOT / "grid_filenames.csv")
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Configuration error", str(error))
            return
        self.parameters = values
        self.convergence_steps.clear()
        self.convergence_errors.clear()
        self._draw_convergence()
        self.process_lines = ""
        compiler = "nvc" if self.gpu_check.isChecked() else "gcc"
        missing_tools = [compiler] if not shutil.which(compiler) else []
        if missing_tools:
            QMessageBox.critical(self, "Compiler not found", f"Could not find {compiler} in PATH.")
            return
        try:
            command, arguments, executable_name = self._compile_command(
                compiler,
                init_file if self.use_init.isChecked() else "",
            )
            self.solver_executable = ROOT / executable_name
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Build configuration error", str(error))
            return
        self._set_run_state("compiling", "Compiling solver...")
        self.statusBar().showMessage("Compiling solver")
        self._append_log(f"Compile command: {command} {' '.join(arguments)}")
        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(ROOT))
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        environment = QProcessEnvironment.systemEnvironment()
        include_paths = (str(ROOT), str(ROOT / "src" / "lib"), environment.value("CPATH"))
        environment.insert("CPATH", os.pathsep.join(path for path in include_paths if path))
        self.process.setProcessEnvironment(environment)
        self.process.readyReadStandardOutput.connect(self._read_process_output)
        self.process.finished.connect(self._compile_finished)
        self.process.errorOccurred.connect(self._process_error)
        self.process.start(command, arguments)

    @staticmethod
    def _compile_command(compiler, init_file, platform_name=None):
        platform_name = platform_name or sys.platform
        is_windows = platform_name.startswith("win")
        source_directory = ROOT / "src" / "lib"
        if not source_directory.is_dir():
            raise FileNotFoundError(f"Solver source directory not found: {source_directory}")

        source_files = sorted(source_directory.glob("*.c"))
        main_source = ROOT / "mg_NS_solver.c"
        if not main_source.is_file():
            raise FileNotFoundError(f"Solver entry point not found: {main_source}")
        source_files.append(main_source)
        if init_file:
            init_source = Path(init_file).resolve()
            if not init_source.is_file():
                raise FileNotFoundError(f"Initialization source not found: {init_source}")
            source_files.append(init_source)

        compiler_name = Path(compiler).stem.lower()
        if compiler_name == "nvc":
            flags = ["-acc", "-mp", "-Minfo=accel", "-O3"]
        else:
            flags = ["-O3", "-Wall", "-Wpointer-arith", "-fopenmp"]
            if not is_windows:
                flags.append("-march=native")
            elif compiler_name.endswith("gcc"):
                flags.append("-static")

        if is_windows and init_file:
            init_source = Path(init_file).resolve()
            init_code = init_source.read_text(encoding="utf-8", errors="replace")
            init_code = re.sub(r"/\*.*?\*/|//[^\n]*", "", init_code, flags=re.DOTALL)
            hooks = {
                "initial_conditions": "MEMPHYS_CUSTOM_INITIAL_CONDITIONS",
                "boundary_conditions": "MEMPHYS_CUSTOM_BOUNDARY_CONDITIONS",
            }
            for function_name, macro_name in hooks.items():
                definition = rf"\bvoid\s+{function_name}\s*\([^;{{}}]*\)\s*\{{"
                if re.search(definition, init_code, flags=re.DOTALL):
                    flags.append(f"-D{macro_name}")

        flags.append("-lm")

        executable_name = "memphys_solver.exe" if is_windows else "memphys_solver"
        output_path = ROOT / executable_name
        arguments = [str(path.resolve()) for path in source_files]
        arguments.extend((*flags, "-o", str(output_path)))
        return compiler, arguments, executable_name

    def _compile_finished(self, exit_code, _exit_status):
        if exit_code != 0:
            self._finish_solver(False, f"Compilation failed (exit code {exit_code})")
            return
        self._set_run_state("running", "Solver is running")
        self.statusBar().showMessage("Solver running")
        self._append_log("Compilation successful; starting solver")
        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(ROOT))
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_process_output)
        self.process.finished.connect(self._solver_finished)
        self.process.errorOccurred.connect(self._process_error)
        if self.solver_executable is None:
            self._finish_solver(False, "Solver executable path was not configured")
            return
        self.process.start(str(self.solver_executable), [])

    def _read_process_output(self):
        if not self.process:
            return
        text = bytes(self.process.readAllStandardOutput()).decode(errors="replace")
        if not text:
            return
        self._append_log(text.rstrip())
        self.process_lines += text
        lines = self.process_lines.splitlines(keepends=True)
        self.process_lines = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            self.process_lines = lines.pop()
        for line in lines:
            match = CONVERGENCE_PATTERN.search(line)
            if match:
                self.convergence_steps.append(int(match.group(1)))
                self.convergence_errors.append(float(match.group(2)))
                self._draw_convergence()

    def _draw_convergence(self):
        if not self.axes:
            return
        background = "#252d2a" if getattr(self, "_is_dark_theme", False) else "#ffffff"
        foreground = "#e1e9e5" if getattr(self, "_is_dark_theme", False) else "#202927"
        grid_color = "#53615b" if getattr(self, "_is_dark_theme", False) else "#c8d2ce"
        self.figure.set_facecolor(background)
        self.axes.clear()
        self.axes.set_facecolor(background)
        self.axes.plot(self.convergence_steps, self.convergence_errors, color="#008773", linewidth=1.7)
        self.axes.set_xlabel("Time step", color=foreground)
        self.axes.set_ylabel("Steady-state error", color=foreground)
        self.axes.set_yscale("log")
        self.axes.tick_params(colors=foreground)
        for spine in self.axes.spines.values():
            spine.set_color(grid_color)
        self.axes.grid(True, which="both", color=grid_color, alpha=0.4)
        self.canvas.draw_idle()

    def _solver_finished(self, exit_code, _exit_status):
        self._finish_solver(exit_code == 0, "Solver completed" if exit_code == 0 else f"Solver exited with code {exit_code}")
        if exit_code == 0:
            try:
                from src.utils.output_manager import organize_solver_outputs

                organize_solver_outputs()
            except Exception as error:
                logger.warning(f"Could not organize solver outputs: {error}")

    def _finish_solver(self, success, message):
        self._append_log(message)
        state = "success" if success else "stopped" if message == "Solver stopped by user" else "failed"
        self._set_run_state(state, message)
        self.statusBar().showMessage(message, 10000)
        self.process = None

    def _set_run_state(self, state, message):
        button_text = {
            "idle": "Compile and run solver",
            "compiling": "Compiling...",
            "running": "Solver running...",
            "success": "Run again",
            "failed": "Try again",
            "stopping": "Stopping...",
            "stopped": "Run again",
        }
        self.run_button.setText(button_text[state])
        self.run_button.setProperty("runState", state)
        self.run_button.setEnabled(state in {"idle", "success", "failed", "stopped"})
        self.stop_button.setEnabled(state in {"compiling", "running"})
        self.run_button.style().unpolish(self.run_button)
        self.run_button.style().polish(self.run_button)
        self.run_button.update()
        self.statusBar().showMessage(message)

    def _process_error(self, error):
        if self.process and self.process.state() == QProcess.ProcessState.NotRunning:
            self._finish_solver(False, f"Could not start process: {self.process.errorString()}")

    def _stop_solver(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self._set_run_state("stopping", "Stopping solver...")
            self.process.terminate()
            if not self.process.waitForFinished(3000):
                self.process.kill()
            self._finish_solver(False, "Solver stopped by user")

    def _open_plot(self):
        vtk_path = Path(self.vtk_edit.text().strip())
        if not vtk_path.is_file():
            latest = ROOT / "output"
            candidates = sorted(latest.glob("**/" + vtk_path.name)) if latest.is_dir() else []
            if candidates:
                vtk_path = candidates[-1]
            else:
                QMessageBox.warning(self, "VTK file not found", str(vtk_path))
                return
        script = """import sys, numpy as np, pyvista as pv
from pyvistaqt import BackgroundPlotter
mesh = pv.read(sys.argv[1])
name, cmap = sys.argv[2], sys.argv[3]
if name in ('u', 'v', 'w'):
    data = mesh['velocity'][:, {'u': 0, 'v': 1, 'w': 2}[name]]
elif name == 'velocity magnitude':
    data = np.linalg.norm(mesh['velocity'], axis=1)
else:
    data = mesh['pressure']
plotter = BackgroundPlotter()
plotter.add_mesh(mesh, scalars=data, cmap=cmap, scalar_bar_args={'title': name})
plotter.show()
plotter.app.exec()
"""
        try:
            subprocess.Popen([sys.executable, "-c", script, str(vtk_path), self.plot_variable.currentText(), self.colormap.currentText()])
        except OSError as error:
            QMessageBox.critical(self, "Could not open plot", str(error))

    def _open_paraview(self):
        vtk_path = self.vtk_edit.text().strip()
        executable = shutil.which("paraview")
        if not Path(vtk_path).is_file():
            QMessageBox.warning(self, "VTK file not found", vtk_path)
        elif not executable:
            QMessageBox.information(self, "ParaView not found", "Install ParaView or add it to PATH.")
        else:
            subprocess.Popen([executable, vtk_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _open_configuration(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open configuration", str(ROOT), "MeMPhyS configuration (*.json)")
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as stream:
                data = json.load(stream)
            for spec in PARAMETERS:
                if spec.key not in data.get("parameters", {}):
                    continue
                self._set_widget_value(spec, self.parameter_widgets.get(spec.key), data["parameters"][spec.key])
            self.mesh_count.setValue(max(1, min(10, len(data.get("mesh_files", [])))))
            for edit, mesh in zip(self.mesh_edits, data.get("mesh_files", [])):
                edit.setText(mesh)
            self.use_init.setChecked(bool(data.get("use_init", False)))
            self.init_edit.setText(data.get("init_file", ""))
            self._refresh_conditional_fields()
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.critical(self, "Could not load configuration", str(error))

    def _save_configuration(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save configuration", str(ROOT / "memphys_config.json"), "MeMPhyS configuration (*.json)")
        if not path:
            return
        data = {"parameters": self._collect_parameters(), "mesh_files": [edit.text() for edit in self.mesh_edits], "use_init": self.use_init.isChecked(), "init_file": self.init_edit.text()}
        try:
            with open(path, "w", encoding="utf-8") as stream:
                json.dump(data, stream, indent=2)
        except OSError as error:
            QMessageBox.critical(self, "Could not save configuration", str(error))

    @staticmethod
    def _set_widget_value(spec, widget, value):
        if widget is None:
            return
        if spec.kind == "bool":
            widget.setChecked(bool(value))
        elif spec.kind == "choice":
            index = widget.findData(int(value))
            if index >= 0:
                widget.setCurrentIndex(index)
        elif spec.kind in {"int", "float"}:
            widget.setValue(value)
        else:
            widget.setText(str(value))

    def _choose_font(self):
        accepted, font = QFontDialog.getFont(self.font(), self, "Choose application font")
        if accepted:
            QApplication.setFont(font)
            self.settings.setValue("font", font.toString())

    def _reset_font(self):
        font = QApplication.font()
        font.setPointSize(10)
        QApplication.setFont(font)
        self.settings.remove("font")

    def _toggle_theme(self):
        dark = self.settings.value("darkTheme", False, bool)
        self._apply_theme(not dark)

    def _apply_theme(self, dark):
        app = QApplication.instance()
        if dark:
            colors = {
                QPalette.ColorRole.Window: "#202725",
                QPalette.ColorRole.WindowText: "#e1e9e5",
                QPalette.ColorRole.Base: "#252d2a",
                QPalette.ColorRole.AlternateBase: "#303a35",
                QPalette.ColorRole.ToolTipBase: "#303a35",
                QPalette.ColorRole.ToolTipText: "#e1e9e5",
                QPalette.ColorRole.Text: "#e1e9e5",
                QPalette.ColorRole.Button: "#35403b",
                QPalette.ColorRole.ButtonText: "#e1e9e5",
                QPalette.ColorRole.BrightText: "#ffffff",
                QPalette.ColorRole.Highlight: "#008773",
                QPalette.ColorRole.HighlightedText: "#ffffff",
                QPalette.ColorRole.PlaceholderText: "#9ca9a3",
                QPalette.ColorRole.Link: "#70d4bd",
            }
            stylesheet = DARK_STYLE
        else:
            colors = {
                QPalette.ColorRole.Window: "#f3f5f4",
                QPalette.ColorRole.WindowText: "#202927",
                QPalette.ColorRole.Base: "#ffffff",
                QPalette.ColorRole.AlternateBase: "#e8ecea",
                QPalette.ColorRole.ToolTipBase: "#ffffff",
                QPalette.ColorRole.ToolTipText: "#202927",
                QPalette.ColorRole.Text: "#202927",
                QPalette.ColorRole.Button: "#e4eae7",
                QPalette.ColorRole.ButtonText: "#202927",
                QPalette.ColorRole.BrightText: "#ffffff",
                QPalette.ColorRole.Highlight: "#007d6c",
                QPalette.ColorRole.HighlightedText: "#ffffff",
                QPalette.ColorRole.PlaceholderText: "#66736d",
                QPalette.ColorRole.Link: "#006b5c",
            }
            stylesheet = LIGHT_STYLE
        palette = QPalette()
        for role, color in colors.items():
            palette.setColor(role, QColor(color))
        app.setPalette(palette)
        app.setStyleSheet(stylesheet)
        self.settings.setValue("darkTheme", dark)
        self._is_dark_theme = dark
        self._draw_convergence()

    def _restore_preferences(self):
        font_data = self.settings.value("font", "")
        if font_data:
            font = QFont()
            if font.fromString(font_data) and font.family() in QFontDatabase.families():
                QApplication.setFont(font)
            else:
                self.settings.remove("font")
        self._apply_theme(self.settings.value("darkTheme", False, bool))

    def _about(self):
        QMessageBox.about(self, "About MeMPhyS", f"<b>{APP_FULL_NAME}</b><br>Meshless Multi-Physics Solver<br><br>PySide6 desktop interface")

    def _append_log(self, message):
        if not message:
            return
        self.log_view.append(message.replace("\n", "<br>"))

    def timerEvent(self, event):
        for message in logger.drain_gui_queue():
            self._append_log(message)

    def closeEvent(self, event):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            answer = QMessageBox.question(self, "Solver running", "Stop the running solver and exit?")
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self._stop_solver()
        self.settings.setValue("windowSize", self.size())
        logger.info("MeMPhyS GUI closed")
        event.accept()


def main():
    os.chdir(ROOT)
    app = QApplication(sys.argv)
    app.setOrganizationName("MeMPhyS")
    app.setApplicationName("MeMPhyS")
    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.GeneralFont)
    font.setPointSize(10)
    app.setFont(font)
    window = MainWindow()
    saved_size = window.settings.value("windowSize")
    if saved_size:
        window.resize(
            min(saved_size.width(), window.maximumWidth()),
            min(saved_size.height(), window.maximumHeight()),
        )
    window.show()
    return app.exec()