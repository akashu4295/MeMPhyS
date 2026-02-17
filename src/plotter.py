import sys
import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter


def main():
    """
    Usage:
    python plotter.py vtk_path variable cmap dimension
    """

    vtk_path = sys.argv[1]
    var_choice = sys.argv[2]
    cmap = sys.argv[3]
    dimension = sys.argv[4]  # "2" or "3"

    mesh = pv.read(vtk_path)

    plot_data = None
    if var_choice in ("u", "v", "w"):
        idx = {"u": 0, "v": 1, "w": 2}[var_choice]
        plot_data = mesh["velocity"][:, idx]
    elif var_choice == "velocity magnitude":
        vectors = mesh["velocity"]
        plot_data = np.linalg.norm(vectors, axis=1)
    elif var_choice == "p":
        plot_data = mesh["pressure"]
    if plot_data is None:
        raise RuntimeError(f"Variable '{var_choice}' not found")

    plotter = BackgroundPlotter()
    plotter.add_mesh(
        mesh,
        scalars=plot_data,
        cmap=cmap,
        scalar_bar_args={"title": var_choice},
    )

    if dimension == "2":
        plotter.view_xy()
        plotter.enable_parallel_projection()
        plotter.reset_camera()

    plotter.show()

    plotter.app.exec_()


if __name__ == "__main__":
    main()