# src/plotter.py
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





# src/plotter.py
# import socket
# import json
# import threading
# import numpy as np
# import pyvista as pv
# import sys
# import traceback

# from pyvistaqt import BackgroundPlotter
# from qtpy import QtWidgets   # ✅ THE correct import

# HOST = "127.0.0.1"
# PORT = 50555


# class PlotServer:
#     def __init__(self):
#         self.plotter = None
#         self.lock = threading.Lock()

#     def ensure_plotter(self):
#         if self.plotter is None or self.plotter.app_window is None:
#             self.plotter = BackgroundPlotter(show=False)
#         return self.plotter

#     def update_plot(self, cfg):
#         with self.lock:
#             plotter = self.ensure_plotter()

#             mesh = pv.read(cfg["vtk_path"])
#             var = cfg["var"]
#             cmap = cfg["cmap"]
#             dim = cfg["dim"]

#             if var in ("u", "v", "w"):
#                 idx = {"u": 0, "v": 1, "w": 2}[var]
#                 data = mesh["velocity"][:, idx]
#             elif var == "velocity magnitude":
#                 data = np.linalg.norm(mesh["velocity"], axis=1)
#             elif var == "p":
#                 data = mesh["pressure"]
#             else:
#                 raise ValueError(f"Unknown variable '{var}'")

#             plotter.clear()
#             plotter.add_mesh(
#                 mesh,
#                 scalars=data,
#                 cmap=cmap,
#                 scalar_bar_args={"title": var},
#             )

#             if dim == "2":
#                 plotter.view_xy()
#                 plotter.enable_parallel_projection()
#                 plotter.reset_camera()

#             plotter.show()


# def socket_server(server):
#     try:
#         sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#         sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
#         sock.bind((HOST, PORT))
#         sock.listen(1)

#         print(f"[PlotServer] Listening on {HOST}:{PORT}", flush=True)

#         while True:
#             conn, _ = sock.accept()
#             try:
#                 msg = conn.recv(16384).decode()
#                 cfg = json.loads(msg)

#                 if cfg.get("cmd") == "exit":
#                     conn.close()
#                     break

#                 server.update_plot(cfg)

#             except Exception:
#                 traceback.print_exc()
#             finally:
#                 conn.close()

#         sock.close()
#         QtWidgets.QApplication.quit()

#     except Exception:
#         print("[PlotServer] Socket startup failed:")
#         traceback.print_exc()


# def main():
#     print("[PlotServer] Starting Qt application", flush=True)

#     # ✅ Correct Qt application creation
#     app = QtWidgets.QApplication.instance()
#     if app is None:
#         app = QtWidgets.QApplication(sys.argv)

#     server = PlotServer()

#     threading.Thread(
#         target=socket_server,
#         args=(server,),
#         daemon=True,
#     ).start()

#     app.exec_()


# if __name__ == "__main__":
#     main()



# # gui functions
# def ensure_plotter_server():
#     proc = STATE.get("plotter_process")
#     if proc is None or proc.poll() is not None:
#         script = os.path.join("src", "plotter.py")
#         proc = subprocess.Popen(
#             [sys.executable, script],
#             stdout=subprocess.PIPE,   # TEMP: capture for debugging
#             stderr=subprocess.PIPE,
#         )
#         STATE["plotter_process"] = proc

#     # ---- Wait until server is actually listening ----
#     for _ in range(20):  # ~2 seconds max
#         try:
#             sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#             sock.connect((HOST, PORT))
#             sock.close()
#             return
#         except ConnectionRefusedError:
#             time.sleep(0.1)

#     raise RuntimeError("Plotter server did not start")

# def update_plot():
#     try:
#         ensure_plotter_server()

#         cfg = {
#             "vtk_path": dpg.get_value("contour_vtk_path") or DEFAULT_VTK,
#             "var": dpg.get_value("contour_var"),
#             "cmap": dpg.get_value("contour_cmap") or "viridis",
#             "dim": dpg.get_value("param_domain_dimensions"),
#         }

#         sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#         sock.connect((HOST, PORT))
#         sock.sendall(json.dumps(cfg).encode())
#         sock.close()

#     except Exception as e:
#         append_log(f"Plotter connection error: {e}")

# def shutdown_plotter():
#     try:
#         sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#         sock.connect((HOST, PORT))
#         sock.sendall(json.dumps({"cmd": "exit"}).encode())
#         sock.close()
#     except:
#         pass

