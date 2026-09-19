#!/usr/bin/env python3
"""
Render smooth-shaded animations of velocity magnitude and pressure from a
sequence of MeMPhyS Solution_NNNNNN.vtk files (legacy ASCII VTK, written by
write_vtk() in src/c_header_files/write_functions.c).

No pyvista/vtk dependency — parses the exact format MeMPhyS writes
(POINTS / CELLS / CELL_TYPES / POINT_DATA with VECTORS velocity +
SCALARS pressure) directly with numpy. Produces two separate animation
files (velocity, pressure) with continuous (Gouraud-shaded) coloring —
no contour banding.

Usage:
    python3 make_animation.py
    python3 make_animation.py --dir . --pattern "Solution_*.vtk" \
        --out-velocity velocity.mp4 --out-pressure pressure.mp4 \
        --fps 5 --stride 1 --pclip 5 95

Output is .mp4 if ffmpeg is available, otherwise falls back to .gif
automatically.
"""
import argparse
import glob
import os
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter


def read_n_values(f, n, dtype=float):
    vals = []
    while len(vals) < n:
        line = f.readline()
        if not line:
            raise EOFError("Unexpected end of file while reading values")
        vals.extend(line.split())
    return np.array(vals[:n], dtype=dtype)


def read_geometry(path):
    """Read POINTS + triangle CELLS once — geometry is constant across frames."""
    with open(path, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("POINTS"):
                n = int(line.split()[1])
                pts = read_n_values(f, n * 3).reshape(n, 3)
            elif line.startswith("CELLS"):
                parts = line.split()
                ncells, total = int(parts[1]), int(parts[2])
                raw = read_n_values(f, total, dtype=int)
                nper = raw[0]
                cells = raw.reshape(ncells, nper + 1)[:, 1:]
                return pts, cells, n
            line = f.readline()
    raise ValueError(f"Could not find POINTS/CELLS in {path}")


def read_fields(path, num_points):
    """Read velocity (u, v) + pressure for one frame. Skips geometry."""
    with open(path, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("VECTORS velocity"):
                vec = read_n_values(f, num_points * 3).reshape(num_points, 3)
                u, v = vec[:, 0], vec[:, 1]
            elif line.startswith("SCALARS pressure"):
                f.readline()  # LOOKUP_TABLE line
                p = read_n_values(f, num_points)
                return u, v, p
            line = f.readline()
    raise ValueError(f"Could not find VECTORS velocity / SCALARS pressure in {path}")


def find_ffmpeg():
    """Locate a WORKING ffmpeg binary: system PATH first — but actually
    test-run it, since a module can put a broken binary in PATH (e.g.
    missing a shared library) and `shutil.which` alone won't catch that,
    it only checks the file exists. Falls back to the imageio-ffmpeg
    bundled binary (`pip install --user imageio-ffmpeg`) if the system
    one is missing or fails to run. Returns None if neither works."""
    import subprocess
    path = shutil.which("ffmpeg")
    if path:
        try:
            subprocess.run([path, "-version"], capture_output=True, check=True, timeout=10)
            return path
        except Exception as e:
            print(f"System ffmpeg found at {path} but failed to run ({e}) "
                  f"— trying imageio-ffmpeg fallback")
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        return None


def make_writer(out_path, fps):
    """Pick ffmpeg (mp4) if available (system or imageio-ffmpeg bundled),
    else fall back to Pillow (gif)."""
    ffmpeg_path = find_ffmpeg()
    if ffmpeg_path:
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg_path
        if not out_path.endswith(".mp4"):
            out_path = os.path.splitext(out_path)[0] + ".mp4"
        return out_path, FFMpegWriter(fps=fps)
    else:
        if not out_path.endswith(".gif"):
            out_path = os.path.splitext(out_path)[0] + ".gif"
        return out_path, PillowWriter(fps=fps)


def render(tri, frames, values_per_frame, vmin, vmax, cmap, label, out_path, fps, files):
    fig, ax = plt.subplots(figsize=(8, 6))
    tpc = ax.tripcolor(tri, values_per_frame[0], shading="gouraud",
                        cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_aspect("equal")
    ax.set_title("")
    fig.colorbar(tpc, ax=ax, label=label)

    def update(i):
        tpc.set_array(values_per_frame[i])
        step = os.path.basename(files[i])
        ax.set_title(f"{label} — {step}  ({i+1}/{frames})")
        print(f"  [{label}] rendering frame {i+1}/{frames}: {step}")
        return (tpc,)

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    out_path, writer = make_writer(out_path, fps)
    anim.save(out_path, writer=writer)
    plt.close(fig)
    print(f"Saved {label} animation to {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".", help="directory containing the VTK files")
    ap.add_argument("--pattern", default="Solution_*.vtk", help="glob pattern for VTK files")
    ap.add_argument("--out-velocity", default="velocity_animation.mp4")
    ap.add_argument("--out-pressure", default="pressure_animation.mp4")
    ap.add_argument("--fps", type=int, default=5, help="frames per second")
    ap.add_argument("--stride", type=int, default=1, help="use every Nth file (for quick previews)")
    ap.add_argument("--pclip", type=float, nargs=2, default=(5, 95),
                     help="percentile clipping for color limits (tighter = shorter, punchier range)")
    ap.add_argument("--cmap-velocity", default="RdBu_r", help="colormap for velocity (default: blue-white-red diverging)")
    ap.add_argument("--cmap-pressure", default="RdBu_r", help="colormap for pressure (default: blue-white-red diverging)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, args.pattern)))[::args.stride]
    if not files:
        sys.exit(f"No files matched {args.pattern!r} in {args.dir!r}")
    print(f"Found {len(files)} VTK files (stride={args.stride})")

    print("Reading geometry from first frame...")
    points, cells, num_points = read_geometry(files[0])
    x, y = points[:, 0], points[:, 1]
    tri = mtri.Triangulation(x, y, triangles=cells)

    # Single pass: read every frame's fields once, cache for both animations.
    print("Reading fields for all frames...")
    vmag_frames, p_frames = [], []
    for i, fpath in enumerate(files):
        u, v, p = read_fields(fpath, num_points)
        vmag_frames.append(np.hypot(u, v))
        p_frames.append(p)
        if i % max(1, len(files) // 10) == 0:
            print(f"  {i+1}/{len(files)}")

    vmag_all = np.concatenate(vmag_frames)
    p_all = np.concatenate(p_frames)

    lo, hi = args.pclip
    vmag_min, vmag_max = np.nanpercentile(vmag_all, [lo, hi])
    p_min, p_max = np.nanpercentile(p_all, [lo, hi])
    print(f"Velocity magnitude color range: [{vmag_min:.4g}, {vmag_max:.4g}] ({lo}-{hi} percentile)")
    print(f"Pressure color range: [{p_min:.4g}, {p_max:.4g}] ({lo}-{hi} percentile)")

    render(tri, len(files), vmag_frames, vmag_min, vmag_max,
           args.cmap_velocity, "Velocity magnitude", args.out_velocity, args.fps, files)

    render(tri, len(files), p_frames, p_min, p_max,
           args.cmap_pressure, "Pressure", args.out_pressure, args.fps, files)


if __name__ == "__main__":
    main()
