#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sh2dwave_wrapper


def get_env_float(name, default):
    return float(os.getenv(name, default))


def get_env_int(name, default):
    return int(float(os.getenv(name, default)))


def get_env_str(name, default):
    return os.getenv(name, default)


def load_model_file(path, nx, nz):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Model file not found: {path}")

    ext = os.path.splitext(path)[1].lower()

    if ext == ".npy":
        model = np.load(path)

    elif ext in [".txt", ".dat", ".csv"]:
        model = np.loadtxt(path)

    else:
        raise ValueError(f"Unsupported file format: {ext}")

    # Ensure correct shape
    if model.shape != (nx, nz):
        raise ValueError(
            f"Model shape {model.shape} does not match (nx, nz)=({nx}, {nz})"
        )

    return model.astype(np.float64)


# -----------------------------
# Read runtime parameters
# -----------------------------
mode = get_env_str("MODE", "wave").strip().lower()

xmax = get_env_float("XMAX", 500.0)
zmax = get_env_float("ZMAX", 300.0)

nx = get_env_int("NX", 201)
nz = get_env_int("NZ", 201)

dx = xmax / (nx - 1)
dz = zmax / (nz - 1)

xsrc = get_env_float("XSRC", 250.0)
zsrc = get_env_float("ZSRC", 5.0)

xr = get_env_float("XR", 330.0)
zr = get_env_float("ZR", 150.0)

cfl = get_env_float("CFL", 0.5)
tmax = get_env_float("TMAX", 1.50)

f0 = get_env_float("F0", 40.0)
t0 = get_env_float("T0", 4.0 / f0)

accuracy = get_env_int("ACCURACY", 2)

use_absorb = get_env_int("USE_ABSORB", 0)
w = get_env_int("W", 60)
a = get_env_float("A", 0.0053)

mskip = get_env_int("MSKIP", 5)
fps = get_env_int("FPS", 20)
interval = get_env_int("INTERVAL", 50)

out_movie = get_env_str("OUT_MOVIE", "cpp_wave_animation.mp4")
out_seis = get_env_str("OUT_SEIS", "seismogram.txt")

vs_file = get_env_str("VS_FILE", "vs_model.npy")
rho_file = get_env_str("RHO_FILE", "rho_model.npy")



# -----------------------------
# Read external model files
# -----------------------------
if not os.path.exists(vs_file):
    raise FileNotFoundError(f"Vs model file not found: {vs_file}")

if not os.path.exists(rho_file):
    raise FileNotFoundError(f"Rho model file not found: {rho_file}")

vs = load_model_file(vs_file, nx, nz)
rho = load_model_file(rho_file, nx, nz)

if vs.shape != (nx, nz):
    raise ValueError(f"Vs model shape {vs.shape} does not match (nx, nz)=({nx}, {nz})")

if rho.shape != (nx, nz):
    raise ValueError(f"Rho model shape {rho.shape} does not match (nx, nz)=({nx}, {nz})")



# -----------------------------
# Time setup
# -----------------------------
dt = cfl * dx / np.max(vs)
nt = int(tmax / dt)

print(f"MODE={mode}")
print(f"nx={nx}, nz={nz}, nt={nt}")
print(f"dx={dx:.6f}, dz={dz:.6f}, dt={dt:.6f}")
print(f"Vs file = {vs_file}")
print(f"Rho file = {rho_file}")
print(f"use_absorb={use_absorb}, w={w}, a={a}, accuracy={accuracy}")




# -----------------------------
# MODE = seis
# -----------------------------
if mode == "seis":
    isrc = int(xsrc / dx)
    jsrc = int(zsrc / dz)
    ir = int(xr / dx)
    jr = int(zr / dz)

    if not (0 <= isrc < nx and 0 <= jsrc < nz):
        raise ValueError("Source is outside the model")

    if not (0 <= ir < nx and 0 <= jr < nz):
        raise ValueError("Receiver is outside the model")

    time, seis = sh2dwave_wrapper.py_sh_seis(
        nx, nz, nt, dx, dz, dt,
        f0, t0, isrc, jsrc, ir, jr,
        vs, rho,
        use_absorb, w, a,
        accuracy
    )

    np.savetxt(out_seis, np.column_stack((time, seis)), fmt="%.8e")
    print(f"Seismogram saved to {out_seis}")

    plt.figure(figsize=(8, 4))
    plt.plot(time, seis, linewidth=1.2)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.title("Synthetic Seismogram")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# -----------------------------
# MODE = wave
# -----------------------------
elif mode == "wave":
    total_solution_space, vy = sh2dwave_wrapper.py_wave_propagate(
        nx, nz, nt, dx, dz, dt,
        f0, t0, int(xsrc), int(zsrc),
        vs, rho,
        use_absorb, w, a,
        accuracy
    )

    total_solution_space = np.asarray(total_solution_space, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)

    max_amp = np.max(np.abs(total_solution_space))
    if max_amp == 0.0:
        max_amp = 1.0

    vscale = vy.T - np.mean(vy)
    vden = np.max(np.abs(vscale))
    if vden == 0.0:
        vden = 1.0
    vscale = -vscale / vden

    fig, ax1 = plt.subplots(figsize=(8, 5))

    def SH_2D(frame):
        k = frame * mskip
        if k >= nt:
            k = nt - 1

        ax1.clear()

        wave_field = total_solution_space[:, :, k].T / max_amp
        data = wave_field + 0.0025 * vscale[:, ::-1]

        ax1.imshow(
            vs.T,
            cmap="winter",
            extent=[0, xmax, zmax, 0],
            alpha=0.7
        )

        im = ax1.imshow(
            data,
            cmap="seismic",
            extent=[0, xmax, zmax, 0],
            vmin=-0.08,
            vmax=0.08,
            alpha=0.6
        )

        ax1.set_xlim([0, xmax])
        ax1.set_ylim([zmax, 0])
        ax1.set_xlabel("Distance (m)", fontsize=12)
        ax1.set_ylabel("Depth (m)", fontsize=12)
        ax1.set_title(f"2D SH at {k * dt * 1000:.2f} ms")

        return [im]

    nframes = max(1, nt // mskip)

    anim = animation.FuncAnimation(
        fig,
        SH_2D,
        frames=nframes,
        interval=interval,
        blit=False
    )

    anim.save(out_movie, writer="ffmpeg", fps=fps)
    print(f"Wavefield movie saved to {out_movie}")
    plt.show()

else:
    raise ValueError("MODE must be either 'wave' or 'seis'")