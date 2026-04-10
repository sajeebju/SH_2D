#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def build_layered_model(nx, nz, dx, dz, layers):
    """
    Build a layered 2D model.

    Parameters
    ----------
    nx, nz : int
        Number of grid points in x and z
    dx, dz : float
        Grid spacing in x and z
    layers : list of tuples
        Each tuple is (thickness_m, vs_mps, rho_kgm3)

    Returns
    -------
    vs, rho : 2D numpy arrays of shape (nx, nz)
    """
    vs = np.zeros((nx, nz), dtype=np.float64)
    rho = np.zeros((nx, nz), dtype=np.float64)

    z_top = 0.0

    for thickness, vs_val, rho_val in layers:
        z_bot = z_top + thickness

        for j in range(nz):
            z = j * dz
            if z_top <= z < z_bot:
                vs[:, j] = vs_val
                rho[:, j] = rho_val

        z_top = z_bot

    # Fill any remaining depth with the last layer
    last_vs = layers[-1][1]
    last_rho = layers[-1][2]
    for j in range(nz):
        z = j * dz
        if z >= z_top:
            vs[:, j] = last_vs
            rho[:, j] = last_rho

    return vs, rho


def plot_model(model, title, outfile, xmax, zmax):
    plt.figure(figsize=(7, 4.5))
    plt.imshow(model.T, cmap="viridis", extent=[0, xmax, zmax, 0], aspect="auto")
    plt.colorbar(label=title)
    plt.xlabel("Distance (m)")
    plt.ylabel("Depth (m)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def main():
    # -----------------------------
    # Grid definition
    # -----------------------------
    xmax = 500.0
    zmax = 300.0
    nx = 201
    nz = 201

    dx = xmax / (nx - 1)
    dz = zmax / (nz - 1)

    # -----------------------------
    # Layered model definition
    # (thickness_m, Vs_m/s, rho_kg/m^3)
    # -----------------------------
    layers = [
        (200.0, 580.0, 1000.0),
        (100.0, 480.0, 1000.0),
    ]

    # -----------------------------
    # Build model
    # -----------------------------
    vs, rho = build_layered_model(nx, nz, dx, dz, layers)

    # -----------------------------
    # Save model files
    # -----------------------------
    np.save("vs_model.npy", vs)
    np.save("rho_model.npy", rho)

    # Optional text export
    np.savetxt("vs_model.txt", vs)
    np.savetxt("rho_model.txt", rho)

    # -----------------------------
    # Save preview plots
    # -----------------------------
    plot_model(vs, "Vs Model (m/s)", "vs_model.png", xmax, zmax)
    plot_model(rho, "Density Model (kg/m^3)", "rho_model.png", xmax, zmax)

    print("Model files written:")
    print("  vs_model.npy")
    print("  rho_model.npy")
    print("  vs_model.txt")
    print("  rho_model.txt")
    print("  vs_model.png")
    print("  rho_model.png")


if __name__ == "__main__":
    main()
