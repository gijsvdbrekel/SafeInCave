import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
import safeincave.PostProcessingTools as post
import safeincave as sf
import meshio
import json


hour = 60 * 60
MPa  = 1e6

GEOMETRY_TYPE = "irregular"  # or "regular"
PROBE_INDEX   = 0            # which auto-probe to use (0..N-1)




# =============== HULPFUNCTIES ===================


class WallProfileData:
    def __init__(self, operation_folder, scale=1.0):
        # Read displacements from xdmf
        points, self.time_list, u_field = post.read_node_vector(
            os.path.join(operation_folder, "u", "u.xdmf")
        )

        # Read msh node coordinates
        reader_msh = meshio.read(os.path.join(operation_folder, "mesh", "geom.msh"))
        points_msh = reader_msh.points

        # Find line elements = cavern wall
        wall_idx = None
        if hasattr(reader_msh, "cells_dict"):   # old meshio
            wall_idx = np.unique(reader_msh.cells_dict["line"].flatten())
        elif isinstance(reader_msh.cells, dict):  # dict-format
            if "line" in reader_msh.cells:
                wall_idx = np.unique(reader_msh.cells["line"].flatten())
        else:   # list of CellBlock
            for cb in reader_msh.cells:
                if getattr(cb, "type", None) == "line":
                    wall_idx = np.unique(cb.data.flatten())
                    break
        if wall_idx is None:
            raise ValueError("No 'line' cells found in mesh.")

        # Map msh indices to xdmf indices
        mapping = post.build_mapping(points_msh, points)
        for i in range(len(wall_idx)):
            wall_idx[i] = mapping[wall_idx[i]]

        self.wall_points = points[wall_idx]
        self.wall_u      = u_field[:, wall_idx]

        # sort by z
        sorted_idx = np.argsort(self.wall_points[:, 2])
        self.wall_points = self.wall_points[sorted_idx]
        self.wall_u      = self.wall_u[:, sorted_idx]
        self.scale       = scale


def auto_generate_probes(wall_profile: WallProfileData, n_bend_probes: int):
    """
    Same logic as in your other script:
      always top, mid, bottom + n_bend_probes extra 'bend' points.
    """
    pts = wall_profile.wall_points
    x = pts[:, 0]
    z = pts[:, 2]

    # base probes
    idx_bottom = np.argmin(z)
    idx_top    = np.argmax(z)
    z_mid      = 0.5 * (z[idx_bottom] + z[idx_top])
    idx_mid    = np.argmin(np.abs(z - z_mid))

    probes_idx = [idx_top, idx_mid, idx_bottom]

    # curvature-based bend points
    dx  = np.gradient(x)
    dz_ = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz_)

    denom = (dx*dx + dz_*dz_)**1.5
    denom[denom == 0.0] = np.inf
    curvature = np.abs(dx * ddz - dz_ * ddx) / denom
    curvature[np.isnan(curvature)] = 0.0

    idx_all = np.argsort(curvature)[::-1]

    def too_close(new_i, existing_idx, min_gap=5):
        return any(abs(new_i - ei) < min_gap for ei in existing_idx)

    for i in idx_all:
        if len(probes_idx) >= 3 + n_bend_probes:
            break
        if i in probes_idx:
            continue
        if too_close(i, probes_idx):
            continue
        probes_idx.append(i)

    # sort from top to bottom in z
    probes_idx = sorted(probes_idx, key=lambda k: z[k], reverse=True)
    return pts[probes_idx]


def compute_triaxial_stresses(p_Pa, q_Pa):
    """
    Neem aan: axisymmetrische toestand met sigma1 = axial, sigma3 = radial (sigma2 = sigma3).
    p = (sigma1 + 2 sigma3)/3
    q = |sigma1 - sigma3|  (von Mises)
    """
    p = p_Pa
    q = q_Pa
    sigma1 = p + 2.0/3.0 * q      # Pa
    sigma3 = p - 1.0/3.0 * q      # Pa
    return sigma1, sigma3

def read_probe_series(operation_folder, probe_point):
    """
    Lees alle relevante grootheden op een gegeven probe (dichtstbijzijnde cel).
    """
    # --- mean & von Mises stress ---
    pts, t_list, p_elems = post.read_cell_scalar(
        os.path.join(operation_folder, "p_elems", "p_elems.xdmf")
    )
    _,   _,     q_elems = post.read_cell_scalar(
        os.path.join(operation_folder, "q_elems", "q_elems.xdmf")
    )

    # index van cel dichtst bij probe
    idx = post.find_closest_point(probe_point, pts)

    p_t = p_elems[:, idx]       # Pa (negatief in compressie in SafeInCave)
    q_t = q_elems[:, idx]       # Pa (>=0)

    # maak compressie positief
    p_t_c = -p_t
    q_t_c = q_t

    # sigma_axial, sigma_radial (Pa)
    sigma_axial_Pa, sigma_radial_Pa = compute_triaxial_stresses(p_t_c, q_t_c)

    # naar MPa
    sigma_axial = sigma_axial_Pa / MPa
    sigma_radial = sigma_radial_Pa / MPa

    # --- totale rek tensorkomponent ---
    # aanname: eps_tot wordt als 3x3 tensor opgeslagen
    pts_eps, t_eps, eps_tot = post.read_cell_tensor(
        os.path.join(operation_folder, "eps_tot", "eps_tot.xdmf")
    )
    assert np.allclose(t_list, t_eps), "time_list mismatch between p_elems and eps_tot"

    # neem dezelfde cel-index
    eps_cell = eps_tot[:, idx, :, :]   # shape (Nt, 3, 3)

    # assen: z = axial, x = radial
    eps_axial  = -eps_cell[:, 2, 2] * 100.0
    eps_radial = -eps_cell[:, 0, 0] * 100.0
  

    # --- hardening-parameter en yield-functie ---
    _, _, alpha = post.read_cell_scalar(
        os.path.join(operation_folder, "alpha", "alpha.xdmf")
    )
    _, _, Fvp   = post.read_cell_scalar(
        os.path.join(operation_folder, "Fvp", "Fvp.xdmf")
    )

    alpha_t = alpha[:, idx]
    Fvp_t   = Fvp[:, idx]

    # --- invarianten voor p-q diagram ---
    # p (MPa, compressie > 0)
    p_MPa = p_t_c / MPa
    # J2: q is von Mises = sqrt(3 J2)  => sqrt(J2) = q / sqrt(3)
    sqrtJ2_MPa = q_t_c / (np.sqrt(3.0) * MPa)
    # I1 = 3 p  (MPa)
    I1_MPa = 3.0 * p_MPa

    t_hours = t_list / hour

    return {
        "t_hours":      t_hours,
        "sigma_axial":  sigma_axial,
        "sigma_radial": sigma_radial,
        "eps_axial":    eps_axial,
        "eps_radial":   eps_radial,
        "diff_stress":  sigma_axial - sigma_radial,
        "alpha":        alpha_t,
        "Fvp":          Fvp_t,
        "I1":           I1_MPa,
        "sqrtJ2":       sqrtJ2_MPa,
    }

def build_layout():
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 1.0, 0.15])

    ax_sig   = fig.add_subplot(gs[0, 0])
    ax_eps   = fig.add_subplot(gs[0, 1])
    ax_sig_eps = fig.add_subplot(gs[0, 2])

    ax_pq    = fig.add_subplot(gs[1, 0])
    ax_alpha = fig.add_subplot(gs[1, 1])
    ax_empty = fig.add_subplot(gs[1, 2])  # placeholder

    ax_slider = fig.add_subplot(gs[2, :])

    ax_empty.axis("off")

    return fig, ax_sig, ax_eps, ax_sig_eps, ax_pq, ax_alpha, ax_slider

# =============== HOOFDSCRIPT ===================

def main():
    global PROBE_INDEX

    output_folder    = os.path.join("output", "case_irregular(5)_100days_irregular_original")
    operation_folder = os.path.join(output_folder, "operation")

    # --- auto probes, same as in cavern plots ---
    wall_profile = WallProfileData(operation_folder, scale=1.0)

    if GEOMETRY_TYPE == "irregular":
        n_bend_probes = 3    # top + mid + bottom + 3 bends = 6
    else:  # "regular"
        n_bend_probes = 1    # top + mid + bottom + 1 bend = 4

    probes = auto_generate_probes(wall_profile, n_bend_probes=n_bend_probes)

    # pick which probe to analyse (0..len(probes)-1)
    PROBE_INDEX = min(PROBE_INDEX, probes.shape[0]-1)
    probe_point = probes[PROBE_INDEX]

    print("Using probe", PROBE_INDEX, "at", probe_point)

    series = read_probe_series(operation_folder, probe_point)

    t      = series["t_hours"]
    Nt     = t.size

    fig, ax_sig, ax_eps, ax_sig_eps, ax_pq, ax_alpha, ax_slider = build_layout()

    # --- 1. sigma_axial & sigma_radial vs tijd ---
    line_sig_ax, = ax_sig.plot(t, series["sigma_axial"], label=r"$\sigma_{\mathrm{axial}}$", color="tab:blue")
    line_sig_rd, = ax_sig.plot(t, series["sigma_radial"], label=r"$\sigma_{\mathrm{radial}}$", color="tab:red")
    marker_sig_ax = ax_sig.scatter(t[0], series["sigma_axial"][0], c="k", zorder=5)
    marker_sig_rd = ax_sig.scatter(t[0], series["sigma_radial"][0], c="k", zorder=5)
    ax_sig.set_xlabel("Time (h)")
    ax_sig.set_ylabel("Stress (MPa)")
    ax_sig.legend()

    # --- 2. eps_axial & eps_radial vs tijd ---
    line_eps_ax, = ax_eps.plot(t, series["eps_axial"], label=r"$\varepsilon_{\mathrm{axial}}$", color="tab:blue")
    line_eps_rd, = ax_eps.plot(t, series["eps_radial"], label=r"$\varepsilon_{\mathrm{radial}}$", color="tab:red")
    marker_eps_ax = ax_eps.scatter(t[0], series["eps_axial"][0], c="k", zorder=5)
    marker_eps_rd = ax_eps.scatter(t[0], series["eps_radial"][0], c="k", zorder=5)
    ax_eps.set_xlabel("Time (h)")
    ax_eps.set_ylabel("Total strain (%)")
    ax_eps.legend()

    # --- 3. differential stress vs axial strain ---
    line_sig_eps, = ax_sig_eps.plot(series["eps_axial"], series["diff_stress"], color="tab:blue")
    marker_sig_eps = ax_sig_eps.scatter(series["eps_axial"][0], series["diff_stress"][0], c="k", zorder=5)
    ax_sig_eps.set_xlabel("Axial strain (%)")
    ax_sig_eps.set_ylabel("Differential stress (MPa)")

    # --- 4. p-q diagram (I1 vs sqrt(J2)) ---
    line_pq, = ax_pq.plot(series["I1"], series["sqrtJ2"], color="tab:blue")
    marker_pq = ax_pq.scatter(series["I1"][0], series["sqrtJ2"][0], c="k", zorder=5)
    ax_pq.set_xlabel(r"$I_1$ (MPa)")
    ax_pq.set_ylabel(r"$\sqrt{J_2}$ (MPa)")
    ax_pq.set_title("Stress path in invariant space")

    # --- 5. alpha (en Fvp) vs tijd ---
    line_alpha, = ax_alpha.plot(t, series["alpha"], label="Hardening param.", color="tab:blue")
    ax_alpha.set_xlabel("Time (h)")
    ax_alpha.set_ylabel("Hardening param.")

    ax_alpha2 = ax_alpha.twinx()
    line_Fvp, = ax_alpha2.plot(t, series["Fvp"], label="Yield function", color="tab:red")
    ax_alpha2.set_ylabel("Yield function")

    # eenvoudige gecombineerde legend
    lines_alpha = [line_alpha, line_Fvp]
    labels_alpha = [l.get_label() for l in lines_alpha]
    ax_alpha.legend(lines_alpha, labels_alpha, loc="best")

    # --- slider ---
    slider = Slider(
        ax=ax_slider,
        label="Time (h)",
        valmin=t[0],
        valmax=t[-1],
        valinit=t[0],
        valstep=(t[1]-t[0]) if Nt > 1 else 1.0,
    )

    # mapping van tijd naar index
    def time_to_index(t_val):
        return int(np.clip(np.searchsorted(t, t_val), 0, Nt-1))

    def update(val):
        ti = time_to_index(slider.val)

        # update markers in alle figuren
        marker_sig_ax.set_offsets([t[ti], series["sigma_axial"][ti]])
        marker_sig_rd.set_offsets([t[ti], series["sigma_radial"][ti]])

        marker_eps_ax.set_offsets([t[ti], series["eps_axial"][ti]])
        marker_eps_rd.set_offsets([t[ti], series["eps_radial"][ti]])

        marker_sig_eps.set_offsets([series["eps_axial"][ti], series["diff_stress"][ti]])

        marker_pq.set_offsets([series["I1"][ti], series["sqrtJ2"][ti]])

        # alpha & Fvp kun je ook markeren als je wilt; hier laten we alleen de lijnen.
        fig.canvas.draw_idle()

    slider.on_changed(update)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
