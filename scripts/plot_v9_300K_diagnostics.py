#!/usr/bin/env python3
"""Generate provenance-audited canonical v3 300 K diagnostic figures."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from arrhenius_fracture.v9_canonical_four_class_elastic_fem import (
    CanonicalElasticFourClassFEMCondition,
)
from scripts.run_v9_quiet_tail_kernel import CLASSES


CLASS_ORDER = ("Peak", "DBTT", "weak-T", "ceramic")
COLORS = {"Peak": "#c43c39", "DBTT": "#3569a8", "weak-T": "#2d8a57", "ceramic": "#8a5aa5"}
MARKERS = {"Peak": "o", "DBTT": "s", "weak-T": "^", "ceramic": "D"}
QUANTILES = (("N10", 0.10536051565782628), ("N50", 0.6931471805599453),
             ("N90", 2.302585092994046))


def csv_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def numeric_key(row):
    return row["material_class"], float(row["sigma_a_MPa"]), float(row["N"])


def active_generation(path):
    generation = json.loads((path / "ACTIVE.json").read_text())["generation"]
    return path / generation


def load_checkpoint(path):
    generation = active_generation(path)
    summary = json.loads((generation / "summary.json").read_text())
    metadata = json.loads((generation / "state_metadata.json").read_text())
    with np.load(generation / "state_arrays.npz", allow_pickle=False) as archive:
        arrays = {key: np.asarray(archive[key]).copy() for key in archive.files}
    return summary, metadata, arrays, generation.name


def checkpoint_path(root, material, stress, N, anchors):
    if math.isclose(stress, anchors[material], rel_tol=0.0, abs_tol=1e-12):
        base = root / "quiet_tail_checkpoints" / material
    else:
        base = root / "quiet_tail_checkpoints" / material / "descent" / f"sigmaA_{stress:g}MPa"
    return base / f"N_{N:.17g}"


def audit_tables(root):
    quiet = csv_rows(root / "300K_quiet_tail_state.csv")
    survival = csv_rows(root / "300K_stable_birth_survival.csv")
    descent = csv_rows(root / "300K_endurance_stress_descent.csv")
    gates = csv_rows(root / "300K_energy_gate_envelope.csv")
    Q = {numeric_key(row): row for row in quiet}
    S = {numeric_key(row): row for row in survival}
    D = {numeric_key(row): row for row in descent}
    if len(Q) != len(quiet) or len(S) != len(survival) or len(D) != len(descent):
        raise RuntimeError("duplicate numerical condition keys remain in v3 tables")
    if set(Q) != set(S):
        raise RuntimeError("quiet-tail and survival physical checkpoint keys differ")
    mismatch = [key for key in Q if float(Q[key]["H_cleave"]) != float(S[key]["H_cleave"])]
    if mismatch:
        raise RuntimeError(f"survival table mixes projected and direct H at {mismatch}")
    descent_mismatch = [key for key in D if float(D[key]["H_cleave"]) != float(Q[key]["H_cleave"])]
    if descent_mismatch:
        raise RuntimeError(f"descent table differs from direct checkpoint H at {descent_mismatch}")
    anchors = {name: stress for name, (_option, stress) in CLASSES.items()}
    checkpoint_records = {}
    for key, row in Q.items():
        material, stress, N = key
        path = checkpoint_path(root, material, stress, N, anchors)
        summary, metadata, arrays, generation = load_checkpoint(path)
        direct_H = float(row["H_cleave"])
        capsule_H = float(metadata["condition_capsule"]["birth"]["cumulative_cleavage_hazard"])
        if float(summary["H_cleave"]) != direct_H or capsule_H != direct_H:
            raise RuntimeError(f"atomic checkpoint H mismatch at {key}")
        if float(summary["cycle_hazard"]) != float(row["cycle_hazard"]):
            raise RuntimeError(f"atomic cycle hazard mismatch at {key}")
        checkpoint_records[key] = {
            "path": path, "summary": summary, "metadata": metadata,
            "arrays": arrays, "generation": generation,
        }
    return quiet, survival, descent, gates, checkpoint_records


def save_figure(fig, output, stem):
    fig.savefig(output / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def style_axes(ax):
    ax.grid(True, which="both", alpha=0.22, linewidth=0.6)
    ax.tick_params(labelsize=8)


def direct_quantiles(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["material_class"], float(row["sigma_a_MPa"])), []).append(row)
    result = []
    for (material, stress), group in groups.items():
        values = sorted((float(row["N"]), float(row["H_cleave"])) for row in group)
        for label, threshold in QUANTILES:
            for (n0, h0), (n1, h1) in zip(values, values[1:]):
                if h0 <= threshold <= h1 and h0 > 0 and h1 > h0:
                    fraction = (math.log(threshold)-math.log(h0))/(math.log(h1)-math.log(h0))
                    N = math.exp(math.log(n0)+fraction*(math.log(n1)-math.log(n0)))
                    result.append((material, stress, label, N, n0, n1))
                    break
    return result


def plot_sn(root, output, survival):
    estimates = csv_rows(root / "300K_practical_survival_stress_estimates.csv")
    quantiles = direct_quantiles(survival)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        color = COLORS[material]
        for label, marker in (("N10", "o"), ("N50", "s"), ("N90", "^")):
            points = [(N, stress) for mat, stress, q, N, _n0, _n1 in quantiles
                      if mat == material and q == label]
            if points:
                ax.scatter(*zip(*points), marker=marker, s=36, color=color,
                           edgecolor="black", linewidth=0.45, zorder=3)
        for row in estimates:
            if row["material_class"] != material:
                continue
            label = {0.9: "N10", 0.5: "N50", 0.1: "N90"}[
                float(row["target_S_stable_birth_1e14"])]
            marker = {"N10": "o", "N50": "s", "N90": "^"}[label]
            ax.scatter(1e14, float(row["estimated_sigma_a_MPa"]), marker=marker,
                       s=55, facecolor="none", edgecolor=color, linewidth=1.5, zorder=4)
        ax.set_xscale("log"); style_axes(ax); ax.set_title(material, fontsize=10, weight="bold")
        ax.set_ylabel("stress amplitude, MPa", fontsize=9)
    for ax in axes[-1]: ax.set_xlabel("cycles to stable birth / finite-horizon N", fontsize=9)
    legend = [
        Line2D([], [], marker="o", linestyle="none", color="black", label="N10"),
        Line2D([], [], marker="s", linestyle="none", color="black", label="N50"),
        Line2D([], [], marker="^", linestyle="none", color="black", label="N90"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
               color="0.25", label=r"interpolated at $10^{14}$"),
    ]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=4, frameon=False, fontsize=8)
    fig.suptitle("Canonical v3, 300 K: probabilistic S–N diagnostics", y=0.995, fontsize=12)
    fig.subplots_adjust(top=0.86)
    save_figure(fig, output, "01_probabilistic_SN_300K")


def descent_groups(descent):
    groups = {}
    for row in descent:
        groups.setdefault((row["material_class"], float(row["sigma_a_MPa"])), []).append(row)
    return {key: sorted(rows, key=lambda row: float(row["N"])) for key, rows in groups.items()}


def plot_hazards(output, descent, quiet):
    groups = descent_groups(descent)
    quiet_map = {numeric_key(row): row for row in quiet}
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        material_groups = [(key, rows) for key, rows in sorted(groups.items())
                           if key[0] == material]
        for index, ((_mat, stress), rows) in enumerate(material_groups):
            N = np.asarray([float(row["N"]) for row in rows])
            H = np.asarray([float(row["H_cleave"]) for row in rows])
            ax.plot(N, H, marker="o", label=f"{stress:g} MPa",
                    color=COLORS[material], alpha=0.55+0.4*index, linewidth=1.5)
        ax.set_xscale("log"); ax.set_yscale("log"); style_axes(ax)
        ax.set_title(material, fontsize=10, weight="bold"); ax.legend(fontsize=8, frameon=False)
        ax.set_ylabel(r"direct $H_{\rm cleave}(N)$", fontsize=9)
    for ax in axes[-1]: ax.set_xlabel("cycles N", fontsize=9)
    fig.suptitle("Cumulative cleavage action from actual v3 checkpoints", fontsize=12)
    save_figure(fig, output, "02_cumulative_cleavage_hazard")

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        material_groups = [(key, rows) for key, rows in sorted(groups.items())
                           if key[0] == material]
        for index, ((_mat, stress), rows) in enumerate(material_groups):
            N = np.asarray([float(row["N"]) for row in rows])
            h = np.asarray([float(quiet_map[(material, stress, float(row["N"]))]["cycle_hazard"])
                            for row in rows])
            ax.plot(N, h, marker="o", label=f"{stress:g} MPa",
                    color=COLORS[material], alpha=0.55+0.4*index, linewidth=1.5)
        ax.set_xscale("log"); ax.set_yscale("log"); style_axes(ax)
        ax.set_title(material, fontsize=10, weight="bold"); ax.legend(fontsize=8, frameon=False)
        ax.set_ylabel("direct instantaneous cycle hazard", fontsize=9)
    for ax in axes[-1]: ax.set_xlabel("cycles N", fontsize=9)
    fig.suptitle("Instantaneous cleavage hazard (not H/N)", fontsize=12)
    save_figure(fig, output, "03_instantaneous_cycle_hazard")


def array_sum(arrays, name):
    return float(np.sum(np.asarray(arrays[f"kernel_{name}"], float)))


def plot_mpz(output, records):
    selected = {"Peak": 631.0, "DBTT": 648.0, "weak-T": 310.0, "ceramic": 332.0}
    data = {}
    for material, stress in selected.items():
        rows = []
        for N in (1e8, 1e10, 1e12, 1e14):
            record = records[(material, stress, N)]
            a = record["arrays"]; m = record["metadata"]
            scalars = m["condition_capsule"]["birth"]["mpz"]["scalars"]
            rows.append({
                "N": N,
                "backstress": float(np.max(np.abs(a["kernel_sigma_back_by_system_Pa"]))),
                "shielding": float(m["kernel_scalars"]["signed_active_K_shield_Pa_sqrt_m"]),
                "tip_radius": float(m["kernel_scalars"]["tip_radius_m"]),
                "mobile_positive": array_sum(a, "mobile_positive"),
                "mobile_negative": array_sum(a, "mobile_negative"),
                "retained_positive": array_sum(a, "retained_positive"),
                "retained_negative": array_sum(a, "retained_negative"),
                "slip": array_sum(a, "accumulated_slip_positive") + array_sum(a, "accumulated_slip_negative"),
                "wake": sum(array_sum(a, name) for name in (
                    "wake_mobile_positive", "wake_mobile_negative",
                    "wake_retained_positive", "wake_retained_negative")),
                "advance": float(scalars["advance_total_m"]),
            })
        data[material] = rows
    panels = [
        ("backstress", "Taylor backstress, GPa", 1e-9, "log"),
        ("shielding", r"$|$signed shielding$|$, MPa$\sqrt{m}$", 1e-6, "log"),
        ("tip_radius", r"tip radius, $\mu$m", 1e6, "linear"),
        ("mobile", "signed-channel mobile content", 1.0, "log"),
        ("retained", "signed-channel retained content", 1.0, "log"),
        ("slip", "accumulated slip", 1.0, "log"),
        ("wake", "wake mobile + retained", 1.0, "linear"),
        ("advance", r"advection / advance, $\mu$m", 1e6, "linear"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(10, 11), sharex=True)
    for ax, (key, ylabel, scale, yscale) in zip(axes.flat, panels):
        for material, rows in data.items():
            N = [row["N"] for row in rows]; color = COLORS[material]
            if key == "mobile":
                ax.plot(N, [row["mobile_positive"] for row in rows], color=color, marker="o", label=material)
                ax.plot(N, [row["mobile_negative"] for row in rows], color=color, marker="x", linestyle="--")
            elif key == "retained":
                ax.plot(N, [row["retained_positive"] for row in rows], color=color, marker="o", label=material)
                ax.plot(N, [row["retained_negative"] for row in rows], color=color, marker="x", linestyle="--")
            elif key == "shielding":
                ax.plot(N, [abs(row[key])*scale for row in rows], color=color,
                        marker="o", label=material)
            else:
                ax.plot(N, [row[key]*scale for row in rows], color=color, marker="o", label=material)
        ax.set_xscale("log")
        if yscale == "log":
            ax.set_yscale(yscale)
        ax.set_ylabel(ylabel, fontsize=8); style_axes(ax)
        if key == "shielding":
            ax.text(0.02, 0.04, "all values negative (direct shielding sign)",
                    transform=ax.transAxes, fontsize=7)
        if key in ("wake", "advance"):
            ax.text(0.5, 0.5, "identically zero at all four checkpoints",
                    transform=ax.transAxes, ha="center", va="center", fontsize=8)
    for ax in axes[-1]: ax.set_xlabel("cycles N", fontsize=9)
    axes[0, 0].legend(fontsize=8, ncol=2, frameon=False)
    axes[1, 1].text(0.02, 0.04, "+ solid circles; − dashed crosses",
                    transform=axes[1, 1].transAxes, fontsize=7)
    axes[2, 0].text(0.02, 0.04, "+ solid circles; − dashed crosses",
                    transform=axes[2, 0].transAxes, fontsize=7)
    fig.suptitle("Direct signed-MPZ state evolution at lower survival brackets", fontsize=12)
    save_figure(fig, output, "04_physical_MPZ_state_evolution")


def plot_phase_response(output, records):
    selected = {"Peak": 631.0, "DBTT": 648.0, "weak-T": 310.0, "ceramic": 332.0}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
    cmap = plt.get_cmap("viridis")
    checkpoints = (1e8, 1e10, 1e12, 1e14)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        stress = selected[material]; twin = ax.twinx()
        for index, N in enumerate(checkpoints):
            arrays = records[(material, stress, N)]["arrays"]
            phase = arrays["kernel_phase"]
            ax.plot(phase, arrays["kernel_cleavage_log_rate_effective_s"],
                    color=cmap(index/3), linewidth=1.35, label=f"$10^{{{int(math.log10(N))}}}$")
        arrays = records[(material, stress, checkpoints[-1])]["arrays"]
        twin.plot(arrays["kernel_phase"], arrays["kernel_root_opening_Pa"]*1e-9,
                  color="0.25", linestyle="--", linewidth=1.0)
        ax.set_title(f"{material}, {stress:g} MPa", fontsize=9, weight="bold")
        ax.set_ylabel("effective cleavage log-rate, s$^{-1}$", fontsize=8)
        twin.set_ylabel("root opening, GPa", fontsize=8, color="0.25")
        style_axes(ax); twin.tick_params(labelsize=7, colors="0.25")
    for ax in axes[-1]: ax.set_xlabel(r"cycle phase $\phi$", fontsize=9)
    axes[0, 0].legend(fontsize=7, ncol=2, frameon=False)
    fig.suptitle("Phase-resolved cleavage response from atomic checkpoints", fontsize=12)
    save_figure(fig, output, "05_phase_resolved_cleavage_response")


def plot_gate(output, descent):
    groups = descent_groups(descent)
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 8.2), sharex=True)
    for (material, stress), rows in sorted(groups.items()):
        N = np.asarray([float(row["N"]) for row in rows]); color = COLORS[material]
        label = f"{material} {stress:g} MPa"
        axes[0].plot(N, [float(row["conditional_attempt_admission_probability"]) for row in rows],
                     color=color, marker=MARKERS[material], label=label, alpha=0.75)
        axes[0].plot(N, [float(row["rejected_phase_xi_fraction"]) for row in rows],
                     color=color, linestyle="--", alpha=0.75)
        axes[1].fill_between(N,
            [float(row["minimum_admitted_length_m"])*1e6 for row in rows],
            [float(row["maximum_admitted_length_m"])*1e6 for row in rows],
            color=color, alpha=0.10)
        axes[1].plot(N, [float(row["maximum_admitted_length_m"])*1e6 for row in rows],
                     color=color, linewidth=1.0, alpha=0.75)
        axes[2].plot(N, [float(row["minimum_energy_gate_margin_J_per_m"]) for row in rows],
                     color=color, marker=".", alpha=0.75)
    axes[0].set_ylabel("admission probability\n(solid) / rejected fraction (dash)", fontsize=8)
    axes[1].set_ylabel("admitted length range, μm", fontsize=8)
    axes[2].set_ylabel("minimum trial energy\nresidual, J m$^{-1}$", fontsize=8)
    axes[2].set_xlabel("cycles N", fontsize=9); axes[2].axhline(0, color="black", linewidth=0.7)
    for ax in axes: ax.set_xscale("log"); style_axes(ax)
    axes[0].legend(fontsize=7, ncol=2, frameon=False)
    fig.suptitle("Post-first-passage energy-admission gate (not PD damage)", fontsize=12)
    save_figure(fig, output, "06_energy_gate_diagnostics")


def plot_geometry(output, run_args, source_root):
    option, stress = CLASSES["Peak"]
    condition = CanonicalElasticFourClassFEMCondition.from_run_args(
        run_args, option, source_root, stress, 1720)
    mesh = condition.mesh; root = condition.root_xy
    triangles = mesh.elems
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax in axes:
        ax.triplot(mesh.nodes[:, 0]*1e3, mesh.nodes[:, 1]*1e3, triangles,
                   color="#54606c", linewidth=0.22, alpha=0.65)
        ax.set_aspect("equal"); ax.set_xlabel("x, mm"); ax.set_ylabel("y, mm")
    axes[0].set_title("actual K360 elastic FEM mesh", fontsize=10, weight="bold")
    zoom_x = (root[0]-0.08e-3, root[0]+0.28e-3)
    zoom_y = (root[1]-0.20e-3, root[1]+0.20e-3)
    axes[1].set_xlim(np.asarray(zoom_x)*1e3); axes[1].set_ylim(np.asarray(zoom_y)*1e3)
    axes[1].set_title("notch-root refinement", fontsize=10, weight="bold")
    axes[1].scatter([root[0]*1e3], [root[1]*1e3], color="#c43c39", s=28, zorder=5)
    args = json.loads(Path(run_args).read_text())
    corridor = plt.Rectangle((root[0]*1e3, (root[1]-args["path_refine_half_height_m"])*1e3),
                             args["path_refine_length_m"]*1e3,
                             2*args["path_refine_half_height_m"]*1e3,
                             fill=False, edgecolor="#e38b2c", linewidth=1.3, linestyle="--")
    axes[1].add_patch(corridor)
    length_m = float(condition.birth.mpz.state.cfg.length_m)
    axes[1].annotate("root-local signed MPZ", xy=(root[0]*1e3, root[1]*1e3),
                     xytext=((root[0]+0.02e-3)*1e3, (root[1]-0.15e-3)*1e3),
                     arrowprops=dict(arrowstyle="->", color="#2d8a57"), fontsize=8)
    axes[1].arrow(root[0]*1e3, root[1]*1e3, length_m*1e3, 0,
                  width=0.001, head_width=0.018, head_length=0.012,
                  color="#2d8a57", length_includes_head=True)
    axes[1].annotate("prospective +x crack direction",
                     xy=((root[0]+0.12e-3)*1e3, root[1]*1e3),
                     xytext=((root[0]+0.12e-3)*1e3, (root[1]+0.075e-3)*1e3),
                     fontsize=8, color="#2d8a57")
    annotation = ("elliptical edge notch\n"
                  r"$a=360\,\mu$m; nominal $\rho=45\,\mu$m" "\n"
                  r"$b=127.279\,\mu$m" "\n"
                  f"mesh-resolved root radius = {condition.root_radius_initial_m*1e6:.3f} μm\n"
                  "orange dashed: local refinement corridor")
    axes[1].text(0.02, 0.98, annotation, transform=axes[1].transAxes,
                 va="top", fontsize=8, bbox=dict(facecolor="white", alpha=0.88, edgecolor="0.8"))
    fig.suptitle("Canonical v3 initial K360 geometry and root-local constitutive region", fontsize=12)
    save_figure(fig, output, "07_initial_geometry_mesh")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--run-args", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False})
    quiet, survival, descent, gates, records = audit_tables(args.run_root)
    plot_sn(args.run_root, args.output, survival)
    plot_hazards(args.output, descent, quiet)
    plot_mpz(args.output, records)
    plot_phase_response(args.output, records)
    plot_gate(args.output, descent)
    plot_geometry(args.output, args.run_args, args.source_root)
    audit = {
        "schema": "V9_CANONICAL_300K_FIGURE_DATA_AUDIT_1",
        "commit_expected": "5d6b1eef331ca88a4a55c8a34fc7d988f2882ac9",
        "direct_checkpoint_conditions": len(quiet),
        "descent_checkpoint_conditions": len(descent),
        "gate_rows": len(gates),
        "atomic_checkpoints_verified": len(records),
        "H_relationship": "quiet_tail_state == atomic summary == atomic birth capsule == corrected survival",
        "cycle_hazard_semantics": "direct phase-integrated checkpoint kernel quantity; never H/N",
        "survival_semantics": "derived from direct checkpoint H and evaluated gate admission",
        "practical_stress_semantics": "bracketed finite_1e14 interpolation, not endurance limit",
        "stationary_projection_in_physical_rows": False,
    }
    (args.output / "figure_data_audit.json").write_text(json.dumps(audit, indent=2)+"\n")


if __name__ == "__main__":
    main()
