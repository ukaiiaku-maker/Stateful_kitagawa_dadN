#!/usr/bin/env python3
"""Read-only rendering of six representative Peak PD terminal states."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


ROOT = Path("runs/sn_v9_shared_root_m1_peak")
SPECS = (
    (2500, "B001", ROOT / "sn_matrix_v1/front_2500/B001", "captured"),
    (2500, "B000", ROOT / "sn_matrix_v1/front_2500/B000", "censored"),
    (2250, "B001", ROOT / "sn_matrix_v1/front_2250/B001", "captured"),
    (2250, "B000", ROOT / "sn_matrix_v1/front_2250/B000", "censored"),
    (1500, "B004", ROOT / "conditioned_ensemble_1500_v1/B004", "captured"),
    (1500, "B005", ROOT / "conditioned_ensemble_1500_v1/B005", "censored"),
)


def resolve_case(branch: Path, stress: int) -> Path:
    matches = list(branch.rglob(f"sigmaA_{stress}MPa/summary.json"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one saved {stress} MPa case under {branch}, found {len(matches)}")
    return matches[0].parent


def load_state(case: Path):
    required = (case / "summary.json", case / "pd_state_final.npz")
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    summary = json.loads(required[0].read_text())
    state = np.load(required[1], allow_pickle=False)
    return summary, state


def panel(ax, summary, z, title):
    xy = np.asarray(z["xy"], float); bonds = np.asarray(z["bonds"], int)
    damage = np.asarray(z["bond_damage"], float)
    seg = xy[bonds]
    light = damage < .95
    lc = LineCollection(seg[light] * 1e6, array=damage[light], cmap="viridis",
                        norm=plt.Normalize(0, 1), linewidths=.45, alpha=.72)
    ax.add_collection(lc)
    broken = damage >= .95
    if np.any(broken):
        ax.add_collection(LineCollection(seg[broken] * 1e6, colors="black", linewidths=1.7))
    boundary = np.asarray(z["feature_surface_xy"], float)
    if boundary.ndim == 2 and len(boundary):
        ax.plot(boundary[:, 0] * 1e6, boundary[:, 1] * 1e6, color="#555555", lw=1.2)
    attempt = summary.get("last_marked_attempt") or {}
    node = attempt.get("selected_node", summary.get("stable_crack_birth_node"))
    if node is not None and 0 <= int(node) < len(xy):
        ax.scatter(*(xy[int(node)] * 1e6), marker="*", s=95, c="#d95f02", edgecolor="white", lw=.5, zorder=8)
    seed = int(np.asarray(z["primary_seed_node"]).item())
    if 0 <= seed < len(xy):
        ax.scatter(*(xy[seed] * 1e6), s=70, facecolors="none", edgecolors="#e7298a", lw=1.7, zorder=9)
    if "front_backbone_bonds" in z.files:
        connected = np.asarray(z["front_backbone_bonds"], bool)
        if np.any(connected):
            ax.add_collection(LineCollection(seg[connected] * 1e6, colors="#00bfc4", linewidths=2.3))
    path = np.asarray(z["active_front_path_xy"], float)
    if path.ndim == 2 and len(path) > 1:
        ax.plot(path[:, 0] * 1e6, path[:, 1] * 1e6, "-o", color="#7b3294", lw=1.7, ms=3)
    ax.scatter(150, 0, marker="x", c="red", s=40, zorder=10)
    ax.set_xlim(40, 650); ax.set_ylim(-330, 330); ax.set_aspect("equal")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_title(title, fontsize=10, weight="bold")
    captured = summary.get("cycles_front_capture") is not None
    stable = summary.get("cycles_first_stable")
    exposure = float(summary["cycles_total"]) - float(stable) if stable is not None else np.nan
    prob = attempt.get("selected_probability", attempt.get("mark_probability", np.nan))
    note = (f"N={summary['cycles_total']:.4g}; ΔNstable={exposure:.3g}\n"
            f"node={node}; pmark={prob:.3g}; Dmax={summary['pd_bond_damage_final_max']:.3f}; "
            f"broken={summary['pd_broken_bonds_final']}")
    if captured: note += f"\nNfront={summary['cycles_front_capture']:.4g}"
    ax.text(.02, .02, note, transform=ax.transAxes, fontsize=7.2, va="bottom",
            bbox=dict(facecolor="white", alpha=.84, edgecolor="none", pad=2.5))
    ax.grid(False)
    return lc


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv); out = args.out; end = out / "end_states"; data = out / "figure_data"
    end.mkdir(parents=True, exist_ok=True); data.mkdir(parents=True, exist_ok=True)
    (out / "current").mkdir(parents=True, exist_ok=True)
    loaded = []
    for stress, branch_id, branch, status in SPECS:
        try:
            case = resolve_case(branch, stress); summary, z = load_state(case)
        except (FileNotFoundError, RuntimeError) as exc:
            print(f"SKIP optional branch {branch}: {exc}"); continue
        title = f"{stress} MPa — {'captured' if status == 'captured' else ('long censor' if stress == 1500 else '1e6 censor')}"
        fig, ax = plt.subplots(figsize=(8, 6), dpi=200); panel(ax, summary, z, title)
        target = end / f"{stress}_{branch_id}_{status}.png"
        fig.savefig(target, dpi=200, bbox_inches="tight", pil_kwargs={"compress_level": 9}); plt.close(fig)
        z.close(); loaded.append((stress, branch_id, case, status, target))
    if len(loaded) == 6:
        fig, axes = plt.subplots(3, 2, figsize=(10.5, 12), constrained_layout=True)
        labels = "abcdef"; mappable = None
        for i, ((stress, bid, case, status, _), ax) in enumerate(zip(loaded, axes.flat)):
            summary, z = load_state(case)
            qualifier = "captured" if status == "captured" else ("long censor" if stress == 1500 else "1e6 censor")
            mappable = panel(ax, summary, z, f"({labels[i]}) {stress} MPa — {qualifier}"); z.close()
        cb = fig.colorbar(mappable, ax=axes, shrink=.72, pad=.02); cb.set_label("bond damage")
        fig.suptitle("Representative Peak PD End States", fontsize=15, weight="bold")
        fig.savefig(out / "current/07_representative_end_states_montage.png", dpi=190,
                    bbox_inches="tight", pil_kwargs={"compress_level": 9})
        fig.savefig(out / "current/07_representative_end_states_montage.pdf", bbox_inches="tight")
        plt.close(fig)
    with (data / "end_state_sources.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(("sigma_a_MPa","branch_id","observation","source_path","image_path"))
        for stress,bid,case,status,target in loaded: w.writerow((stress,bid,status,case,target))


if __name__ == "__main__": main()
