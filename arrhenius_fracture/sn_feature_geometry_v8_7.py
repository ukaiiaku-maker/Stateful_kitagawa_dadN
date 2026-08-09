"""Generalized fixed edge-feature geometry for Stateful-PD v8.7.

This module preserves the v8.3 geometry API while adding two resolved feature
families and a refined prospective crack corridor:

* ``ellipse``: the existing half-elliptical edge notch;
* ``rounded_v``: a V-notch with independent depth, root radius, and included
  opening angle.

The global FEM remains intact and the geometry is fixed during production runs.
The added corridor supplies a representative local/front spacing for large
features whose PD patch also contains a coarse outer region.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial import Delaunay

from .mesh import (
    TriMesh,
    BoundaryData,
    _estimate_hbar,
    _estimate_hbar_tip,
    _precompute_element_data,
)


@dataclass
class BluntNotchGeometry:
    Lx: float = 2.0e-3
    Ly: float = 4.0e-3
    depth_a: float = 0.15e-3
    half_height_b: float = 0.30e-3
    feature_type: str = "ellipse"
    root_radius_m: float | None = None
    opening_angle_deg: float = 60.0
    path_refine_length_m: float = 0.24e-3
    path_refine_half_height_m: float = 60e-6

    def __post_init__(self) -> None:
        self.feature_type = str(self.feature_type).strip().lower()
        if self.feature_type not in {"ellipse", "rounded_v"}:
            raise ValueError(f"unsupported feature_type={self.feature_type!r}")
        if self.depth_a <= 0 or self.Lx <= self.depth_a or self.Ly <= 0:
            raise ValueError("invalid feature/domain dimensions")
        if self.feature_type == "ellipse" and self.half_height_b <= 0:
            raise ValueError("ellipse half-height must be positive")
        if self.feature_type == "rounded_v":
            rho = self.root_radius
            theta = math.radians(0.5 * self.opening_angle_deg)
            if not (0 < theta < 0.5 * math.pi):
                raise ValueError("rounded-V included angle must lie between 0 and 180 degrees")
            if not (0 < rho < self.depth_a):
                raise ValueError("rounded-V root radius must be positive and smaller than depth")
            if self.virtual_apex_x <= self.depth_a:
                raise ValueError("rounded-V virtual apex must lie beyond the rounded root")

    @property
    def root_xy(self) -> tuple[float, float]:
        return (float(self.depth_a), 0.0)

    @property
    def root_radius(self) -> float:
        if self.feature_type == "ellipse":
            return float(self.half_height_b**2 / max(self.depth_a, 1e-30))
        if self.root_radius_m is not None and self.root_radius_m > 0:
            return float(self.root_radius_m)
        return float(self.half_height_b**2 / max(self.depth_a, 1e-30))

    @property
    def half_angle_rad(self) -> float:
        return math.radians(0.5 * float(self.opening_angle_deg))

    @property
    def circle_center_x(self) -> float:
        return float(self.depth_a - self.root_radius)

    @property
    def virtual_apex_x(self) -> float:
        theta = self.half_angle_rad
        return float(self.circle_center_x + self.root_radius / max(math.sin(theta), 1e-30))

    @property
    def tangent_half_height(self) -> float:
        return float(self.root_radius * math.cos(self.half_angle_rad))

    @property
    def mouth_half_height(self) -> float:
        if self.feature_type == "ellipse":
            return float(self.half_height_b)
        return float(self.virtual_apex_x * math.tan(self.half_angle_rad))

    def profile_x(self, y: np.ndarray | float) -> np.ndarray:
        yy = np.asarray(y, dtype=float)
        ay = np.abs(yy)
        out = np.full_like(yy, -np.inf, dtype=float)
        if self.feature_type == "ellipse":
            inside = ay <= self.half_height_b * (1.0 + 1e-12)
            q = np.clip(1.0 - (ay[inside] / max(self.half_height_b, 1e-30)) ** 2, 0.0, 1.0)
            out[inside] = self.depth_a * np.sqrt(q)
            return out

        rho = self.root_radius
        yt = self.tangent_half_height
        ym = self.mouth_half_height
        root_arc = ay <= yt * (1.0 + 1e-12)
        out[root_arc] = self.circle_center_x + np.sqrt(
            np.maximum(rho * rho - ay[root_arc] ** 2, 0.0)
        )
        flank = (ay > yt) & (ay <= ym * (1.0 + 1e-12))
        out[flank] = self.virtual_apex_x - ay[flank] / max(math.tan(self.half_angle_rad), 1e-30)
        return out

    def boundary_points(self, spacing: float) -> np.ndarray:
        spacing = max(float(spacing), 1e-12)
        if self.feature_type == "ellipse":
            n = max(48, int(math.pi * max(self.depth_a, self.half_height_b) / spacing))
            th = np.linspace(-0.5 * math.pi, 0.5 * math.pi, n)
            return np.column_stack([
                self.depth_a * np.cos(th),
                self.half_height_b * np.sin(th),
            ])

        rho = self.root_radius
        theta = self.half_angle_rad
        phi_t = 0.5 * math.pi - theta
        n_arc = max(33, int(2.0 * phi_t * rho / spacing) + 1)
        if n_arc % 2 == 0:
            n_arc += 1
        phi = np.linspace(-phi_t, phi_t, n_arc)
        arc = np.column_stack([
            self.circle_center_x + rho * np.cos(phi),
            rho * np.sin(phi),
        ])
        y_t = self.tangent_half_height
        y_m = self.mouth_half_height
        flank_len = math.hypot(self.virtual_apex_x - (self.circle_center_x + rho * math.sin(theta)), y_m - y_t)
        n_flank = max(12, int(flank_len / spacing) + 1)
        y_lower = np.linspace(-y_m, -y_t, n_flank, endpoint=False)
        y_upper = np.linspace(y_t, y_m, n_flank)
        lower = np.column_stack([self.profile_x(y_lower), y_lower])
        upper = np.column_stack([self.profile_x(y_upper), y_upper])
        return np.vstack([lower, arc, upper])


def _inside_void(points: np.ndarray, geom: BluntNotchGeometry, margin: float = 0.0) -> np.ndarray:
    p = np.asarray(points, float)
    x = p[:, 0]
    y = p[:, 1]
    xb = geom.profile_x(y)
    if margin:
        xb = xb + float(margin)
    return (x >= -1e-15) & np.isfinite(xb) & (x < xb - 1e-12)


def _triangle_crosses_void(nodes: np.ndarray, elems: np.ndarray, geom: BluntNotchGeometry) -> np.ndarray:
    pts = nodes[elems]
    samples = [pts.mean(axis=1)]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        pi, pj = pts[:, i, :], pts[:, j, :]
        for t in (0.2, 0.4, 0.5, 0.6, 0.8):
            samples.append((1.0 - t) * pi + t * pj)
    bad = np.zeros(len(elems), dtype=bool)
    for sample in samples:
        bad |= _inside_void(sample, geom)
    return bad


def _path_refinement_points(geom: BluntNotchGeometry, h: float) -> np.ndarray:
    length = max(float(geom.path_refine_length_m), 0.0)
    half = max(float(geom.path_refine_half_height_m), 0.0)
    if length <= 0 or half <= 0:
        return np.empty((0, 2), float)
    x0 = geom.depth_a + 0.35 * h
    x1 = min(geom.Lx - 0.5 * h, geom.depth_a + length)
    if x1 <= x0:
        return np.empty((0, 2), float)
    xs = np.arange(x0, x1 + 0.25 * h, h)
    ys = np.arange(-half, half + 0.25 * h, h)
    rows = []
    for j, y in enumerate(ys):
        offset = 0.5 * h if j % 2 else 0.0
        xx = xs + offset
        xx = xx[xx <= x1]
        if xx.size:
            rows.append(np.column_stack([xx, np.full_like(xx, y)]))
    if not rows:
        return np.empty((0, 2), float)
    q = np.vstack(rows)
    inside = (
        (q[:, 0] >= 0.0) & (q[:, 0] <= geom.Lx) &
        (q[:, 1] >= -geom.Ly / 2.0) & (q[:, 1] <= geom.Ly / 2.0)
    )
    q = q[inside]
    return q[~_inside_void(q, geom)]



def _merge_points_min_spacing(groups: list[np.ndarray], min_spacing: float) -> np.ndarray:
    """Merge prioritized point groups while suppressing near-duplicates.

    Feature-surface and corridor points are supplied before root rings and the
    coarse background.  A small spatial hash removes ring/corridor/background
    intersections that would otherwise create an artificially small nearest-
    neighbor spacing and a mesh-dependent nonlocal horizon ratio.
    """
    h = max(float(min_spacing), 1e-15)
    cell = h
    accepted: list[np.ndarray] = []
    buckets: dict[tuple[int, int], list[int]] = {}
    h2 = h * h
    for group in groups:
        for point in np.asarray(group, float):
            key = (math.floor(point[0] / cell), math.floor(point[1] / cell))
            keep = True
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    for idx in buckets.get((key[0] + di, key[1] + dj), ()):
                        delta = accepted[idx] - point
                        if float(delta @ delta) < h2:
                            keep = False
                            break
                    if not keep:
                        break
                if not keep:
                    break
            if keep:
                idx = len(accepted)
                accepted.append(point.copy())
                buckets.setdefault(key, []).append(idx)
    return np.asarray(accepted, float)

def make_blunt_edge_notch_mesh(
    geom: BluntNotchGeometry,
    nx: int = 60,
    ny: int = 120,
    jitter: float = 0.10,
    root_h_fine: float = 20e-6,
    seed: int = 42,
):
    """Create a resolved triangular mesh for an ellipse or rounded V-notch."""
    rng = np.random.default_rng(seed)
    xv = np.linspace(0.0, geom.Lx, nx + 1)
    yv = np.linspace(-geom.Ly / 2.0, geom.Ly / 2.0, ny + 1)
    hx = xv[1] - xv[0]
    hy = yv[1] - yv[0]

    pts = []
    for j, y in enumerate(yv):
        for i, x in enumerate(xv):
            xx, yy = x, y
            if 0 < i < nx and 0 < j < ny:
                xx += jitter * hx * rng.uniform(-1.0, 1.0)
                yy += jitter * hy * rng.uniform(-1.0, 1.0)
            pts.append((xx, yy))
    pts = np.asarray(pts, float)
    pts = pts[~_inside_void(pts, geom)]

    arc_len_scale = max(root_h_fine, min(hx, hy) * 0.35)
    feature = geom.boundary_points(arc_len_scale)

    root = np.asarray(geom.root_xy)
    rings = []
    r = max(root_h_fine, 0.25 * min(hx, hy))
    # Keep the original root-centered ring extent; the separate corridor
    # resolves the prospective crack path without filling the entire enlarged
    # patch with circumferential rings.
    rmax = min(0.35e-3, 0.25 * geom.Ly)
    while r <= rmax:
        nr = max(12, int(2 * math.pi * r / max(root_h_fine, 1e-30)))
        th = np.linspace(0, 2 * math.pi, nr, endpoint=False)
        q = root[None, :] + r * np.column_stack([np.cos(th), np.sin(th)])
        inside_domain = (
            (q[:, 0] >= 0.0) & (q[:, 0] <= geom.Lx) &
            (q[:, 1] >= -geom.Ly / 2.0) & (q[:, 1] <= geom.Ly / 2.0)
        )
        q = q[inside_domain]
        q = q[~_inside_void(q, geom)]
        rings.append(q)
        r *= 1.45

    corridor = _path_refinement_points(geom, max(root_h_fine, 1e-12))
    ring_points = np.vstack(rings) if rings else np.empty((0, 2), float)
    # Priority is physical feature surface, crack-path corridor, root rings,
    # then coarse background.  The 0.62h exclusion removes accidental point
    # pairs much closer than the intended local resolution while retaining a
    # sufficiently irregular triangulation.
    nodes = _merge_points_min_spacing(
        [feature, corridor, ring_points, pts],
        0.62 * max(root_h_fine, 1e-12),
    )

    tri = Delaunay(nodes)
    elems = tri.simplices.copy()
    cent = nodes[elems].mean(axis=1)
    in_box = (
        (cent[:, 0] >= 0.0) & (cent[:, 0] <= geom.Lx) &
        (cent[:, 1] >= -geom.Ly / 2.0) & (cent[:, 1] <= geom.Ly / 2.0)
    )
    elems = elems[in_box]
    elems = elems[~_triangle_crosses_void(nodes, elems, geom)]

    used = np.unique(elems.ravel())
    remap = -np.ones(len(nodes), dtype=int)
    remap[used] = np.arange(len(used))
    nodes = nodes[used]
    elems = remap[elems]

    area_e, dNdx_e, B_e = _precompute_element_data(nodes, elems)
    hbar = _estimate_hbar(nodes, elems)
    hbar_tip = _estimate_hbar_tip(nodes, elems, geom.depth_a, 0.0)
    mesh = TriMesh(
        nodes=nodes,
        elems=elems,
        nn=len(nodes),
        ne=len(elems),
        ndof=2 * len(nodes),
        hbar=hbar,
        area_e=area_e,
        dNdx_e=dNdx_e,
        B_e=B_e,
        hbar_tip=hbar_tip,
    )

    x, y = nodes[:, 0], nodes[:, 1]
    tol_b = max(0.75 * hbar, 1.5 * root_h_fine)
    top = np.where(np.abs(y - geom.Ly / 2.0) < tol_b)[0]
    bot = np.where(np.abs(y + geom.Ly / 2.0) < tol_b)[0]
    left_bot = int(np.argmin((x - 0.0) ** 2 + (y + geom.Ly / 2.0) ** 2))
    right_bot = int(np.argmin((x - geom.Lx) ** 2 + (y + geom.Ly / 2.0) ** 2))
    bnd = BoundaryData(
        top_nodes=top,
        bot_nodes=bot,
        left_bot=left_bot,
        right_bot=right_bot,
        notch_nodes=np.array([], dtype=int),
    )
    return mesh, bnd, geom.root_xy


def identify_feature_surface_nodes(mesh: TriMesh, geom: BluntNotchGeometry, n_sample: int = 400) -> np.ndarray:
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    xb = geom.profile_x(y)
    scale = max(geom.depth_a, geom.mouth_half_height, 1e-12)
    tol = max(1e-12, 2e-7 * scale)
    mask = np.isfinite(xb) & (x >= -1e-12) & (np.abs(x - xb) <= tol)
    idx = np.where(mask)[0]
    if len(idx) < 5:
        tol = max(1e-10, 2e-4 * scale)
        idx = np.where(np.isfinite(xb) & (x >= -1e-12) & (np.abs(x - xb) <= tol))[0]
    if len(idx) < 5:
        samples = geom.boundary_points(scale / max(n_sample, 40))
        from scipy.spatial import cKDTree
        tree = cKDTree(mesh.nodes)
        _, nearest = tree.query(samples, k=1)
        idx = np.unique(nearest)
    return idx[np.argsort(mesh.nodes[idx, 1])]


def feature_tangent_normal(mesh: TriMesh, feature_nodes: np.ndarray):
    idx = np.asarray(feature_nodes, dtype=int)
    q = mesh.nodes[idx]
    npt = len(idx)
    t = np.zeros((npt, 2), dtype=float)
    if npt < 2:
        return t, t.copy()
    t[0] = q[1] - q[0]
    t[-1] = q[-1] - q[-2]
    if npt > 2:
        t[1:-1] = q[2:] - q[:-2]
    length = np.linalg.norm(t, axis=1)
    t /= np.maximum(length[:, None], 1e-30)
    n = np.column_stack([-t[:, 1], t[:, 0]])
    # The removed edge feature always lies locally toward decreasing x.
    flip = n[:, 0] > 0.0
    n[flip] *= -1.0
    return t, n


def rebuild_mesh_geometry(mesh: TriMesh, root_xy=None) -> TriMesh:
    area_e, dNdx_e, B_e = _precompute_element_data(mesh.nodes, mesh.elems)
    mesh.area_e = area_e
    mesh.dNdx_e = dNdx_e
    mesh.B_e = B_e
    mesh.hbar = _estimate_hbar(mesh.nodes, mesh.elems)
    if root_xy is None:
        root_xy = mesh.nodes[np.argmax(mesh.nodes[:, 0])]
    mesh.hbar_tip = _estimate_hbar_tip(mesh.nodes, mesh.elems, float(root_xy[0]), float(root_xy[1]))
    return mesh


def local_root_xy(mesh: TriMesh, feature_nodes: np.ndarray) -> tuple[float, float]:
    idx = np.asarray(feature_nodes, dtype=int)
    j = int(idx[np.argmax(mesh.nodes[idx, 0])])
    return float(mesh.nodes[j, 0]), float(mesh.nodes[j, 1])


def local_root_radius(mesh: TriMesh, feature_nodes: np.ndarray) -> float:
    idx = np.asarray(feature_nodes, dtype=int)
    q = mesh.nodes[idx]
    ir = int(np.argmax(q[:, 0]))
    if ir <= 0 or ir >= len(q) - 1:
        return float("nan")
    p1, p2, p3 = q[ir - 1], q[ir], q[ir + 1]
    a = np.linalg.norm(p2 - p1)
    b = np.linalg.norm(p3 - p2)
    c = np.linalg.norm(p3 - p1)
    v21 = p2 - p1
    v31 = p3 - p1
    area2 = abs(v21[0] * v31[1] - v21[1] * v31[0])
    if area2 <= 1e-30:
        return float("inf")
    return float(a * b * c / (2.0 * area2))


def apply_local_ale_surface_update(
    mesh: TriMesh,
    feature_nodes: np.ndarray,
    feature_normal_displacement: np.ndarray,
    decay_length: float,
    fixed_nodes: np.ndarray | None = None,
    max_move: float | None = None,
    min_area_fraction: float = 0.15,
):
    """Compatibility ALE helper; production v8.7 optimization keeps geometry fixed."""
    from scipy.spatial import cKDTree

    idx = np.asarray(feature_nodes, dtype=int)
    dh = np.asarray(feature_normal_displacement, dtype=float)
    if len(idx) == 0 or np.max(np.abs(dh)) <= 0.0:
        return 0.0
    _, normal = feature_tangent_normal(mesh, idx)
    disp_feature = normal * dh[:, None]
    if max_move is not None and np.isfinite(max_move) and max_move > 0:
        m = float(np.max(np.linalg.norm(disp_feature, axis=1)))
        if m > max_move:
            disp_feature *= max_move / max(m, 1e-30)

    tree = cKDTree(mesh.nodes[idx])
    dist, nearest = tree.query(mesh.nodes, k=1)
    decay = max(float(decay_length), 1e-30)
    weight = np.exp(-dist / decay)
    displacement = disp_feature[nearest] * weight[:, None]
    if fixed_nodes is not None:
        displacement[np.asarray(fixed_nodes, dtype=int)] = 0.0

    nodes0 = mesh.nodes.copy()
    min_area0 = float(np.min(mesh.area_e))
    scale = 1.0
    while scale >= 1e-5:
        mesh.nodes[:] = nodes0 + scale * displacement
        try:
            rebuild_mesh_geometry(mesh)
            if np.all(np.isfinite(mesh.area_e)) and np.min(mesh.area_e) >= min_area_fraction * min_area0:
                return float(scale)
        except Exception:
            pass
        scale *= 0.5
    mesh.nodes[:] = nodes0
    rebuild_mesh_geometry(mesh)
    return 0.0
