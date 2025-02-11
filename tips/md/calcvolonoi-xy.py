"""
Voronoi analysis of a gro file with neighbor count coloring.

- 分子内原子番号は1からとし、指定された原子の座標をそのまま対象とする。
- z_threshold のデフォルトはボックスの z 寸法の半分。
- pbcはデフォルトでON。
- Voronoi分割はpbcを反映して生成する。
- 隣接数が指定値（デフォルト6）以上のセルの面積の総和を出力する。
- さらに、通常のカラーマッピングとは別に、隣接数が n_threshold 以上の点を赤、
  それ以外を灰色で描画した画像も出力する。
- "distance" モードの場合は、隣接と判断された各原子間の残基IDと距離もログ出力する。
"""

import argparse
import datetime
import sys

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from scipy.spatial import Voronoi


def polygon_area(polygon):
    """Calculate area of a polygon using the shoelace formula."""
    area = 0.0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return abs(area) / 2.0


def read_gro_file(filename, verbose=False):
    """
    Read a gro file and return:
      - atoms: list of tuples (resid, res_name, atom_name, atom_id, x, y, z)
      - box: tuple (lx, ly, lz)
    """
    atoms = []
    box = (0.0, 0.0, 0.0)
    with open(filename, "r") as f:
        lines = f.readlines()

    box_line = lines[-1].split()
    if len(box_line) < 3:
        print("Error: ボックスサイズが正しく読み取れません.",
              file=sys.stderr)
        sys.exit(1)
    lx, ly, lz = map(float, box_line[:3])
    box = (lx, ly, lz)
    if verbose:
        print(f"Box dimensions: lx={lx}, ly={ly}, lz={lz}")

    for idx, line in enumerate(lines[2:-1], start=3):
        try:
            resid = int(line[0:5].strip())
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_id = int(line[15:20].strip())
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
            atoms.append((resid, res_name, atom_name, atom_id, x, y, z))
        except ValueError as e:
            print(f"Error parsing line {idx}: {line.strip()}",
                  file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)
    if verbose:
        print(f"Total atoms read: {len(atoms)}")
    return atoms, box


def make_2d_points_for_pbc(points, box, replicate=1, verbose=False):
    """
    Replicate 2D points for PBC in the x,y plane.
    Returns:
      - extended_points: Mx2 ndarray
      - mapping_index: array indicating which original point each replicated
        point corresponds to
      - center_indices: indices of the central (unshifted) points
    """
    lx, ly = box[:2]
    extended_points = []
    mapping_index = []
    center_indices = []
    shifts = []
    for ix in range(-replicate, replicate + 1):
        for iy in range(-replicate, replicate + 1):
            shifts.append((ix * lx, iy * ly))
    if verbose:
        print(f"Shifts for replication (total {len(shifts)}): {shifts}")
    for shift_idx, (sx, sy) in enumerate(shifts):
        for i, (x, y) in enumerate(points):
            extended_points.append([x + sx, y + sy])
            mapping_index.append(i)
            if sx == 0 and sy == 0:
                center_indices.append(len(extended_points) - 1)
    extended_points = np.array(extended_points)
    mapping_index = np.array(mapping_index)
    center_indices = np.array(center_indices)
    if verbose:
        print(f"Total extended points: {len(extended_points)}")
        print(f"Central cell indices (should be {len(points)}): {center_indices}")
    return extended_points, mapping_index, center_indices


def compute_voronoi_2d_pbc(points_2d, box, verbose=False):
    """
    Compute a 2D Voronoi tessellation with PBC.
    """
    lx, ly = box[:2]
    extended_points, mapping_idx, center_indices = make_2d_points_for_pbc(
        points_2d, (lx, ly), replicate=1, verbose=verbose)
    if verbose:
        print("Computing Voronoi tessellation with PBC...")
    vor = Voronoi(extended_points)
    if verbose:
        print("Voronoi tessellation completed.")
    return vor, center_indices, mapping_idx


def compute_voronoi_2d(points_2d, verbose=False):
    """
    Compute a 2D Voronoi tessellation without PBC.
    """
    if verbose:
        print("Computing Voronoi tessellation without PBC...")
    vor = Voronoi(points_2d)
    if verbose:
        print("Voronoi tessellation completed.")
    return vor


def compute_neighbor_counts_by_distance(coords, dist_thresh, pbc, box,
                                        verbose=False):
    """
    Compute neighbor counts based on a distance threshold.
    coords: (N,3) ndarray (each row is a coordinate)
    Returns neighbor_counts (excluding self) and the distance matrix.
    """
    if len(coords) == 0:
        return np.array([]), None
    delta = coords[:, None, :] - coords[None, :, :]  # shape (N, N, 3)
    if pbc:
        box_vec = np.array(box[:3])
        delta = delta - np.round(delta / box_vec) * box_vec
    dist_matrix = np.linalg.norm(delta, axis=2)
    neighbor_counts = (dist_matrix <= dist_thresh).sum(axis=1) - 1
    if verbose:
        print("Distance matrix:\n", dist_matrix)
        print("Neighbor counts (by distance):", neighbor_counts)
    return neighbor_counts, dist_matrix


def compute_neighbor_counts_by_voronoi(coords_2d, verbose=False):
    """
    Compute neighbor counts based on Voronoi tessellation.
    coords_2d: (N,2) ndarray (x,y only)
    """
    vor = compute_voronoi_2d(coords_2d, verbose=verbose)
    neighbors_dict = {i: set() for i in range(len(coords_2d))}
    for ridge in vor.ridge_points:
        i, j = ridge
        if i < len(coords_2d) and j < len(coords_2d):
            neighbors_dict[i].add(j)
            neighbors_dict[j].add(i)
    neighbor_counts = np.array([len(neighbors_dict[i])
                                for i in range(len(coords_2d))])
    if verbose:
        print("Neighbor counts (by voronoi):", neighbor_counts)
    return neighbor_counts, vor


def plot_voronoi_colored(vor, points, neighbor_counts, vmin, vmax, ax,
                         center_indices=None, verbose=False):
    """
    Plot Voronoi cells colored by neighbor counts.
    points: (N,2) ndarray (x,y only)
    If center_indices is provided (from pbc replication), iterate over the
    original points.
    The target points are overplotted as small darkened dots.
    """
    if center_indices is not None:
        indices = range(len(neighbor_counts))
    else:
        indices = range(len(vor.point_region))
    regions = []
    colors = []
    for j, _ in enumerate(indices):
        idx = center_indices[j] if center_indices is not None else j
        region_index = vor.point_region[idx]
        region = vor.regions[region_index]
        if -1 in region or len(region) == 0:
            continue
        polygon = [vor.vertices[i] for i in region]
        regions.append(polygon)
        colors.append(neighbor_counts[j])
    if verbose:
        print("Plotting colored Voronoi cells...")
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    coll = PolyCollection(regions, array=np.array(colors), cmap="jet",
                          norm=norm, edgecolors="k")
    ax.add_collection(coll)
    point_norm = mcolors.Normalize(
        vmin=vmin if vmin is not None else np.min(neighbor_counts),
        vmax=vmax if vmax is not None else np.max(neighbor_counts))
    cmap = plt.get_cmap("jet")
    point_colors = cmap(point_norm(neighbor_counts))
    dark_colors = [(c[0] * 0.5, c[1] * 0.5, c[2] * 0.5, c[3])
                   for c in point_colors]
    ax.scatter(points[:, 0], points[:, 1], c=dark_colors, s=20,
               edgecolors="none")
    ax.set_xlim(points[:, 0].min() - 1, points[:, 0].max() + 1)
    ax.set_ylim(points[:, 1].min() - 1, points[:, 1].max() + 1)
    ax.set_aspect("equal", adjustable="box")
    return coll


def plot_points_only(points, neighbor_counts, vmin, vmax, ax, verbose=False):
    """
    Plot only the target points with darkened colors.
    """
    norm = mcolors.Normalize(
        vmin=vmin if vmin is not None else np.min(neighbor_counts),
        vmax=vmax if vmax is not None else np.max(neighbor_counts))
    cmap = plt.get_cmap("jet")
    point_colors = cmap(norm(neighbor_counts))
    dark_colors = [(c[0] * 0.5, c[1] * 0.5, c[2] * 0.5, c[3])
                   for c in point_colors]
    sc = ax.scatter(points[:, 0], points[:, 1], c=dark_colors, s=20,
                    edgecolors="none")
    ax.set_xlim(points[:, 0].min() - 1, points[:, 0].max() + 1)
    ax.set_ylim(points[:, 1].min() - 1, points[:, 1].max() + 1)
    ax.set_aspect("equal", adjustable="box")
    if verbose:
        print("Plotting points only.")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sc, sm


def plot_voronoi_red_gray(vor, points, neighbor_counts, n_threshold, ax,
                          center_indices=None, verbose=False):
    """
    Plot Voronoi cells such that cells with neighbor count >= n_threshold
    are colored red and others gray.
    """
    if center_indices is not None:
        indices = range(len(neighbor_counts))
    else:
        indices = range(len(vor.point_region))
    regions = []
    colors = []
    for j, _ in enumerate(indices):
        idx = center_indices[j] if center_indices is not None else j
        region_index = vor.point_region[idx]
        region = vor.regions[region_index]
        if -1 in region or len(region) == 0:
            continue
        polygon = [vor.vertices[i] for i in region]
        regions.append(polygon)
        if neighbor_counts[j] >= n_threshold:
            colors.append("red")
        else:
            colors.append("gray")
    if verbose:
        print("Plotting red/gray Voronoi cells...")
    coll = PolyCollection(regions, facecolors=colors, edgecolors="k")
    ax.add_collection(coll)
    ax.scatter(points[:, 0], points[:, 1],
               c=["red" if nc >= n_threshold else "gray"
                  for nc in neighbor_counts],
               s=20, edgecolors="none")
    ax.set_xlim(points[:, 0].min() - 1, points[:, 0].max() + 1)
    ax.set_ylim(points[:, 1].min() - 1, points[:, 1].max() + 1)
    ax.set_aspect("equal", adjustable="box")
    return coll


def plot_points_red_gray(points, neighbor_counts, n_threshold, ax, verbose=False):
    """
    Plot only the target points in red if neighbor count >= n_threshold,
    otherwise gray.
    """
    point_colors = ["red" if nc >= n_threshold else "gray"
                    for nc in neighbor_counts]
    sc = ax.scatter(points[:, 0], points[:, 1], c=point_colors, s=20,
                    edgecolors="none")
    ax.set_xlim(points[:, 0].min() - 1, points[:, 0].max() + 1)
    ax.set_ylim(points[:, 1].min() - 1, points[:, 1].max() + 1)
    ax.set_aspect("equal", adjustable="box")
    if verbose:
        print("Plotting red/gray points only.")
    return sc


def get_selected_atoms(atoms, target_atom_names=None, exclude_resnames=None,
                       molecule_filter=None, verbose=False):
    """
    If molecule_filter is specified, select atoms based on molecule name and
    in-group atom indices (1-indexed). Returns a list of selected atom tuples.
    """
    if molecule_filter is not None:
        mol_name = molecule_filter[0]
        try:
            target_intragroup_indices = set(map(int, molecule_filter[1:]))
        except ValueError:
            print("Error: 分子内の原子番号は整数で指定してください.",
                  file=sys.stderr)
            sys.exit(1)
        if verbose:
            print(f"Filtering atoms by molecule: Name = {mol_name}, "
                  f"In-group indices = {target_intragroup_indices}")
        groups = {}
        for a in atoms:
            resid, res_name, atom_name, atom_id, x, y, z = a
            if res_name == mol_name:
                groups.setdefault(resid, []).append(a)
        selected_atoms = []
        for resid, atom_list in groups.items():
            for idx, a in enumerate(atom_list, start=1):
                if idx in target_intragroup_indices:
                    selected_atoms.append(a)
        if verbose:
            print(f"Total selected atoms: {len(selected_atoms)}")
        return selected_atoms
    else:
        if target_atom_names:
            selected_atoms = [a for a in atoms if a[2] in target_atom_names]
            if verbose:
                print(f"Filtering atoms by target names: {target_atom_names}")
        else:
            selected_atoms = atoms
            if verbose:
                print("No target atom names specified; using all atoms.")
        if exclude_resnames:
            selected_atoms = [a for a in selected_atoms if a[1] not in exclude_resnames]
            if verbose:
                print(f"Excluding residue names: {exclude_resnames}")
        return selected_atoms


def log_adjacency_details(filename, coords, dist_matrix, thresh, verbose=False):
    """
    For each selected atom, log its residue ID and the residue IDs of its
    neighbors (distance <= thresh) along with the distance.
    coords: list of tuples (x, y, z, resid)
    dist_matrix: distance matrix corresponding to coords
    """
    with open(filename, "w") as f:
        for i, (x, y, z, resid) in enumerate(coords):
            f.write(f"Residue {resid}:\n")
            for j, (nx, ny, nz, nresid) in enumerate(coords):
                if i != j and dist_matrix[i, j] <= thresh:
                    f.write(f"  Neighbor Residue {nresid} at distance "
                            f"{dist_matrix[i, j]:.3f} nm\n")
    if verbose:
        print(f"Adjacency details written to {filename}")


def main():
    parser = argparse.ArgumentParser(
        description=("Voronoi analysis of a gro file with neighbor count "
                     "coloring."))
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to the gro file.")
    parser.add_argument("-a", "--atoms", type=str, nargs='+',
                        help="Target atom names for analysis.")
    parser.add_argument("-e", "--exclude", type=str, nargs='+',
                        help="Residue names to exclude from analysis.")
    parser.add_argument("-z", "--z_threshold", type=float, default=None,
                        help=("Z-coordinate threshold to separate layers "
                              "(default: half of box z dimension)."))
    parser.add_argument("-nopbc", action="store_true",
                        help="Do not consider periodic boundary conditions.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose debug output.")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Minimum value for color scaling.")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Maximum value for color scaling.")
    parser.add_argument("-m", "--molecule", type=str, nargs='+',
                        help=("Specify molecule name and in-group atom indices "
                              "(e.g., -m DPPC 4 5)."))
    parser.add_argument("--dist_thresh", type=float, default=1.0,
                        help=("Distance threshold (nm) for neighbor determination "
                              "(used in 'distance' mode)."))
    parser.add_argument("--neighbor_mode", type=str,
                        choices=["distance", "voronoi"], default="distance",
                        help=("Method to determine neighbor count: 'distance' "
                              "(default) or 'voronoi'."))
    parser.add_argument("--n_threshold", type=int, default=6,
                        help=("Neighbor count threshold. The total area of "
                              "Voronoi cells with neighbor count >= this value "
                              "will be output."))
    args = parser.parse_args()

    pbc = not args.nopbc  # pbc is ON by default
    verbose = args.verbose
    if verbose:
        print(f"Periodic Boundary Conditions: {'Enabled' if pbc else 'Disabled'}")
    atoms, box = read_gro_file(args.input, verbose=verbose)
    if args.z_threshold is None:
        args.z_threshold = box[2] / 2.0
        if verbose:
            print(f"z_threshold not specified; using half of box z dimension: {args.z_threshold}")
    if args.molecule is None:
        print("Error: 分子内原子番号での解析(-mオプション)が必須です.",
              file=sys.stderr)
        sys.exit(1)
    selected_atoms = get_selected_atoms(
        atoms, target_atom_names=args.atoms, exclude_resnames=args.exclude,
        molecule_filter=args.molecule, verbose=verbose
    )
    # Separate into upper (z > threshold) and lower (z <= threshold)
    upper_atoms = [a for a in selected_atoms if a[6] > args.z_threshold]
    lower_atoms = [a for a in selected_atoms if a[6] <= args.z_threshold]
    if verbose:
        print(f"Number of selected atoms: {len(selected_atoms)}")
        print(f"Upper group atoms (z > {args.z_threshold}): {len(upper_atoms)}")
        print(f"Lower group atoms (z <= {args.z_threshold}): {len(lower_atoms)}")

    def atoms_to_coords(atom_list):
        coords = []
        for a in atom_list:
            resid = a[0]
            x, y, z = a[4], a[5], a[6]
            coords.append((x, y, z, resid))
        return coords

    upper_coords = atoms_to_coords(upper_atoms)
    lower_coords = atoms_to_coords(lower_atoms)

    upper_array = np.array([[x, y, z] for (x, y, z, _) in upper_coords])
    lower_array = np.array([[x, y, z] for (x, y, z, _) in lower_coords])

    if args.neighbor_mode == "distance":
        upper_neighbor_counts, upper_dist_matrix = (
            compute_neighbor_counts_by_distance(upper_array, args.dist_thresh,
                                                pbc, box, verbose=verbose)
        )
        lower_neighbor_counts, lower_dist_matrix = (
            compute_neighbor_counts_by_distance(lower_array, args.dist_thresh,
                                                pbc, box, verbose=verbose)
        )
    else:
        upper_neighbor_counts, _ = compute_neighbor_counts_by_voronoi(
            upper_array[:, :2], verbose=verbose)
        lower_neighbor_counts, _ = compute_neighbor_counts_by_voronoi(
            lower_array[:, :2], verbose=verbose)
    if verbose:
        print("Upper neighbor counts:", upper_neighbor_counts)
        print("Lower neighbor counts:", lower_neighbor_counts)

    if pbc:
        vor_upper, center_upper, _ = compute_voronoi_2d_pbc(
            upper_array[:, :2], box, verbose=verbose)
        vor_lower, center_lower, _ = compute_voronoi_2d_pbc(
            lower_array[:, :2], box, verbose=verbose)
    else:
        vor_upper = compute_voronoi_2d(upper_array[:, :2], verbose=verbose)
        vor_lower = compute_voronoi_2d(lower_array[:, :2], verbose=verbose)
        center_upper = None
        center_lower = None

    # Voronoi plotting for upper layer (jet colormap)
    if len(upper_array) > 0:
        fig_u, ax_u = plt.subplots(figsize=(8, 6))
        coll_u = plot_voronoi_colored(
            vor_upper, upper_array[:, :2], upper_neighbor_counts,
            args.vmin, args.vmax, ax_u,
            center_indices=(center_upper if pbc else None), verbose=verbose)
        ax_u.set_title(
            "Upper Layer Voronoi Colored by Neighbor Count\n(Method: {})"
            .format(args.neighbor_mode))
        ax_u.set_xlim(0, box[0])
        ax_u.set_ylim(0, box[1])
        cbar = plt.colorbar(coll_u, ax=ax_u)
        cbar.set_label("Neighbor Count")
        plt.tight_layout()
        timestamp = datetime.datetime.now().strftime("%m%d_%H%M%S")
        filename_upper = (
            f"voronoi_upper_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving upper layer Voronoi diagram to {filename_upper}")
        plt.savefig(filename_upper, dpi=300)
        plt.close()

        fig_u2, ax_u2 = plt.subplots(figsize=(8, 6))
        sc_u, sm_u = plot_points_only(
            upper_array[:, :2], upper_neighbor_counts, args.vmin,
            args.vmax, ax_u2, verbose=verbose)
        ax_u2.set_title("Upper Layer Points Colored by Neighbor Count")
        ax_u2.set_xlim(0, box[0])
        ax_u2.set_ylim(0, box[1])
        cbar_u = plt.colorbar(sm_u, ax=ax_u2)
        cbar_u.set_label("Neighbor Count")
        plt.tight_layout()
        filename_upper_points = (
            f"points_upper_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving upper layer points-only diagram to {filename_upper_points}")
        plt.savefig(filename_upper_points, dpi=300)
        plt.close()
    else:
        if verbose:
            print("No upper layer atoms found.")

    # Voronoi plotting for lower layer (jet colormap)
    if len(lower_array) > 0:
        fig_l, ax_l = plt.subplots(figsize=(8, 6))
        coll_l = plot_voronoi_colored(
            vor_lower, lower_array[:, :2], lower_neighbor_counts,
            args.vmin, args.vmax, ax_l,
            center_indices=(center_lower if pbc else None), verbose=verbose)
        ax_l.set_title(
            "Lower Layer Voronoi Colored by Neighbor Count\n(Method: {})"
            .format(args.neighbor_mode))
        ax_l.set_xlim(0, box[0])
        ax_l.set_ylim(0, box[1])
        cbar_l = plt.colorbar(coll_l, ax=ax_l)
        cbar_l.set_label("Neighbor Count")
        plt.tight_layout()
        filename_lower = (
            f"voronoi_lower_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving lower layer Voronoi diagram to {filename_lower}")
        plt.savefig(filename_lower, dpi=300)
        plt.close()

        fig_l2, ax_l2 = plt.subplots(figsize=(8, 6))
        sc_l, sm_l = plot_points_only(
            lower_array[:, :2], lower_neighbor_counts, args.vmin,
            args.vmax, ax_l2, verbose=verbose)
        ax_l2.set_title("Lower Layer Points Colored by Neighbor Count")
        ax_l2.set_xlim(0, box[0])
        ax_l2.set_ylim(0, box[1])
        cbar_l2 = plt.colorbar(sm_l, ax=ax_l2)
        cbar_l2.set_label("Neighbor Count")
        plt.tight_layout()
        filename_lower_points = (
            f"points_lower_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving lower layer points-only diagram to {filename_lower_points}")
        plt.savefig(filename_lower_points, dpi=300)
        plt.close()
    else:
        if verbose:
            print("No lower layer atoms found.")

    # Additional output: red/gray images.
    # For red/gray: neighbor count >= n_threshold -> red, else gray.
    if len(upper_array) > 0:
        fig_ur, ax_ur = plt.subplots(figsize=(8, 6))
        coll_ur = plot_voronoi_red_gray(
            vor_upper, upper_array[:, :2], upper_neighbor_counts,
            args.n_threshold, ax_ur,
            center_indices=(center_upper if pbc else None), verbose=verbose)
        ax_ur.set_title("Upper Layer Voronoi (Red/Gray)")
        ax_ur.set_xlim(0, box[0])
        ax_ur.set_ylim(0, box[1])
        plt.tight_layout()
        filename_upper_rg = (
            f"voronoi_upper_redgray_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving upper layer Voronoi red/gray diagram to {filename_upper_rg}")
        plt.savefig(filename_upper_rg, dpi=300)
        plt.close()

        fig_ur2, ax_ur2 = plt.subplots(figsize=(8, 6))
        sc_ur = plot_points_red_gray(
            upper_array[:, :2], upper_neighbor_counts,
            args.n_threshold, ax_ur2, verbose=verbose)
        ax_ur2.set_title("Upper Layer Points (Red/Gray)")
        ax_ur2.set_xlim(0, box[0])
        ax_ur2.set_ylim(0, box[1])
        plt.tight_layout()
        filename_upper_points_rg = (
            f"points_upper_redgray_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving upper layer points red/gray diagram to {filename_upper_points_rg}")
        plt.savefig(filename_upper_points_rg, dpi=300)
        plt.close()
    if len(lower_array) > 0:
        fig_lr, ax_lr = plt.subplots(figsize=(8, 6))
        coll_lr = plot_voronoi_red_gray(
            vor_lower, lower_array[:, :2], lower_neighbor_counts,
            args.n_threshold, ax_lr,
            center_indices=(center_lower if pbc else None), verbose=verbose)
        ax_lr.set_title("Lower Layer Voronoi (Red/Gray)")
        ax_lr.set_xlim(0, box[0])
        ax_lr.set_ylim(0, box[1])
        plt.tight_layout()
        filename_lower_rg = (
            f"voronoi_lower_redgray_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving lower layer Voronoi red/gray diagram to {filename_lower_rg}")
        plt.savefig(filename_lower_rg, dpi=300)
        plt.close()

        fig_lr2, ax_lr2 = plt.subplots(figsize=(8, 6))
        sc_lr = plot_points_red_gray(
            lower_array[:, :2], lower_neighbor_counts,
            args.n_threshold, ax_lr2, verbose=verbose)
        ax_lr2.set_title("Lower Layer Points (Red/Gray)")
        ax_lr2.set_xlim(0, box[0])
        ax_lr2.set_ylim(0, box[1])
        plt.tight_layout()
        filename_lower_points_rg = (
            f"points_lower_redgray_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{args.neighbor_mode}_{timestamp}.png"
        )
        if verbose:
            print(f"Saving lower layer points red/gray diagram to {filename_lower_points_rg}")
        plt.savefig(filename_lower_points_rg, dpi=300)
        plt.close()

    # Compute total area of Voronoi cells with neighbor count >= n_threshold.
    def compute_total_area_exceeding(vor, neighbor_counts, n_threshold,
                                     indices=None, verbose=False):
        total_area = 0.0
        if indices is None:
            for i in range(len(neighbor_counts)):
                region_index = vor.point_region[i]
                region = vor.regions[region_index]
                if -1 in region or len(region) == 0:
                    continue
                if neighbor_counts[i] >= n_threshold:
                    poly = [vor.vertices[j] for j in region]
                    total_area += polygon_area(poly)
        else:
            for j, ext_index in enumerate(indices):
                region_index = vor.point_region[ext_index]
                region = vor.regions[region_index]
                if -1 in region or len(region) == 0:
                    continue
                if neighbor_counts[j] >= n_threshold:
                    poly = [vor.vertices[k] for k in region]
                    total_area += polygon_area(poly)
        return total_area

    n_thresh = args.n_threshold
    if len(upper_array) > 0:
        total_area_upper = compute_total_area_exceeding(
            vor_upper, upper_neighbor_counts, n_thresh,
            indices=(center_upper if pbc else None), verbose=verbose)
        if verbose:
            print(f"Total area of upper layer cells with neighbor count >= "
                  f"{n_thresh}: {total_area_upper}")
    else:
        total_area_upper = 0.0
    if len(lower_array) > 0:
        total_area_lower = compute_total_area_exceeding(
            vor_lower, lower_neighbor_counts, n_thresh,
            indices=(center_lower if pbc else None), verbose=verbose)
        if verbose:
            print(f"Total area of lower layer cells with neighbor count >= "
                  f"{n_thresh}: {total_area_lower}")
    else:
        total_area_lower = 0.0
    total_area_all = total_area_upper + total_area_lower
    print(f"Total area of Voronoi cells with neighbor count >= {n_thresh}: "
          f"{total_area_all}")

    # If distance mode, log adjacency details.
    if args.neighbor_mode == "distance":
        upper_adj_log = (
            f"adjacency_details_upper_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{datetime.datetime.now().strftime('%m%d_%H%M%S')}.txt"
        )
        log_adjacency_details(upper_adj_log, upper_coords, upper_dist_matrix,
                              args.dist_thresh, verbose=verbose)
        lower_adj_log = (
            f"adjacency_details_lower_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
            f"{datetime.datetime.now().strftime('%m%d_%H%M%S')}.txt"
        )
        log_adjacency_details(lower_adj_log, lower_coords, lower_dist_matrix,
                              args.dist_thresh, verbose=verbose)

    log_filename_upper = (
        f"neighbor_count_log_upper_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
        f"{args.neighbor_mode}_{timestamp}.txt"
    )
    with open(log_filename_upper, "w") as f:
        f.write("Upper Layer:\n")
        f.write("Resid\tX\tY\tZ\tNeighborCount\n")
        for (x, y, z, resid), count in zip(upper_coords, upper_neighbor_counts):
            f.write(f"{resid}\t{x}\t{y}\t{z}\t{count}\n")
        f.write(f"Total area of cells with neighbor count >= {n_thresh}: "
                f"{total_area_upper}\n")
    log_filename_lower = (
        f"neighbor_count_log_lower_{args.molecule[0]}_{'_'.join(args.molecule[1:])}_"
        f"{args.neighbor_mode}_{timestamp}.txt"
    )
    with open(log_filename_lower, "w") as f:
        f.write("Lower Layer:\n")
        f.write("Resid\tX\tY\tZ\tNeighborCount\n")
        for (x, y, z, resid), count in zip(lower_coords, lower_neighbor_counts):
            f.write(f"{resid}\t{x}\t{y}\t{z}\t{count}\n")
        f.write(f"Total area of cells with neighbor count >= {n_thresh}: "
                f"{total_area_lower}\n")
    if verbose:
        print("Log files written.")


if __name__ == "__main__":
    main()

'''
python calcvolonoi-xy.py -i DPPC128-1746_441_53_record100.gro -m M0000 4 8 -v --dist_thresh 0.75
'''
