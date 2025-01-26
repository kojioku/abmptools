import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import argparse
import datetime
import sys

def read_gro_file(filename, verbose=False):
    """
    groファイルを読み取り、以下を返す:
      - atoms: [(resid, res_name, atom_name, atom_id, x, y, z), ...] のリスト
      - box: (lx, ly, lz)
    """
    atoms = []
    box = (0.0, 0.0, 0.0)
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 最終行にボックスサイズがあると仮定
    # 例: "  23.29635  23.29635  17.80809"
    box_line = lines[-1].split()
    if len(box_line) < 3:
        print("Error: ボックスサイズが正しく読み取れません。", file=sys.stderr)
        sys.exit(1)
    lx, ly, lz = map(float, box_line[:3])
    box = (lx, ly, lz)
    if verbose:
        print(f"Box dimensions: lx={lx}, ly={ly}, lz={lz}")

    # 粒子行は 3 行目から (タイトル行を 1 行目とした場合) 最終行(ボックス行)の手前まで
    for idx, line in enumerate(lines[2:-1], start=3):
        # groファイルの典型的な列配置を想定
        # 例:
        #   1M0000  A1    1   1.0   1.0   5.0
        try:
            resid = int(line[0:5].strip())  # 1M0000など
            res_name = line[5:10].strip()    # A1など
            atom_name = line[10:15].strip()  # A1など
            atom_id = int(line[15:20].strip())  # 1, 2, ...
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
            atoms.append((resid, res_name, atom_name, atom_id, x, y, z))
        except ValueError as e:
            print(f"Error parsing line {idx}: {line.strip()}", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)

    if verbose:
        print(f"Total atoms read: {len(atoms)}")
    return atoms, box


def make_2d_points_for_pbc(points, box, replicate=1, verbose=False):
    """
    周期境界を考慮するために、中心セルを含むように points を
    x,y 平面で -replicate ~ +replicate の範囲で複製する。
    points: Nx2 の ndarray
    box: (lx, ly)
    replicate: int
    return: (拡張された点集合, インデックス対応, 中央セルのインデックス)
        - extended_points: Mx2 の ndarray
        - mapping_index: M 要素のリスト(拡張された各点が「元の何番目か」を示す)
        - center_indices: 中央セルに対応する拡張点のインデックス
    """
    lx, ly = box[:2]
    extended_points = []
    mapping_index = []
    center_indices = []
    shifts = []
    for ix in range(-replicate, replicate+1):
        for iy in range(-replicate, replicate+1):
            shifts.append((ix * lx, iy * ly))

    if verbose:
        print(f"Shifts for replication (total {len(shifts)}):")
        for shift in shifts:
            print(f"  Shift: ({shift[0]}, {shift[1]})")

    for shift_idx, (sx, sy) in enumerate(shifts):
        for i, (x, y) in enumerate(points):
            extended_x = x + sx
            extended_y = y + sy
            extended_points.append([extended_x, extended_y])
            mapping_index.append(i)
            # 中央セルのシフトの場合、shift=(0,0)のみ
            if sx == 0 and sy == 0:
                center_indices.append(len(extended_points) - 1)
            if verbose and shift_idx == 0:
                print(f"  Adding central point {i}: ({extended_x}, {extended_y})")

    extended_points = np.array(extended_points)
    mapping_index = np.array(mapping_index)
    center_indices = np.array(center_indices)

    if verbose:
        print(f"Total extended points: {len(extended_points)}")
        print(f"Central cell indices (should be {len(points)}): {center_indices}")
    return extended_points, mapping_index, center_indices


def compute_voronoi_2d_pbc(points_2d, box, verbose=False):
    """
    2D の点 (points_2d) に対して、周期境界を考慮したボロノイ分割を行う。
    周囲に 1 タイルずつ複製して通常の Voronoi を計算し、
    中心セル領域のみ取り出す。
    return: Voronoi結果, 中央セルのインデックス, マッピングインデックス
    """
    lx, ly = box[:2]
    # 周囲に 1 タイル複製して点拡張
    extended_points, mapping_idx, center_indices = make_2d_points_for_pbc(points_2d, (lx, ly), replicate=1, verbose=verbose)
    if verbose:
        print("Computing Voronoi tessellation with PBC...")
    vor = Voronoi(extended_points)
    if verbose:
        print("Voronoi tessellation completed.")
    return vor, center_indices, mapping_idx


def compute_voronoi_2d(points_2d, verbose=False):
    """
    2D の点 (points_2d) に対して、ボロノイ分割を行う。
    return: Voronoi結果
    """
    if verbose:
        print("Computing Voronoi tessellation without PBC...")
    vor = Voronoi(points_2d)
    if verbose:
        print("Voronoi tessellation completed.")
    return vor


def compute_neighbor_counts(coords_2d, box, pbc, verbose=False, vmin=None, vmax=None):
    if len(coords_2d) < 2:
        if verbose:
            print("Not enough points to compute neighbors.")
        return coords_2d, np.zeros(len(coords_2d), dtype=int), None

    if pbc:
        vor, center_indices, mapping_idx = compute_voronoi_2d_pbc(coords_2d[:, :2], box, verbose=verbose)
        # center_indices は拡張点のインデックス
        # mapping_idx[center_indices] が元の点のインデックス
        neighbors_dict = {i: set() for i in range(len(coords_2d))}

        if verbose:
            print("Processing ridge points for PBC...")

        for ridge_idx, ridge_points in enumerate(vor.ridge_points):
            i, j = ridge_points
            # Check if either point is in the central cell
            if i in center_indices or j in center_indices:
                orig_i = mapping_idx[i]
                orig_j = mapping_idx[j]
                if orig_i != orig_j:
                    neighbors_dict[orig_i].add(orig_j)
                    neighbors_dict[orig_j].add(orig_i)
                    if verbose:
                        print(f"Ridge {ridge_idx}: Point {orig_i} <-> Point {orig_j}")

        if verbose:
            print("Completed processing ridge points.")

        neighbor_counts = np.array([len(neighbors_dict[i]) for i in range(len(coords_2d))])

        if verbose:
            print("Neighbor counts (PBC):")
            for i, count in enumerate(neighbor_counts):
                print(f"  Point {i}: {count} neighbors")

        return coords_2d, neighbor_counts, vor
    else:
        vor = compute_voronoi_2d(coords_2d[:, :2], verbose=verbose)
        neighbors_dict = {i: set() for i in range(len(coords_2d))}

        if verbose:
            print("Processing ridge points without PBC...")

        for ridge_idx, ridge_points in enumerate(vor.ridge_points):
            i, j = ridge_points
            if i < len(coords_2d) and j < len(coords_2d):
                neighbors_dict[i].add(j)
                neighbors_dict[j].add(i)
                if verbose:
                    print(f"Ridge {ridge_idx}: Point {i} <-> Point {j}")

        neighbor_counts = np.array([len(neighbors_dict[i]) for i in range(len(coords_2d))])

        if verbose:
            print("Neighbor counts (no PBC):")
            for i, count in enumerate(neighbor_counts):
                print(f"  Point {i}: {count} neighbors")

        return coords_2d, neighbor_counts, vor


def get_neighbor_count_same_type(atoms, box, target_atom_names=None, exclude_resnames=None, z_threshold=8.0, pbc=True, verbose=False):
    """
    与えられた全 atoms のうち、指定された原子を基準にして
    2D Voronoi を計算し、「ボロノイセルが隣接している同種原子の数」を数える。

    return:
      - coords_2d_upper: shape=(N_upper,3) の ndarray
      - neighbor_counts_upper: shape=(N_upper,) の ndarray
      - vor_upper: Voronoiオブジェクト
      - coords_2d_lower: shape=(N_lower,3) の ndarray
      - neighbor_counts_lower: shape=(N_lower,) の ndarray
      - vor_lower: Voronoiオブジェクト
      - molecule_centers: 全分子の重心情報
    """
    if target_atom_names:
        target_atoms = [a for a in atoms if a[2] in target_atom_names]
        if verbose:
            print(f"Filtering atoms by target names: {target_atom_names}")
    else:
        target_atoms = atoms
        if verbose:
            print("No target atom names specified; using all atoms.")

    if exclude_resnames:
        target_atoms = [a for a in target_atoms if a[1] not in exclude_resnames]
        if verbose:
            print(f"Excluding residue names: {exclude_resnames}")

    if verbose:
        print(f"Total target atoms after filtering: {len(target_atoms)}")

    # 分子ごとに重心を計算
    molecules = {}
    for atom in target_atoms:
        resid = atom[0]
        if resid not in molecules:
            molecules[resid] = []
        molecules[resid].append(atom)

    molecule_centers = []
    for resid, atoms_in_mol in molecules.items():
        x_mean = np.mean([a[4] for a in atoms_in_mol])
        y_mean = np.mean([a[5] for a in atoms_in_mol])
        z_mean = np.mean([a[6] for a in atoms_in_mol])
        molecule_centers.append((x_mean, y_mean, z_mean, resid))

    if verbose:
        print(f"Computed centers for {len(molecule_centers)} molecules.")

    # 中心座標を基準にして2D座標を取得
    coords_2d_upper = np.array([(x, y, resid) for x, y, z, resid in molecule_centers if z > z_threshold])
    coords_2d_lower = np.array([(x, y, resid) for x, y, z, resid in molecule_centers if z <= z_threshold])

    if verbose:
        print(f"Upper layer points (z > {z_threshold}): {len(coords_2d_upper)}")
        print(f"Lower layer points (z <= {z_threshold}): {len(coords_2d_lower)}")

    def compute_neighbor_counts_wrapper(coords_2d, box, pbc):
        return compute_neighbor_counts(coords_2d, box, pbc, verbose=verbose)

    coords_2d_upper, neighbor_counts_upper, vor_upper = compute_neighbor_counts_wrapper(coords_2d_upper, box, pbc)
    coords_2d_lower, neighbor_counts_lower, vor_lower = compute_neighbor_counts_wrapper(coords_2d_lower, box, pbc)

    return (coords_2d_upper, neighbor_counts_upper, vor_upper), (coords_2d_lower, neighbor_counts_lower, vor_lower), molecule_centers


def main():
    parser = argparse.ArgumentParser(description="Voronoi analysis of a gro file.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the gro file.")
    parser.add_argument("-a", "--atoms", type=str, nargs='+', help="Target atom names for Voronoi analysis.")
    parser.add_argument("-e", "--exclude", type=str, nargs='+', help="Residue names to exclude from analysis.")
    parser.add_argument("-z", "--z_threshold", type=float, default=8.0, help="Z-coordinate threshold to separate layers.")
    parser.add_argument("-nopbc", action='store_true', help="Do not consider periodic boundary conditions.")
    parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose debug output.")
    parser.add_argument("--vmin", type=float, default=None, help="Minimum value for color scaling.")
    parser.add_argument("--vmax", type=float, default=None, help="Maximum value for color scaling.")
    args = parser.parse_args()

    pbc = not args.nopbc
    verbose = args.verbose
    if verbose:
        print(f"Periodic Boundary Conditions: {'Enabled' if pbc else 'Disabled'}")

    atoms, box = read_gro_file(args.input, verbose=verbose)
    (coords_2d_upper, neighbor_counts_upper, vor_upper), (coords_2d_lower, neighbor_counts_lower, vor_lower), molecule_centers = get_neighbor_count_same_type(
        atoms, box, args.atoms, args.exclude, args.z_threshold, pbc, verbose=verbose)

    if len(coords_2d_upper) == 0 and len(coords_2d_lower) == 0:
        print("対象の原子は見つかりませんでした。")
        return

    # タイムスタンプを生成
    timestamp = datetime.datetime.now().strftime("%m%d_%H%M%S")

    # ボロノイ図を描画して保存
    plt.figure(figsize=(12, 6))

    if len(coords_2d_upper) > 0:
        plt.subplot(1, 2, 1)
        voronoi_plot_2d(vor_upper, ax=plt.gca(), show_vertices=False, line_colors='k', line_width=1, line_alpha=0.6, point_size=2)
        sc_upper = plt.scatter(
            coords_2d_upper[:, 0], coords_2d_upper[:, 1],
            c=neighbor_counts_upper, cmap='jet', s=100, edgecolors='k',
            vmin=args.vmin, vmax=args.vmax
        )
        plt.colorbar(sc_upper, label="Number of same-type neighbors")
        plt.title("Upper layer Voronoi neighbor count")
        plt.xlim(0, box[0])
        plt.ylim(0, box[1])
        plt.gca().set_aspect('equal', adjustable='box')

    if len(coords_2d_lower) > 0:
        plt.subplot(1, 2, 2)
        voronoi_plot_2d(vor_lower, ax=plt.gca(), show_vertices=False, line_colors='k', line_width=1, line_alpha=0.6, point_size=2)
        sc_lower = plt.scatter(
            coords_2d_lower[:, 0], coords_2d_lower[:, 1],
            c=neighbor_counts_lower, cmap='jet', s=100, edgecolors='k',
            vmin=args.vmin, vmax=args.vmax
        )
        plt.colorbar(sc_lower, label="Number of same-type neighbors")
        plt.title("Lower layer Voronoi neighbor count")
        plt.xlim(0, box[0])
        plt.ylim(0, box[1])
        plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout()

    # ファイル名を生成
    atom_names_str = "_".join(args.atoms) if args.atoms else "all"
    exclude_str = "_".join(args.exclude) if args.exclude else "none"
    pbc_str = "nopbc" if args.nopbc else "pbc"
    voronoi_filename = f"voronoi_diagram_{atom_names_str}_exclude_{exclude_str}_z{args.z_threshold}_{pbc_str}_{timestamp}.png"
    if verbose:
        print(f"Saving Voronoi diagram to {voronoi_filename}")
    plt.savefig(voronoi_filename, dpi=300)
    plt.close()

    # 隣接数のプロットを描画して保存
    plt.figure(figsize=(12, 6))

    if len(coords_2d_upper) > 0:
        plt.subplot(1, 2, 1)
        sc_upper = plt.scatter(
            coords_2d_upper[:, 0], coords_2d_upper[:, 1],
            c=neighbor_counts_upper, cmap='jet', s=100, edgecolors='k',
            vmin=args.vmin, vmax=args.vmax
        )
        plt.colorbar(sc_upper, label="Number of same-type neighbors")
        plt.title("Upper layer Voronoi neighbor count")
        plt.xlim(0, box[0])
        plt.ylim(0, box[1])
        plt.gca().set_aspect('equal', adjustable='box')

    if len(coords_2d_lower) > 0:
        plt.subplot(1, 2, 2)
        sc_lower = plt.scatter(
            coords_2d_lower[:, 0], coords_2d_lower[:, 1],
            c=neighbor_counts_lower, cmap='jet', s=100, edgecolors='k',
            vmin=args.vmin, vmax=args.vmax
        )
        plt.colorbar(sc_lower, label="Number of same-type neighbors")
        plt.title("Lower layer Voronoi neighbor count")
        plt.xlim(0, box[0])
        plt.ylim(0, box[1])
        plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout()

    neighbor_count_filename = f"neighbor_count_{atom_names_str}_exclude_{exclude_str}_z{args.z_threshold}_{pbc_str}_{timestamp}.png"
    if verbose:
        print(f"Saving neighbor count plot to {neighbor_count_filename}")
    plt.savefig(neighbor_count_filename, dpi=300)
    plt.close()

    # ログファイルに隣接数を出力
    if len(coords_2d_upper) > 0:
        log_filename_upper = f"neighbor_count_log_upper_{atom_names_str}_exclude_{exclude_str}_z{args.z_threshold}_{pbc_str}_{timestamp}.txt"
        if verbose:
            print(f"Writing neighbor counts to {log_filename_upper}")
        with open(log_filename_upper, 'w') as log_file:
            log_file.write(f"Periodic Boundary Conditions: {pbc}\n")
            log_file.write("Resid\tX\tY\tNeighborCount\n")
            for (x, y, resid), count in zip(coords_2d_upper, neighbor_counts_upper):
                log_file.write(f"{int(resid)}\t{x}\t{y}\t{count}\n")

    if len(coords_2d_lower) > 0:
        log_filename_lower = f"neighbor_count_log_lower_{atom_names_str}_exclude_{exclude_str}_z{args.z_threshold}_{pbc_str}_{timestamp}.txt"
        if verbose:
            print(f"Writing neighbor counts to {log_filename_lower}")
        with open(log_filename_lower, 'w') as log_file:
            log_file.write(f"Periodic Boundary Conditions: {pbc}\n")
            log_file.write("Resid\tX\tY\tNeighborCount\n")
            for (x, y, resid), count in zip(coords_2d_lower, neighbor_counts_lower):
                log_file.write(f"{int(resid)}\t{x}\t{y}\t{count}\n")


if __name__ == "__main__":
    main()

'''
python calcvolonoi-xy.py -i DPPC128-1746_441_53_record100.gro -e M0128 -z 8.0 -v
python calcvolonoi-xy.py -i DPPC128-1746_441_53_record100.gro -a Na1 Qa4 -z 8.0 -v
'''
