import argparse
import os
import shutil
from Bio import PDB

def split_residue(pdb_file, residue_number, output_prefix, chain_id=None):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    io = PDB.PDBIO()

    # 残基を除いた構造（receptor）
    class ResidueExcluder(PDB.Select):
        def accept_residue(self, residue):
            if chain_id:
                return not (residue.get_id()[1] == residue_number and residue.get_parent().id == chain_id)
            else:
                return residue.get_id()[1] != residue_number

    io.set_structure(structure)
    receptor_file = f"{output_prefix}_receptor_{residue_number}.pdb"
    io.save(receptor_file, ResidueExcluder())
    print(f"Receptor file saved: {receptor_file}")

    # 残基のみの構造（ligand）
    class ResidueSelector(PDB.Select):
        def accept_residue(self, residue):
            if chain_id:
                return residue.get_id()[1] == residue_number and residue.get_parent().id == chain_id
            else:
                return residue.get_id()[1] == residue_number

    io.set_structure(structure)
    ligand_file = f"{output_prefix}_ligand_{residue_number}.pdb"
    io.save(ligand_file, ResidueSelector())
    print(f"Ligand file saved: {ligand_file}")

    return receptor_file, ligand_file

def process_pdb_file(pdb_file, residue_number, chain_id=None):
    pdb_basename = os.path.splitext(os.path.basename(pdb_file))[0]
    output_folder = pdb_basename

    # Create a directory for the PDB file
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output directory created: {output_folder}")

    # Copy the original PDB file
    shutil.copy(pdb_file, output_folder)
    print(f"Original PDB file copied to: {output_folder}")

    # Generate the output files
    receptor_file, ligand_file = split_residue(pdb_file, residue_number, pdb_basename, chain_id)

    # Move the generated PDB files
    shutil.move(receptor_file, output_folder)
    shutil.move(ligand_file, output_folder)
    print(f"Generated PDB files moved to: {output_folder}")

    return output_folder

def process_folder(input_folder, residue_number, chain_id=None):
    pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
    output_folders = []

    for pdb_file in pdb_files:
        pdb_path = os.path.join(input_folder, pdb_file)
        output_folder = process_pdb_file(pdb_path, residue_number, chain_id)
        output_folders.append(output_folder)

    return output_folders

def main():
    parser = argparse.ArgumentParser(description='Process PDB files to exclude and isolate a specific residue.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input PDB file or folder containing PDB files')
    parser.add_argument('-r', '--residue', type=int, required=True, help='Residue number to exclude and isolate')
    parser.add_argument('-c', '--chain', type=str, help='Chain ID of the residue to exclude and isolate')

    args = parser.parse_args()

    if os.path.isdir(args.input):
        output_folders = process_folder(args.input, args.residue, args.chain)
    else:
        output_folders = [process_pdb_file(args.input, args.residue, args.chain)]

    # Write the output folder names to dirnames.in
    with open('dirnames.in', 'w') as f:
        for folder in output_folders:
            f.write(f"{folder}\n")
    print("Output folder names written to dirnames.in")

if __name__ == '__main__':
    main()

