"""Build a PVP (polyvinylpyrrolidone) 5-mer 3D SDF for the amorphous sample.

Each repeat unit carries a 2-pyrrolidinone (gamma-lactam) side group whose
nitrogen is tertiary (no N-H) -> an amide C=O H-bond ACCEPTOR with no donor.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

n = 5
smiles = "CC(N1CCCC1=O)" + "CC(N1CCCC1=O)" * (n - 1) + "C"
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToMolFile(mol, "pvp_5mer.sdf")
print("wrote pvp_5mer.sdf:", mol.GetNumAtoms(), "atoms")
