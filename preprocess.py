import glob
import os
import shutil
import sys
import subprocess
from typing import Dict, TextIO
from Bio.PDB import PDBParser, PDBIO, Select, is_aa

# --- Define Sets for Standard Base Residue Names (lowercase) ---
# These sets now contain only the core base names
RNA_BASE_NAMES_LOWER = {'a', 'c', 'g', 'u', 'ra', 'rc', 'rg', 'ru'}
DNA_BASE_NAMES_LOWER = {'da', 'dc', 'dg', 'dt'}
WATER_RESIDUES_LOWER = {"hoh", "wat", "h2o", "tip", "tip3", "tip3p", "tip4p", "spc"}
# ----------------------------------------------------------

# --- Helper Functions for Nucleic Acid Checks ---
def is_dna_residue(resname_lower, base_dna_names):
    """Checks if a residue name matches a base DNA name optionally followed by digits."""
    for base_name in base_dna_names:
        if resname_lower.startswith(base_name):
            suffix = resname_lower[len(base_name):]
            # Check if the suffix is empty OR consists only of digits
            if not suffix or suffix.isdigit():
                return True
    return False

def is_rna_residue(resname_lower, base_rna_names):
    """Checks if a residue name matches a base RNA name optionally followed by digits."""
    for base_name in base_rna_names:
        if resname_lower.startswith(base_name):
            suffix = resname_lower[len(base_name):]
            # Check if the suffix is empty OR consists only of digits
            if not suffix or suffix.isdigit():
                return True
    return False
# ------------------------------------------------


class NonSmallMoleculeSelector(Select):
    """Selects biological macromolecules (Protein/DNA/RNA), water, and ions."""
    def accept_residue(self, residue):
        resname_lower: str = residue.resname.strip().lower()
        resname_upper: str = residue.resname.strip()

        is_protein = is_aa(resname_upper, standard=True)

        # Use helper functions for nucleic acid checks
        is_dna_acid = is_dna_residue(resname_lower, DNA_BASE_NAMES_LOWER)
        is_rna_acid = is_rna_residue(resname_lower, RNA_BASE_NAMES_LOWER)

        is_water = resname_lower in WATER_RESIDUES_LOWER
        is_ion = (len(residue) == 1 and not is_water)

        return is_protein or is_rna_acid or is_dna_acid or is_water or is_ion

class SmallMoleculeSelector(Select):
    """Selects small molecule ligands (excluding water, ions, and standard protein/NA)."""
    def accept_residue(self, residue):
        resname_lower: str = residue.resname.strip().lower()
        resname_upper: str = residue.resname.strip()

        is_protein = is_aa(resname_upper, standard=True)

        # Use helper functions for nucleic acid checks
        is_dna_acid = is_dna_residue(resname_lower, DNA_BASE_NAMES_LOWER)
        is_rna_acid = is_rna_residue(resname_lower, RNA_BASE_NAMES_LOWER)

        is_water = resname_lower in WATER_RESIDUES_LOWER
        is_ion = (len(residue) == 1 and not is_water)

        return not (is_protein or is_rna_acid or is_dna_acid or is_water or is_ion)

def process_pdb_file(input_file, output_prefix):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", input_file)
    first_model = structure[0]  # 只处理第一个MODEL 

    # 生成输出文件名
    output_main = f"{output_prefix}.nosmall.pdb"
    output_ligand = f"{output_prefix}.small.pdb"

    # 检查是否存在小分子
    has_ligand = any(SmallMoleculeSelector().accept_residue(res) 
                    for res in first_model.get_residues())
    has_nonsmall = any(NonSmallMoleculeSelector().accept_residue(res) 
                    for res in first_model.get_residues())
    # 保存主文件（必生成）
    io = PDBIO()
    io.set_structure(first_model)
    
    if has_nonsmall:
        io.save(output_main, NonSmallMoleculeSelector())
        print(f"bioMol filename: {output_main}")
    else:
        print("为检测到大分子")

    # 仅当存在小分子时生成配体文件
    if has_ligand:
        io.save(output_ligand, SmallMoleculeSelector())
        print(f"chemMol filename: {output_ligand}")
    else:
        print("未检测到小分子配体")
        

def convert_to_pdbs(input_dir: str):
    supported_extensions = (
        '.smi', '.sdf', '.mol', '.mol2', '.inchi', '.inchikey', '.xyz', '.cif', '.mcif', '.pqr',
        '.gout', '.gau', '.g03', '.g09', '.g16', '.gamout', '.gms', '.nwo', '.orca', '.gro',
        '.cml', '.fa', '.fasta'
    )
    for f in os.listdir(input_dir):
        f = os.path.join(input_dir, f)
        if not os.path.isfile(f):
            continue
        if f.endswith('.pdb') and not input_dir == os.getcwd():
            basename = os.path.basename(f)
            shutil.copy(f, basename)
            continue
        if not f.endswith(supported_extensions):
            continue
        basename = os.path.basename(f)
        id, ext = os.path.splitext(basename)
        output = id + '.pdb'
        print(f'converting {basename} to {output}')
        subprocess.run(['obabel', f, '-o', 'pdb', '-O', output])


def split_pdbs():
    for f in os.listdir('.'):
        if not os.path.isfile(f):
            continue
        if not f.endswith('.pdb'):
            continue
        if f.endswith(".small.pdb"):
            continue
        if f.endswith(".nosmall.pdb"):
            continue
        if f.endswith(".ligand.pdb"):
            continue
        id, ext = os.path.splitext(os.path.basename(f))
        process_pdb_file(f, id)


def get_small_molecules(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)
    
    small_molecules = []
    
    # 遍历结构中的每个模型、链和残基
    for model in structure:
        for chain in model:
            for residue in chain:
                resname_lower = residue.resname.lower()
                resname_upper = residue.resname.upper()
                is_protein = is_aa(resname_upper, standard=True)
                if is_protein:
                    continue

                is_dna_acid = is_dna_residue(resname_lower, DNA_BASE_NAMES_LOWER)
                if is_dna_acid:
                    continue

                is_rna_acid = is_rna_residue(resname_lower, RNA_BASE_NAMES_LOWER)
                if is_rna_acid:
                    continue

                is_water = resname_lower in WATER_RESIDUES_LOWER
                if is_water:
                    continue
                                
                # 其他情况认为是小分子
                small_molecules.append(residue.resname)
    
    # 去重并统计数量
    unique_small_molecules = list(set(small_molecules))
    count = len(unique_small_molecules)
    
    return count, unique_small_molecules


def split_ligands(file: str):
    count, molecules = get_small_molecules(file)
    basename = os.path.basename(file)
    id, ext = os.path.splitext(basename)
    id = id[:-len(".small")]
    with open(file, 'r') as infile:
        outfiles: Dict[str, TextIO] = {}
        try:
            for molecule in molecules:
                f = open(f'{id}.{molecule}.ligand{ext}', 'w')
                outfiles[molecule] = f
            
            # Process the input and write to all output files
            for line in infile:
                for molecule in outfiles:
                    if molecule in line:
                        outfiles[molecule].write(line)
                        break

        finally:
            # Make sure to close all output files
            for molecule in outfiles:
                outfiles[molecule].close()

def splib_all_ligands():
    for f in os.listdir('.'):
        if not os.path.isfile(f):
            continue
        if not f.endswith('.small.pdb'):
            continue
        split_ligands(f)


def merge_pdbs():
    # Find all PDB files in the directory
    pdb_files = sorted(glob.glob("*.nosmall.pdb"))
    if not pdb_files:
        return
    output_file = 'combined.nosmall.pdb'
    with open(output_file, 'w') as out_f:
        for i, pdb_file in enumerate(pdb_files):
            # Read the content directly for more control
            with open(pdb_file, 'r') as in_f:
                lines = in_f.readlines()
                
            # Write all lines except END if not the last file
            for line in lines:
                if i < len(pdb_files) - 1 and line.strip() == "END":
                    out_f.write("TER\n")
                else:
                    out_f.write(line)
            
            # Ensure there's a separator between files
            if i < len(pdb_files) - 1 and not lines[-1].strip() in ["END", "TER"]:
                out_f.write("TER\n")

def preprocess(input_dir: str):
    convert_to_pdbs(input_dir)
    split_pdbs()
    splib_all_ligands()
    merge_pdbs()

if __name__ == "__main__":
    input_dir = sys.argv[1]
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        output_dir = input_dir
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)
    preprocess(input_dir)