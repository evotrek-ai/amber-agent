import os
import sys
import subprocess
from typing import Dict, TextIO
from Bio.PDB import PDBParser, PDBIO, Select

class NonSmallMoleculeSelector(Select):
    """选择生物大分子（蛋白/DNA/RNA）、水和离子"""
    def accept_residue(self, residue):
        resname = residue.resname.strip().lower()
        # 判断蛋白质（含CA原子）、DNA/RNA标准残基
        return (residue.has_id("CA") or 
                resname in ["da", "dc", "dg", "dt", "a", "c", "g", "u"] or
                resname == "hoh" or  # 水分子
                (len(residue) == 1 and resname != "hoh"))  # 单原子离子 

class SmallMoleculeSelector(Select):
    """选择小分子配体（排除水、离子和标准残基）"""
    def accept_residue(self, residue):
        resname = residue.resname.strip().lower()
        # 排除标准残基、水和离子
        return not (residue.has_id("CA") or 
                   resname in ["da", "dc", "dg", "dt", "a", "c", "g", "u"] or
                   resname == "hoh" or 
                   (len(residue) == 1 and resname != "hoh"))

def process_pdb_file(input_file, output_prefix):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", input_file)
    first_model = structure[0]  # 只处理第一个MODEL 

    # 生成输出文件名
    output_main = f"{output_prefix}_nosmall.pdb"
    output_ligand = f"{output_prefix}_small.pdb"

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
        

def convert_to_pdbs(input_dir: str, output_dir: str):
    supported_extensions = (
        '.smi', '.sdf', '.mol', '.mol2', '.inchi', '.inchikey', '.xyz', '.cif', '.mcif', '.pqr',
        '.gout', '.gau', '.g03', '.g09', '.g16', '.gamout', '.gms', '.nwo', '.orca', '.gro',
        '.cml', '.fa', '.fasta'
    )
    for f in os.listdir(input_dir):
        if not os.path.isfile(f):
            continue
        if not f.endswith(supported_extensions):
            continue
        basename = os.path.basename(f)
        id, ext = os.path.splitext(basename)
        output = id + '.pdb'
        output_path = os.path.join(output_dir, output)
        print(f'converting {basename} to {output}')
        subprocess.run(['obabel', f, '-o', 'pdb', '-O', output_path])


def split_pdbs(dir: str):
    for f in os.listdir(input_dir):
        if not os.path.isfile(f):
            continue
        if not f.endswith('.pdb'):
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
                # 检查是否为蛋白质残基
                if residue.has_id('CA'):  # 蛋白质残基通常有CA原子
                    continue
                # 检查是否为RNA或DNA残基
                elif residue.resname in ["A", "U", "C", "G", "DA", "DC", "DG", "DT"]:
                    continue
                # 检查是否为水分子
                elif residue.resname == "HOH":
                    continue
                # 其他情况认为是小分子
                else:
                    small_molecules.append(residue.resname)
    
    # 去重并统计数量
    unique_small_molecules = list(set(small_molecules))
    count = len(unique_small_molecules)
    
    return count, unique_small_molecules


def split_ligands(file: str):
    count, molecules = get_small_molecules(file)
    basename = os.path.basename(file)
    id, ext = os.path.splitext(basename)
    with open(file, 'r') as infile:
        outfiles: Dict[str, TextIO] = {}
        try:
            for molecule in molecules:
                f = open(f'{id}_{molecule}.{ext}', 'w')
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
        

def splib_all_ligands(dir: str):
    for f in os.listdir(dir):
        if not os.path.isfile(f):
            continue
        if not f.endswith('_small.pdb'):
            continue
        split_ligands(f)


def merge_pdbs(dir: str):
    subprocess.run(['obabel', '-ipdb', '*_nosmall.pdb', '-O', 'combined_output.pdb'])


def preprocess(input_dir: str, output_dir: str):
    convert_to_pdbs(input_dir, output_dir)
    split_pdbs(output_dir)
    splib_all_ligands(output_dir)
    merge_pdbs(output_dir)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        output_dir = input_dir
    os.makedirs(output_dir, exist_ok=True)
    preprocess(input_dir, output_dir)