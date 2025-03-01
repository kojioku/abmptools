import os
import glob
import subprocess

# dirnames.in ファイルからディレクトリ名を読み込む
with open('dirnames.in', 'r') as file:
    directories = [line.strip() for line in file.readlines()]

# setupmd.log ファイルを開く
with open('setupmd.log', 'w') as log_file:
    # 各ディレクトリに対して処理を実行
    for directory in directories:
        os.chdir(directory)

        # "_ligand_*.pdb" ファイルを検索
        ligand_files = glob.glob("*_ligand_*.pdb")
        if not ligand_files:
            raise FileNotFoundError(f"No '*_ligand_*.pdb' files found in the directory {directory}.")

        # "_receptor_*.pdb" ファイルを検索
        receptor_files = glob.glob("*_receptor_*.pdb")
        if not receptor_files:
            raise FileNotFoundError(f"No '_receptor_*.pdb' files found in the directory {directory}.")

        # 最初のファイルを使用
        ligand_file = ligand_files[0]
        receptor_file = receptor_files[0]

        # acpype を実行
        acpype_command = f"acpype -i {ligand_file} -c bcc -k maxcyc=0"
        subprocess.run(acpype_command, shell=True, check=True)

        # .acpype フォルダ内の .frcmod と .mol2 ファイルを検索
        acpype_dir = f"{ligand_file.split('.')[0]}.acpype"
        frcmod_file = glob.glob(f"{acpype_dir}/*.frcmod")[0]
        mol2_file = glob.glob(f"{acpype_dir}/*_bcc_*.mol2")[0]

        # leaprc ファイルの内容 (complex)
        leaprc_complex_content = f"""
# gengeral
verbosity 2
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p
source leaprc.gaff2
source leaprc.gaff
addPdbResMap {{ {{ "LI+" "LI" }} {{ "NA+" "NA" }} {{ "MG2" "MG" }} {{ "CL-" "CL" }} {{ "K+" "K" }} {{ "RB+" "RB" }} {{ "CS+" "CS" }} }}
addPdbAtomMap {{ {{ "LI+" "LI" }} {{ "NA+" "NA" }} {{ "MG2" "MG" }} {{ "CL-" "CL" }} {{ "K+" "K" }} {{ "RB+" "RB" }} {{ "CS+" "CS" }} }}
# read unparamitrized mol param
loadAmberParams {frcmod_file}
# read mol coord
ligand = loadmol2 {mol2_file}
protein = loadpdb {receptor_file}
complex = combine {{ protein ligand }}
addIons2 complex Na+ 0
addIons2 complex Cl- 0
savePdb complex complex.pdb
solvateBox complex TIP3PBOX 10
saveAmberParm complex complex.prmtop complex.inpcrd
quit
"""

        # leaprc ファイルに内容を書き込む
        with open("leaprc_complex", "w") as file:
            file.write(leaprc_complex_content)

        # tleap を実行する (complex)
        subprocess.run("tleap -f leaprc_complex > leap_complex.log 2>&1", shell=True, check=True)
        # tleap を実行する (ligand)

        # leap.log ファイルから 'Exiting LEaP: ' の行を抽出して setupmd.log に書き出す
        for log_type in ["complex"]:
            log_filename = f"leap_{log_type}.log"
            with open(log_filename, 'r') as leap_log:
                for line in leap_log:
                    if 'Exiting LEaP: ' in line:
                        log_file.write(f"{directory} ({log_type}): {line}")

        # 元のディレクトリに戻る
        os.chdir("..")
