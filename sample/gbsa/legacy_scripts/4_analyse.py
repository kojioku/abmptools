import csv
import re
import os
import matplotlib.pyplot as plt
from datetime import datetime

# dirnames.in ファイル名
dirnames_file = 'dirnames.in'

# データを格納するリスト
data = []

# dirnames.in ファイルを読み込む
with open(dirnames_file, 'r') as file:
    directories = [line.strip() for line in file.readlines()]

# 各フォルダ内のログファイルを読み込む
for directory in directories:
    row = [directory]
    for log_type in ['complex', 'ligand', 'receptor']:
        log_file = os.path.join(directory, f'{log_type}.log')
        with open(log_file, 'r') as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if '[STEP4] Compute Single Point Energy for Molecules' in line:
                # 4行下のデータを取得
                target_line = lines[i + 4].strip()
                # 正規表現で数値を抽出（整数の0も含む）
                values = re.findall(r'-?\d+\.\d+|(?<!\d)0(?!\d)', target_line)
                # 第2項目と第9項目を取得してリストに追加
                energy = float(values[1])
                solvation = float(values[8])
                egas = energy - solvation
                row.extend([egas, solvation])
                break
    data.append(row)

# 現在の日時を取得
now = datetime.now()
timestamp = now.strftime("%Y%m%d-%H-%M-%S")

# 出力ファイル名
output_file = f'analysis_results_{timestamp}.csv'

# CSVファイルに書き込む
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    # ヘッダーを書き込む
    writer.writerow(['Folder', 'Complex Egas', 'Complex SOLVLATION', 'Ligand Egas', 'Ligand SOLVLATION', 'Receptor Egas', 'Receptor SOLVLATION', 'Total dGbind'])
    
    # データを書き込む
    for row in data:
        complex_egas = row[1]
        complex_solvation = row[2]
        ligand_egas = row[3]
        ligand_solvation = row[4]
        receptor_egas = row[5]
        receptor_solvation = row[6]
        total_dg_bind = complex_egas + complex_solvation - (ligand_egas + ligand_solvation + receptor_egas + receptor_solvation)
        row.append(total_dg_bind)
        writer.writerow(row)

print(f"Results have been written to {output_file}")

# 棒グラフを描画
folders = [row[0] for row in data]
dg_bind_values = [row[7] for row in data]

plt.figure(figsize=(10, 6))
plt.bar(folders, dg_bind_values, color='blue')
plt.xlabel('Folder')
plt.ylabel('Total dGbind')
plt.title('Total dGbind for each Folder')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# グラフをPNGファイルとして保存
plt.savefig(f'dg_bind_plot_{timestamp}.png')
plt.show()

