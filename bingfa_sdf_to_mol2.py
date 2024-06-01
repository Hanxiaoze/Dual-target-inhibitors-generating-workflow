import os
import glob
from concurrent.futures import ProcessPoolExecutor
from openbabel import openbabel

# 设置输入SDF文件夹和输出MOL2文件夹
input_folder = './ZINC_in-man_sdf_files_3D'
output_folder = './ZINC_in-man_mol2'
os.makedirs(output_folder, exist_ok=True)

def convert_sdf_to_mol2(sdf_file):
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats('sdf', 'mol2')
    
    mol = openbabel.OBMol()
    ob_conversion.ReadFile(mol, sdf_file)
    
    mol2_file = os.path.join(output_folder, os.path.basename(sdf_file).replace('.sdf', '.mol2'))
    ob_conversion.WriteFile(mol, mol2_file)
    print(f"Converted {sdf_file} to MOL2 format.")
    
    return f"Converted {sdf_file} to MOL2 format."

if __name__ == '__main__':
    sdf_files = glob.glob(os.path.join(input_folder, '*.sdf'))
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(convert_sdf_to_mol2, sdf_files))
        
        for result in results:
            print(result)

    print("Conversion completed.")
