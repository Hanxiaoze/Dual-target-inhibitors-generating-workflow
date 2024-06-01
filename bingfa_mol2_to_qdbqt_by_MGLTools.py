import os
from concurrent.futures import ProcessPoolExecutor

def prepare_ligand(input_file):
    ligand_name = os.path.splitext(input_file)[0]
    output_file = ligand_name + ".pdbqt"
    command = f"/home/zhou/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l {input_directory_path}/{input_file} -o {output_directory_path}/{output_file}"
    os.system(command)
    print(f"Converted {input_file} to {output_file}")
    return f"Converted {input_file} to {output_file}"

if __name__ == '__main__':
    input_directory_path = "/home/zhou/ZINC_data_bank/ZINC_in-man___3D_add_H___mol2___108314_files_lib"
    output_directory_path = "/home/zhou/ZINC_data_bank/ZINC_in-man___3D_add_H___pdbqt___108314_files_lib"
    input_files = [file for file in os.listdir(input_directory_path) if os.path.isfile(os.path.join(input_directory_path, file))]

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(prepare_ligand, input_files))

        for result in results:
            print(result)

    print("Conversion completed.")
