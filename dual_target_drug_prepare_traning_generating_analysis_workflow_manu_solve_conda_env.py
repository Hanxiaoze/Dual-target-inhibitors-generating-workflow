# Author: Zhong-Xing Zhou @ Jilin University
# Time: 2023.09.15

# before perform this script, the path in this script should be modified to your file env 

# the times*.ttf is the word font file to use in the matplotlib, please ensure
# them in to the path ~/anaconda3/envs/3D_Scaffold_03/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/

# config_1.txt and config_2.txt are used for Docking-site-pocket setting, please change according to your two target protein systems

import os 
import urllib.request
import subprocess
import sys
import importlib.util

MGLToolsPckgs = '/home/zhou/MGLTools-1.5.6/MGLToolsPckgs'   # install MGLTools and change the path to your env
Vinna_bin_path = '/home/zhou/autodock_vina_1_1_2_linux_x86/bin'    # install AutoDock Vinna and change the path to your env


# ######
# ###### creat root work-folder for dual-target workflow
work_dir = './dual_target_work_dir'

if os.path.exists(work_dir):
    os.system(f"rm -r {work_dir}")  
os.makedirs(work_dir)

absolute_work_dir_path = os.path.abspath(work_dir)

absolute_Mperformer_path = os.path.abspath('./MPerformer-master')

absolute_cG_SchNet_path = os.path.abspath('./cG-SchNet-main')

os.chdir(f"{absolute_work_dir_path}")



def check_and_install_package(package_name):
    # Check if the package is installed
    if importlib.util.find_spec(package_name) is None:
        print(f"{package_name} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
    else:
        print(f"{package_name} is already installed.")

# Check and install requests if not available
check_and_install_package("requests")

def check_conda_env(env_name):
    # check Conda env existence
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        env_list = result.stdout
        # judge the env whether in the env_list
        if env_name in env_list:
            return True
        else:
            return False
    except Exception as e:
        print(f"Error checking Conda environments: {e}")
        return False

def remove_conda_env(env_name):
    # remove Conda env 
    try:
        subprocess.run(
            ["conda", "env", "remove", "-n", env_name],
            check=True
        )
        print(f"Conda environment '{env_name}' removed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error removing Conda environment '{env_name}': {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


def download_file(url, local_directory, local_filename):
    # download the url and save to the local_directory
    file_path = f"{local_directory}/{local_filename}"    
    
    try:
        urllib.request.urlretrieve(url, file_path)
        print(f"File downloaded and saved to {file_path}")
    except Exception as e:
        print(f"Failed to download file. Error: {e}")


# ######
# ###### prepare conda env

# # download Generating SchNet lib from: https://github.com/atomistic-machine-learning/G-SchNet
# # download c-Generating SchNet lib from: https://github.com/atomistic-machine-learning/cG-SchNet
        
# # if check_conda_env('3D_Scaffold_03'):
# #     remove_conda_env('3D_Scaffold_03')
# #     print(f"Conda environment '3D_Scaffold_03' is removed, now is re-creating it.")
# # else:
# #     print(f"Conda environment '3D_Scaffold_03' does not exist, now is creating it.")



# # os.system(f"conda create -n 3D_Scaffold_03 python=3.8 numpy=1.23.5 pytorch=1.5.1 torchvision cudatoolkit=10.2 ase=3.19.0 openbabel=3.1.1 rdkit=2019.09.2.0 requests matplotlib seaborn scikit-learn -c pytorch -c openbabel -c defaults -c conda-forge")
# # os.system(f"source $(conda info --base)/etc/profile.d/conda.sh")
# # os.system(f"conda run -n 3D_Scaffold_03 pip install 'schnetpack==0.3'")
# # os.system(f"cp ../times*.ttf ~/anaconda3/envs/3D_Scaffold_03/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/")



# ######
# ###### download the ZINC in-man subsets from ZINC20 website using a parallel-downloading script
os.system(f"cp ../*.py .")
os.system(f"python ./ZINC_code_acquire.py")   # write ZINC code in ZINC_code_in-man.txt
os.system(f"python ./sort_ZINC_code.py")      # sort the ZINC code in ZINC_code_in-man.txt
os.system(f"python ./bingfa_pachong_in-man_sdf.py")   # parallel-download the ZINC .sdf file according to the ZINC code
os.system(f"python ./discree_and_depart_2D_or_3D.py")  # file clean to avoid 2D sdf file, judge the sdf file is 3D or 2D
os.system(f"python ./bingfa_pachong_duojincheng_2D_smi.py")   # parallel-download the 2D sdf corresponding SMILES file to convert to 3D sdf in later
os.system(f"conda run -n 3D_Scaffold_03 python ./duojincheng_obabel_smi_to_3D_sdf.py")   # parallel-convert the SMILES to the 3D sdf files
os.system(f"cp ./ZINC_in-man_smi_to_3D_sdf_files/*.sdf ./ZINC_in-man_sdf_files_3D/")    # add the 2D_sdf_to_smi_to_3D_sdf to the origin right 3D_sdf
os.system(f"conda run -n 3D_Scaffold_03 python ./bingfa_sdf_to_mol2.py")   # parallel-convert sdf to mol2 file and calculate the Gasteiger charge needed by AutoDock Vinna
os.system(f"python ./check_nan_and_remove_in_mol2.py")   #  check 'nan' in mol2 files and remove the error files

# ######
# ###### prepare AutoDock Vinna input file (.pdbqt)
os.system(f"conda run -n 3D_Scaffold_03 python ./bingfa_mol2_to_pdbqt_by_obabel.py")   # convert the drug-like molecules to .pqbqt files for Vinna
os.system(f"conda run -n 3D_Scaffold_03 python ./filter_small_55.py")    # optional; control the atom number of training molecules
os.system(f"cp ../target_protein_1.pdb .")
os.system(f"cp ../target_protein_2.pdb .")
os.system(f"conda run -n 3D_Scaffold_03 python ./pdb_A_B_C_D_E_chains_pymol.py")   # optional; handle the target protein PDB residue ID
os.system(f"conda run -n 3D_Scaffold_03 python ./pdb_duplicate_residue_remove.py")   # optional; handle the target protein residue multi-conformation in PDB
os.system(f"{MGLToolsPckgs}/AutoDockTools/Utilities24/prepare_receptor4.py -r ../target_protein_1.pdb")  # convert the target-protein PDB to pdbqt
os.system(f"{MGLToolsPckgs}/AutoDockTools/Utilities24/prepare_receptor4.py -r ../target_protein_2.pdb")


# ######
# ###### perform Vinna docking for two target protein
os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_1")
os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_2")

os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_1/log_file")
os.system(f"cp ../extract_affinity.py {absolute_work_dir_path}/Docking_dir_for_target_protein_1/extract_affinity.py")
os.system(f"cp ../sort_energy.py {absolute_work_dir_path}/Docking_dir_for_target_protein_1/sort_energy.py")
os.system(f"cp ../config_1.txt {absolute_work_dir_path}/Docking_dir_for_target_protein_1/config.txt")
os.system(f"cp ./target_protein_1.pdbqt {absolute_work_dir_path}/Docking_dir_for_target_protein_1/target_protein_1.pdbqt")

os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_2/log_file")
os.system(f"cp ../extract_affinity.py {absolute_work_dir_path}/Docking_dir_for_target_protein_2/extract_affinity.py")
os.system(f"cp ../sort_energy.py {absolute_work_dir_path}/Docking_dir_for_target_protein_2/sort_energy.py")
os.system(f"cp ../config_2.txt {absolute_work_dir_path}/Docking_dir_for_target_protein_2/config.txt")
os.system(f"cp ./target_protein_2.pdbqt {absolute_work_dir_path}/Docking_dir_for_target_protein_2/target_protein_2.pdbqt")

pdbqt_files = []
for filename in os.listdir('./ZINC_in-man_pdbqt/'):
        if filename.endswith('.pdbqt'):
            pdbqt_files.append(filename.split('.')[0])

os.chdir(f"{absolute_work_dir_path}/Docking_dir_for_target_protein_1")
for mol_id in pdbqt_files:
    print('\n')
    os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_1/results/{mol_id}")    
    print(f"Now is Docking for {mol_id}.pdbqt with target protein 1")
    os.system(f"{Vinna_bin_path}/vina --config {absolute_work_dir_path}/Docking_dir_for_target_protein_1/config.txt --ligand ../ZINC_in-man_pdbqt/{mol_id}.pdbqt --out {absolute_work_dir_path}/Docking_dir_for_target_protein_1/results/{mol_id}/out.pdbqt --log {absolute_work_dir_path}/Docking_dir_for_target_protein_1/log_file/{mol_id}_log.txt")

os.chdir(f"{absolute_work_dir_path}/Docking_dir_for_target_protein_2")  
for mol_id in pdbqt_files:
    print('\n')
    os.system(f"mkdir -p {absolute_work_dir_path}/Docking_dir_for_target_protein_2/results/{mol_id}")
    print(f"Now is Docking for {mol_id}.pdbqt with target protein 2")
    os.system(f"{Vinna_bin_path}/vina --config {absolute_work_dir_path}/Docking_dir_for_target_protein_2/config.txt --ligand ../ZINC_in-man_pdbqt/{mol_id}.pdbqt --out {absolute_work_dir_path}/Docking_dir_for_target_protein_2/results/{mol_id}/out.pdbqt --log {absolute_work_dir_path}/Docking_dir_for_target_protein_2/log_file/{mol_id}_log.txt")
  
# ######
# ###### extract the Vinna docking results for each drug-like molecule    
os.chdir(f"{absolute_work_dir_path}/Docking_dir_for_target_protein_1/log_file")
os.system(f"conda run -n 3D_Scaffold_03 python ../extract_affinity.py") 
os.system(f"cp ./affinity.out ..")
os.chdir(f"../")
os.system(f"conda run -n 3D_Scaffold_03 python ./sort_energy.py") 

os.chdir(f"{absolute_work_dir_path}/Docking_dir_for_target_protein_2/log_file")
os.system(f"conda run -n 3D_Scaffold_03 python ../extract_affinity.py") 
os.system(f"cp ./affinity.out ..")
os.chdir(f"../")
os.system(f"conda run -n 3D_Scaffold_03 python ./sort_energy.py") 

os.chdir(f"{absolute_work_dir_path}")

# ######
# ###### combine the two target protein docking results with each ZINC drug-like molecules .xyz file
os.system(f"conda run -n 3D_Scaffold_03 python ./duojincheng_pdbqt_to_xyz_by_obabel.py")  # from the pdbqt to xyz_without_nonpolar_hydrogen to save GPU memory
os.system(f"conda run -n 3D_Scaffold_03 python ./xyz_add_property.py")    # the molecules which not succeed in two docking processes will be filtered out in this step
os.system(f"conda run -n 3D_Scaffold_03 python ./xyz_element_filter.py")  # control the element type of training molecules


# ######
# ###### convert the add-property .xyz file to the ase binary .db file format to accelerate the training
os.system(f"conda run -n 3D_Scaffold_03 python ./generate_3D_Scaffold.py --xyz_path {absolute_work_dir_path}/ZINC_in-man_xyz_with_prop_elem_filter")  


# ######
# ######  pre-calculate the atom-pairs distance and save to accelerate the training
os.system(f"conda run -n 3D_Scaffold_03 python ./zzx_preprocess_dataset.py ./world.db")   # pre-calculate the atom-pairs distance
os.system(f"conda run -n 3D_Scaffold_03 python ./read_stat_npz.py ./world/worldgen_statistics.npz")  # show the preprocess_dataset logfile
os.system(f"conda run -n 3D_Scaffold_03 python ./print_db.py")   # print and show the dataset


# ######
# ###### training neuron network
os.chdir(f"{absolute_cG_SchNet_path}")
if os.path.exists("./data/ZINC_in-man/"):
    os.system(f"rm -r ./data/ZINC_in-man")  
os.makedirs("./data/ZINC_in-man/")

if os.path.exists("./models/ZINC_in-man/"):
    os.system(f"rm -r ./models/ZINC_in-man")  

os.system(f"cp {absolute_work_dir_path}/worldgen.db ./data/ZINC_in-man/qm9gen.db")
os.system(f"cp {absolute_work_dir_path}/world_invalid.txt ./data/ZINC_in-man/qm9_invalid.txt")
os.system(f"conda run -n 3D_Scaffold_03 python {absolute_cG_SchNet_path}/gschnet_cond_script.py train gschnet ./data/ZINC_in-man/ ./models/ZINC_in-man/ --conditioning_json_path ./data/cond_Epro_Mpro.json --split 53000 1838 --cuda") 


# ######
# ###### monitor the training performance metrics
os.chdir(f"{absolute_cG_SchNet_path}/models/ZINC_in-man/log")
os.system(f"cp {absolute_work_dir_path}/plot_log_complex.py .")
os.system(f"conda run -n 3D_Scaffold_03 python ./plot_log_complex.py")

os.system(f"tensorboard --logdir=./")


# ######
# ###### recall the optimal nn weight to generate new molecules
newest_file='generated_3000_8_8.mol_dict'    # this 3000 8 8 should be consistent with the below generating conditions:
os.chdir(f"{absolute_cG_SchNet_path}")
os.system(f"conda run -n 3D_Scaffold_03 python {absolute_cG_SchNet_path}/gschnet_cond_script.py generate gschnet ./models/ZINC_in-man/ 3000 --conditioning \"Epro_affinity -8.0; Mpro_affinity -8.0\" --cuda")



# ######
# ###### prepare conda env for MPerformer nn to predict the bond order according to bond distance, bond angle, bond dihedral

# # download MPerformer lib from: https://github.com/FanmengWang/MPerformer
# # download the dependence Uni-Core from: https://github.com/dptech-corp/Uni-Core/releases

# # os.chdir(f"{absolute_Mperformer_path}")
# # if check_conda_env('MPerformer'):
# #     remove_conda_env('MPerformer')
# #     print(f"Conda environment 'MPerformer' is removed, now is re-creating it.")
# # else:
# #     print(f"Conda environment 'MPerformer' does not exist, now is creating it.")

# # os.system(f"conda create -n MPerformer python=3.9 pytorch=2.0 torchvision cudatoolkit=11.7 ase openbabel -c pytorch -c openbabel -c defaults -c conda-forge")
# # os.system(f"source $(conda info --base)/etc/profile.d/conda.sh")
# # download_file('https://github.com/dptech-corp/Uni-Core/releases/download/0.0.3/unicore-0.0.1+cu117torch2.0.0-cp39-cp39-linux_x86_64.whl', './', 'unicore-0.0.1+cu117torch2.0.0-cp39-cp39-linux_x86_64.whl')
# # os.system(f"conda run -n MPerformer pip install ./unicore-0.0.1+cu117torch2.0.0-cp39-cp39-linux_x86_64.whl")
# # os.system(f"conda run -n MPerformer pip install rdkit-pypi==2021.9.4")
# # os.system(f"conda run -n MPerformer pip install dpdata")
# # os.system(f"conda run -n MPerformer pip install torch==2.0 torchvision torchaudio")
# # os.system(f"conda run -n MPerformer pip install pandas")
# # os.system(f"conda run -n MPerformer pip install scikit-learn")
# # os.system(f"conda run -n MPerformer pip install numpy")


######
###### valid and analysis the quality of the nn generated molecules
os.chdir(f"{absolute_Mperformer_path}")
os.system(f"cp {absolute_cG_SchNet_path}/models/ZINC_in-man/generated/generated.mol_dict {absolute_Mperformer_path}/{newest_file}")
newest_dir=newest_file+"_dir"
if os.path.exists(f"{absolute_Mperformer_path}/{newest_dir}"):
    os.system(f"rm -r {absolute_Mperformer_path}/{newest_dir}") 
os.system(f"mkdir {absolute_Mperformer_path}/{newest_dir}")
os.system(f"cp ./{newest_file} {absolute_Mperformer_path}/{newest_dir}/newest.mol_dict")
os.chdir(f"{absolute_Mperformer_path}/{newest_dir}/")
os.system(f"mkdir ./xyz/")
os.system(f"mkdir ./sdf/")
os.system(f"mkdir ./mol2_Q/")
os.system(f"mkdir ./pdbqt/")
os.system(f"cp ../mol_dict_to_xyz.py .")
os.system(f"cp ../sdf_to_mol2Q_obabel.py .")
os.system(f"cp ../bingfa_mol2_to_pdbqt_by_obabel.py .")
os.system(f"conda run -n 3D_Scaffold_03 python ./mol_dict_to_xyz.py")
os.chdir(f"{absolute_Mperformer_path}/")
os.system(f"conda run -n MPerformer python ./predict.py --filename {absolute_Mperformer_path}/{newest_dir}/xyz --outputs_path {absolute_Mperformer_path}/{newest_dir}/sdf")
os.chdir(f"{absolute_Mperformer_path}/{newest_dir}")
os.system(f"conda run -n 3D_Scaffold_03 python ./sdf_to_mol2Q_obabel.py")
os.system(f"conda run -n 3D_Scaffold_03 python ./bingfa_mol2_to_pdbqt_by_obabel.py")
os.chdir(f"{absolute_work_dir_path}")
if os.path.exists(f"{absolute_work_dir_path}/{newest_dir}"):
    os.system(f"rm -r {absolute_work_dir_path}/{newest_dir}") 
os.system(f"mkdir ./{newest_dir}")
os.chdir(f"{newest_dir}")
os.system(f"cp {absolute_Mperformer_path}/{newest_file} .")
os.system(f"cp -r {absolute_Mperformer_path}/{newest_dir}/xyz .")
os.system(f"cp -r {absolute_Mperformer_path}/{newest_dir}/sdf .")
os.system(f"cp -r {absolute_Mperformer_path}/{newest_dir}/mol2_Q .")
os.system(f"cp -r {absolute_Mperformer_path}/{newest_dir}/pdbqt .")
os.system(f"cp ../../E_M_6W_affinity.txt .")    # your training molecular id with two docking affinities
os.system(f"cp -r ../Docking_dir_for_target_protein_1 .")
os.system(f"rm -r ./Docking_dir_for_target_protein_1/results/*")
os.system(f"rm ./Docking_dir_for_target_protein_1/log_file/*")
os.system(f"rm ./Docking_dir_for_target_protein_1/affinity.out")
os.system(f"rm ./Docking_dir_for_target_protein_1/sorted_final_energy.txt")
os.system(f"cp -r ../Docking_dir_for_target_protein_2 .")
os.system(f"rm -r ./Docking_dir_for_target_protein_2/results/*")
os.system(f"rm ./Docking_dir_for_target_protein_2/log_file/*")
os.system(f"rm ./Docking_dir_for_target_protein_2/affinity.out")
os.system(f"rm ./Docking_dir_for_target_protein_2/sorted_final_energy.txt")
os.chdir(f"./Docking_dir_for_target_protein_1")
pdbqt_files = []
for filename in os.listdir('../pdbqt/'):
        if filename.endswith('.pdbqt'):
            pdbqt_files.append(filename.split('.')[0])
for mol_id in pdbqt_files:
    print('\n') 
    print(f"Now is Docking for {mol_id}.pdbqt with target protein 1")   
    os.system(f"mkdir ./results/{mol_id}")
    os.system(f"{Vinna_bin_path}/vina --config ./config.txt --ligand ../pdbqt/{mol_id}.pdbqt --out ./results/{mol_id}/out.pdbqt --log ./log_file/{mol_id}_log.txt")
os.chdir(f"../Docking_dir_for_target_protein_2")
for mol_id in pdbqt_files:
    print('\n')
    print(f"Now is Docking for {mol_id}.pdbqt with target protein 2")    
    os.system(f"mkdir ./results/{mol_id}")
    os.system(f"{Vinna_bin_path}/vina --config ./config.txt --ligand ../pdbqt/{mol_id}.pdbqt --out ./results/{mol_id}/out.pdbqt --log ./log_file/{mol_id}_log.txt")
os.chdir(f"./log_file")
os.system(f"python ../extract_affinity.py")
os.system(f"cp ./affinity.out ..")
os.chdir(f"..")
os.system(f"python sort_energy.py")
os.chdir(f"../Docking_dir_for_target_protein_1")
os.chdir(f"./log_file")
os.system(f"python ../extract_affinity.py")
os.system(f"cp ./affinity.out ..")
os.chdir(f"..")
os.system(f"python sort_energy.py")
os.chdir(f"..")
os.system(f"cp ../E_M_combine.py .")
os.system(f"cp ../sort_sum_energy.py .")
os.system(f"cp ../warm_up_filter.py .")
os.system(f"cp ../plot_1d_prob_hebing.py .")
os.system(f"cp ../plot_3d_surface_prob.py .")
os.system(f"cp ../plot_3d_surface_prob_diff.py .")
os.system(f"cp ../plot_1d_freq.py .")
os.system(f"cp ../plot_3d_surface_freq.py .")
os.system(f"cp ../statistic_and+0.py .")
os.system(f"cp ../statistic_and+1.py .")
os.system(f"cp ../statistic_and+2.py .")
os.system(f"cp ../statistic_sum.py .")
os.system(f"cp ../statistic_func_group_gro.py .")
os.system(f"cp ../statistic_func_group_mol.py .")
os.system(f"cp ../finger_print_plot.py .")
os.system(f"cp ../HDBSCAN.py .")
os.system(f"cp ../Plot_good_mol.py .")
os.system(f"python E_M_combine.py")
os.system(f"python warm_up_filter.py  0.66")    # remove the early part generated molecules to equilibrium the generator, 0.66 is warm_up value
os.system(f"python sort_sum_energy.py")
os.system(f"conda run -n 3D_Scaffold_03 python plot_1d_prob_hebing.py 8 8")    # 8 8 are the generating conditions
os.system(f"conda run -n 3D_Scaffold_03 python plot_3d_surface_prob.py 8 8")   # 8 8 are the generating conditions
os.system(f"conda run -n 3D_Scaffold_03 python plot_3d_surface_prob_diff.py")
os.system(f"conda run -n 3D_Scaffold_03 python plot_1d_freq.py")
os.system(f"conda run -n 3D_Scaffold_03 python plot_3d_surface_freq.py")
os.system(f"conda run -n 3D_Scaffold_03 python statistic_and+0.py")
os.system(f"conda run -n 3D_Scaffold_03 python statistic_and+1.py")
os.system(f"conda run -n 3D_Scaffold_03 python statistic_and+2.py")
os.system(f"conda run -n 3D_Scaffold_03 python statistic_sum.py")
os.system(f"cp ../bingfa_mol2_to_smi.py .")
os.system(f"conda run -n 3D_Scaffold_03 python bingfa_mol2_to_smi.py")
os.system(f"cp ../smi_to_warm_up.py .")
os.system(f"conda run -n 3D_Scaffold_03 python smi_to_warm_up.py 0.66")    #  0.66 is warm_up value
os.system(f"conda run -n 3D_Scaffold_03 python statistic_func_group_gro.py")
os.system(f"conda run -n 3D_Scaffold_03 python statistic_func_group_mol.py")
os.system(f"conda run -n 3D_Scaffold_03 python finger_print_plot.py")
os.system(f"conda run -n 3D_Scaffold_03 python HDBSCAN.py")
os.system(f"conda run -n 3D_Scaffold_03 python Plot_good_mol.py")




