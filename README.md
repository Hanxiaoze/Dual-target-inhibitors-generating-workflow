# Dual-target-inhibitors-generating-workflow

![Dual-target drug generator workflow](https://github.com/Hanxiaoze/Dual-target-inhibitors-generating-workflow/assets/50012246/82870225-e263-449a-93e7-9d0d9dcfe094)


Dependencies
------------
You need download the whole main folder "Dual-target-inhibitors-generating-workflow" (the zip file is not need, which is the datafile for the paper reproducing) with its child folders, all the folders and files need to keep the relative path like this Git reporitory. The Anaconda, MGLTools-1.5.6, and Autodock Vina are needed as dependencies.



Quick Start
------------
The script dual_target_drug_prepare_traning_generating_analysis_workflow.py is the main interface, it will automatically solve and creat Conda env and call other script to run and the whole workflow. The absolute paths of MGLToolsPckgs and Vinna_bin_path are need to change with your installing path. 

If you have solved the Conda env manually, you can just use the script dual_target_drug_prepare_traning_generating_analysis_workflow_manu_solve_conda_env.py to run the whole workflow.

The PDB files target_protein_1.pdb and target_protein_2.pdb are the 3D structures of the two target-proteins which are used as docking receptors, they should be changed to your target-proteins PDB files according to your researching project. The files config_1.txt and config_2.txt are two docking configuration files for target-protein-1 and target-protein-2 respectively, they definite the docking pocket sites of two target-proteins. The files times.ttf, timesi.ttf, timesbi.ttf, timesbd.ttf are the font files to set the matplotlib figure plotting font, the need to be placed to the right place of the corresponding Conda env (~/anaconda3/envs/3D_Scaffold_test/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/).



Reference
--------
1. Reference the Generating SchNet lib from: https://github.com/atomistic-machine-learning/G-SchNet
2. Reference the c-Generating SchNet lib from: https://github.com/atomistic-machine-learning/cG-SchNet
3. Reference the MPerformer lib from: https://github.com/FanmengWang/MPerformer
4. Reference the Uni-Core from: https://github.com/dptech-corp/Uni-Core/releases

