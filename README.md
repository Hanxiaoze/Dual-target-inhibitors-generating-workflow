# Dual-target-inhibitors-generating-workflow

![Dual-target drug generator workflow](https://github.com/Hanxiaoze/Dual-target-inhibitors-generating-workflow/assets/50012246/82870225-e263-449a-93e7-9d0d9dcfe094)

You need download the whole main folder "Dual-target-inhibitors-generating-workflow" (the zip file is not need, which is the datafile for the paper reproduce) with its child folders, all the folders and files need to keep the relative path like this Git reporitory. The Anaconda, MGLTools-1.5.6, and Autodock Vina are needed as dependencies.
The script dual_target_drug_prepare_traning_generating_analysis_workflow.py is the main interface, it will automatically solve and creat Conda env and call other script to run and the whole workflow. The absolute paths of MGLToolsPckgs and Vinna_bin_path are need to change with your installing path. If you have solved the Conda env manually, you can just use the script dual_target_drug_prepare_traning_generating_analysis_workflow_manu_solve_conda_env.py to run the whole workflow.

1. Reference the Generating SchNet lib from: https://github.com/atomistic-machine-learning/G-SchNet
2. Reference the c-Generating SchNet lib from: https://github.com/atomistic-machine-learning/cG-SchNet
3. Reference the MPerformer lib from: https://github.com/FanmengWang/MPerformer
4. Reference the Uni-Core from: https://github.com/dptech-corp/Uni-Core/releases

