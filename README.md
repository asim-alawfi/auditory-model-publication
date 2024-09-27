

These scripts perform bifurcation computations for the auditory model (based on experimental data). They should be run in sequence, as each script saves its results using a `save(...)` command, and the next script loads those results using a `load(...)` command.

After running all the scripts, you can generate the results (figures) by executing the script `Plot_auditory_model_original_case_submissionversion.m`, which will create the necessary figures.

**Essential steps to run the scripts:**

1) It is recommended to use the same version of DDE-Biftool that was used for the computations. You can download DDE-Biftool from: ... 
   
2) Download the repository folder called **Matlab-Supporting-Functions-for-Auditory-Model** from: . 

3) The first script you need to run is `branch_of_sympos_original_case_and_threshold_crossing.m`. Before running it, make sure to check and set the MATLAB path. The path in all these scripts ia: 

   ```matlab
   base = [pwd(), '\..\ddebiftool_snapshot_23October2022\'];
   base2 = [pwd(), '\..\Supporting_function\'];
   addpath([base, 'ddebiftool'], ...
           [base, 'ddebiftool_extra_psol'], ...
           [base, 'ddebiftool_utilities'], ...
           [base, 'ddebiftool_extra_rotsym'], ...
           [base, 'ddebiftool_extra_nmfm'], ...
           [base, 'ddebiftool_extra_symbolic'], ...
           [base, 'ddebiftool_coco']);
   addpath([base2, 'Supporting_functions']);
   ```

   - The folder **ddebiftool_snapshot_23October2022** should be obtained from step (1).
   - The folder **Supporting_function** should be obtained from step (2).

