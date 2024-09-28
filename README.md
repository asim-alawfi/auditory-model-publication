

These scripts perform bifurcation computations for the auditory model (experimental data case). They should be run in sequence, as each script saves its results using a `save(...)` command, and the next script loads those results using a `load(...)` command.

After running all the scripts, you can generate the results (figures) by executing the script `Plot_auditory_model_original_case_submissionversion.m`, which will create the necessary figures.

**Essential steps to run the scripts:**

1) Use the same version of **DDE-Biftool** that was used in the bifurcation analysis for the auditory model. The **DDE-Biftool** package can be obtained from: https://github.com/asim-alawfi/DDE-Biftool-Package.git
   
2) Download the repository folder called **Matlab-Supporting-Functions-for-Auditory-Model** from: https://github.com/asim-alawfi/Matlab-Supporting-Functions-for-Auditory-Model.git

3) The first script that must be run first is `branch_of_sympos_original_case_and_threshold_crossing.m`. "**Note**: At the beginning of this script, the `addpath` is set here and in all other scripts as follows:

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
4) After running the script in step (3) run the other scripts in the following order: 
`branching_towrds_asymmetric_sols_original_case.m`, `sym_breaking_original_case`, `tracking_threshold_crossing_asymmetric`, `identifying_Sols_in_different_regions.m`, and `Plot_auditory_model_original_case_submissionversion.m`


