

# Auditory Model Bifurcation Analysis using DDE-Biftool

This repository contains MATLAB code for performing bifurcation computations in an auditory model using **DDE-Biftool**.

## Requirements

Ensure that the following folders are added to your MATLAB path:

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

### Note: 
You **do not** need to run any scripts from the `ddebiftool_snapshot` or `Supporting_function` folders. They simply need to be in the MATLAB path.

## Script Execution

The scripts are grouped into **sets**, and each set must be run in sequence. Each script saves its results with the `save(...)` command, which are loaded by subsequent scripts using the `load(...)` command.

### Set 1: Bifurcation Analysis Based on Experimental Data
Run the following scripts in order:
1. `branch_of_sympos_original_case_and_threshold_crossing.m`
2. `branching_towrds_asymmetric_sols_original_case.m`
3. `sym_breaking_original_case.m`
4. `tracking_threshold_crossing_asymmetric.m`
5. `identifying_Sols_in_different_regions.m`
6. `Plot_auditory_model_original_case_submissionversion.m`

### Set 2: One-Parameter Computation for Varying Delay and Tone Duration
1. `br_of_POs_different_values_TD_D_in_loop.m`
2. `asymmetric_branches_case12.m`
3. `plot_one_parameter_for_all_caes.m`

### Set 3: Computation in (rp, df)-Parameter Space (Large TD)
1. `touching_theta_case_TD.m`
2. `symmetry_breaking_for_case_TD.m`

### Set 4: Computation in (rp, df)-Parameter Space (D = TD = 0.05)
1. `touching_theta_case_TD_and_D.m`
2. `symmetry_breaking_for_case_TD_and_D.m`

### Set 5: Computation in (rp, df)-Parameter Space (Large D = 0.05)
1. `touching_theta_case_D.m`
2. `symmetry_breaking_for_case_D.m`

