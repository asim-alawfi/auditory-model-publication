

# Auditory model bifurcation analysis with DDE-Biftool 

This repository provides MATLAB code for performing bifurcation analysis on an auditory model using **DDE-Biftool**.

## Requirements

Add the following folders to your MATLAB path:

```matlab
base = [pwd(), filesep(), 'ddebiftool_snapshot_23October2022', filesep()];
addpath([base, 'ddebiftool'], ...
        [base, 'ddebiftool_extra_psol'], ...
        [base, 'ddebiftool_utilities'], ...
        [base, 'ddebiftool_extra_rotsym'], ...
        [base, 'ddebiftool_extra_nmfm'], ...
        [base, 'ddebiftool_extra_symbolic'], ...
        [base, 'ddebiftool_coco'], ...
        [base, 'Supporting_functions']);
```

**Note**: 
1. It is recommended to use the DDE-Biftool version uploaded in `ddebiftool_snapshot_23October2022` for consistency, as this was the version used for the original computations.
2. The `Supporting_functions` folder contains essential functions for the auditory model, such as the right-hand side equations, the sigmoid function, and others.
3. **You don't need to run any scripts from the `ddebiftool_snapshot_23October2022` folder; just ensure they are in your MATLAB path**.

## Running the Scripts

The scripts should be run in a specific sequence to ensure proper execution. This order is based on how results are saved using `save(...)` and subsequently loaded using `load(...)` across scripts.

**Recommended order:**

### Set 1: Bifurcation analysis (based on experimental data) 
1. `branch_of_sympos_original_case_and_threshold_crossing.m`
2. `branching_towrds_asymmetric_sols_original_case.m`
3. `sym_breaking_original_case.m`
4. `tracking_threshold_crossing_asymmetric.m`
5. `identifying_Sols_in_different_regions.m`
6. `plot_results_base_case_publication.m`

### Set 2: One-parameter computations (for tone duration and delay) 
1. `br_of_POs_different_values_TD_D_in_loop.m`
2. `asymmetric_branches_case12.m`

### Set 3: (rp, df)-parameter space computation (case TD=0.05) 
1. `touching_theta_case_TD.m`
2. `symmetry_breaking_for_case_TD.m`

### Set 4: (rp, df)-parameter space computation (case, TD = D = 0.05) 
1. `touching_theta_case_TD_and_D.m`
2. `symmetry_breaking_for_case_TD_and_D.m`

### Set 5: (rp, df)-parameter space computation (case D = 0.05) 
1. `touching_theta_case_D.m`
2. `symmetry_breaking_for_case_D.m`
### 
1. plot_two_parameter_bif_results_for_other_cases.m

