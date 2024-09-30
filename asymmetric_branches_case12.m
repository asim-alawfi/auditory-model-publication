%% Computation of non-symmetric POs for different values of tone duration t_\mathrm{d} and delay (D)  
clear
    base=[pwd(),'\..\ddebiftool_snapshot_23October2022\'];
    base2=[pwd(),'\..\Supporting_function\'];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco']);
addpath([base2,'Supporting_functions'])
%%
load('br_crossing_threshold_try2.mat')
%% Now, branching off from the branch of symmetric solutions to a branch of non-symmetric solutions
% (a note form me:period doubling bifurcation are expected to be detected).
branch2=br_symmetry_wbifs(2);
S2=dde_lincond_struct(size(branch2.point(1).profile,1),'profile',...
    'shift',[1,2],'condprojmat',[1,0,0,0,0,0],'condprojint',[0,0.5]);
plm={@(p)p.parameter(branch2.parameter.free),@(p)dde_psol_lincond(p,S2)};
%%
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.05; in.df,0.01; 0,0.02],'print_residual_info',1};
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_asy,asym_brs2,suc_v]=SetupPsol(funcs_audi,br_symmetry_wbifs(2), p_bif(2,2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,parbds{:},...
  nspoev1args{:},'plot_measure',plm);
%%
figure(334)
hold on
asym_brs2=br_contn(funcs_asy,asym_brs2,276);
asym_brs2=br_rvers(asym_brs2);
asym_brs2=br_contn(funcs_asy,asym_brs2,276);
%%
[asym_brs2_wbifs,nunst_td,bif_td,p_bif_td]=MonitorChange(funcs_asy,asym_brs2,...
    'range',2:length(asym_brs2.point),'printlevel',1,'print_residual_info',1,...
    'min_iterations',5);
%% branching off from the branch of symmetric solutions to a branch of non-symmetric solutions
% (a note form me:period doubling bifurcation are expected to be detected).
branch_d=br_symmetry_wbifs(1);
Sd=dde_lincond_struct(size(branch_d.point(1).profile,1),'profile',...
    'shift',[1,2],'condprojmat',[1,0,0,0,0,0],'condprojint',[0,0.5]);
plm_d={@(p)p.parameter(branch_d.parameter.free),@(p)dde_psol_lincond(p,Sd)};
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.05; in.df,0.01; 0,0.02],'print_residual_info',1};
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_asy,asym_brs1,suc_v1]=SetupPsol(funcs_audi,br_symmetry_wbifs(1), p_bif(2,1),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,parbds{:},...
  nspoev1args{:},'plot_measure',plm_d);
%%
figure(333)
clf;
hold on
asym_brs1=br_contn(funcs_asy,asym_brs1,340);
[asym_brs1_wbifs,nunst_1,bif_1,p_bif_1]=MonitorChange(funcs_asy,asym_brs1,...
    'range',2:length(asym_brs1.point),'printlevel',1,'print_residual_info',1,...
    'min_iterations',5);
%%
asym_brs1_wbifs=br_remove_extracolumns(asym_brs1_wbifs);
asym_brs1=br_remove_extracolumns(asym_brs1);
asym_brs2=br_remove_extracolumns(asym_brs2);
asym_brs2_wbifs=br_remove_extracolumns(asym_brs2_wbifs);
%%
 save('asymmetric_branches_case12.mat')