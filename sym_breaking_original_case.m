%% Original Case: Computation tracking branches with touching theta
clear;
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
 load('branch_of_sympos_original_case_and_threshold_crossing.mat')
 %%  Parameters bound and symmetry condition
parbd={'min_bound',[in.PR,1;in.df,0.01],'max_bound',[in.PR,41; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.01; 0,0.1],'step',0.05,'print_residual_info',1};
condS=eye(xdim);
% We enforce the symmetry on solution $x$ (not on the eigenvectors) 
pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.4,2)','stateproj',condS,'condprojmat',condS);
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
[funcs_sym_breaking_org,Symbk_org_br,succ]=SetupPOEV1(funcs_s2,po2_symmetry_wbifs,po2_symmetry_bifind(1),...
    'contpar',[in.PR,in.df],'dir',in.df,parbd{:},...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
%%
figure(8)
clf;
hold on
Symbk_org_br=br_contn(funcs_sym_breaking_org,Symbk_org_br,450);
Symbk_org_br=br_rvers(Symbk_org_br);
Symbk_org_br=br_contn(funcs_sym_breaking_org,Symbk_org_br,750);
Symbk_org_br=br_remove_extracolumns(Symbk_org_br);
%% Compuate stability 
[Symbk_org_br_with_stab,nunst_bk,dom,triv_defect_bk]=br_stabl(funcs_sym_breaking_org,Symbk_org_br,0,0);
%% Plot results 
 rp_symbk1=arrayfun(@(x)x.parameter(in.PR),Symbk_org_br_with_stab.point);
 df_symbk1=arrayfun(@(x)x.parameter(in.df),Symbk_org_br_with_stab.point);
figure(10)
clf; 
 hold on 
plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b.-',...
    rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
grid on
%%
save('sym_breaking_original_case.mat')
