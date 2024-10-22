clear;
   base=[pwd(),filesep(),'ddebiftool_snapshot_23October2022',filesep()];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco'],...
            [base,'Supporting_functions']);
%% 
load('tracking_threshold_crossing_asymmetric_try2.mat')
%%
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,30; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.01; 0,0.1],'print_residual_info',1};
[~,it]=min(abs(df_x1-0.2));
[po3_branch,sucp]=ChangeBranchParameters(funcs_s2,br_df,it,...
      'degree',6,'intervals',120,'contpar',in.PR,'step',1e-3,'correc',true,parbds{:},...
      'extra_condition',true,'phase_condition',0,...
      'eigmatrix', 'sparse','matrix', 'sparse');
%
po3_branch=br_rvers(po3_branch);
figure(4)
clf;
hold on
po3_branch=br_contn(funcs_s2,po3_branch,700);
po3_branch=br_rvers(po3_branch);
po3_branch=br_contn(funcs_s2,po3_branch,1000);
%%
unstab=GetStability(po3_branch,'funcs',funcs_s2,'exclude_trivial',true);
%%
change_uns=find(diff(unstab));
b_x1=arrayfun(@(x)x.parameter(in.PR),po3_branch.point);
%b_x2=arrayfun(@(x)x.parameter(in.df),po3_branch.point);
ymx_pr3=arrayfun(@(x)max(x.profile(1,:)),po3_branch.point);
ymn_pr3=arrayfun(@(x)min(x.profile(1,:)),po3_branch.point);
figure(4)
clf;
plot(b_x1(unstab==0),ymx_pr3(unstab==0),'ob',b_x1(unstab>=1),ymx_pr3(unstab>=1),'ok')
grid on
%%
[pf_funcs,non_symper,suc_v]=SetupPsol(funcs_audi,po3_branch,change_uns(1),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.01; 0,0.05]);

figure(9)
clf;
non_symper=br_contn(pf_funcs,non_symper,100);
non_symper=br_rvers(non_symper);
non_symper=br_contn(pf_funcs,non_symper,100);
%%
unstab_unsy=GetStability(non_symper,'funcs',pf_funcs,'exclude_trivial',true);
bs_x1=arrayfun(@(x)x.parameter(in.PR),non_symper.point);
%b_x2=arrayfun(@(x)x.parameter(in.df),po3_branch.point);
ymxs_pr3=arrayfun(@(x)max(x.profile(1,:)),non_symper.point);
%%
figure(4)
hold on
plot(bs_x1(unstab_unsy==0),ymxs_pr3(unstab_unsy==0),'og',bs_x1(unstab_unsy>=1),ymxs_pr3(unstab_unsy>=1),'ok')
grid on
%%
[~,it2]=min(abs(bs_x1-5.1));
pt2=non_symper.point(it2);
figure(80)
clf;
plot(pt2.mesh*pt2.period,pt2.profile(1:2,:),'LineWidth',2)
grid on
%%
[~,it3]=min(abs(bs_x1-6));
pt3=non_symper.point(it3);
figure(6)
clf;
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',2)
grid on
%
[~,it4]=min(abs(bs_x1-8.7));
pt4=non_symper.point(it4);
figure(8)
clf;
plot(pt4.mesh*pt4.period,pt4.profile(1:2,:),'LineWidth',2)
grid on
%%
po3_branch=br_remove_extracolumns(po3_branch);
non_symper=br_remove_extracolumns(non_symper);
%%
save('identifying_solutions_in_regions_try2.mat')
%%
