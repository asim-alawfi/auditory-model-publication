clear
   base=[pwd(),filesep(),'ddebiftool_snapshot_23October2022',filesep()];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco'],...
            [base,'Supporting_functions'])
format compact
%load('one_parameter_bif_small_tau.mat')
load('one_parameter_bif_small_tau_df.mat')
%%
parbd={'min_bound',[in.PR,1;in.df,0],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.08; 0,0.08],'print_residual_info',1};
pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.4,2)','stateproj',sproj,'condprojmat',[1,0,0,0,0,0]);
%pfvsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
  %  'condprojint',linspace(0.1,0.5,6)');
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
[funcs_symbrk_tu,symbk_brtu_ini,succ]=SetupPOEV1(funcs_tau2_df,po_symmetry_tinyTau2_df,bifpoints_tinytau2_df(1),...
    'contpar',[in.PR,in.df],'dir',in.df,parbd{:},...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
%%
symbk_brtu=symbk_brtu_ini;
figure(8)
clf;
hold on
symbk_brtu=br_contn(funcs_symbrk_tu,symbk_brtu,500);
symbk_brtu=br_rvers(symbk_brtu);
figure(8)
hold on
symbk_brtu=br_contn(funcs_symbrk_tu,symbk_brtu,950);
%%
symbk_brtu=br_remove_extracolumns(symbk_brtu);
save('sym_breaking_small_tau_upperBranch.mat')

%% The second symmetry-breaking branch (lower branch)
[funcs_symbrk_tu,symbk_brtu_ini2,succ]=SetupPOEV1(funcs_tau2_df,po_symmetry_tinyTau2_df,bifpoints_tinytau2_df(2),...
    'contpar',[in.PR,in.df],'dir',in.df,parbd{:},...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:}); 
symbk_brtu2=symbk_brtu_ini2;
%%
symbk_brtu2.parameter.max_step(:,2)=[0.1;0.1;0.1];
figure(8)
hold on
symbk_brtu2=br_contn(funcs_symbrk_tu,symbk_brtu2,400);
%%
figure(8)
hold on
symbk_brtu2=br_rvers(symbk_brtu2);
symbk_brtu2=br_contn(funcs_symbrk_tu,symbk_brtu2,400);
%%
symbk_brtu2=br_remove_extracolumns(symbk_brtu2);
symbk_brtu=br_remove_extracolumns(symbk_brtu);
% for i=1:length(symbk_brtu.point)

rp_syktinytau2=arrayfun(@(x)x.parameter(in.PR),symbk_brtu2.point);
df_syktinytau2=arrayfun(@(x)x.parameter(in.df),symbk_brtu2.point);
%

figure(70)
clf; 
hold on
plot(rp_syktinytau2,df_syktinytau2,'r.','MarkerSize',5)
grid on
title('tau=0.0.0025')
 %
 save('sym_breaking_small_tau_lowerbranch.mat')
 %%
% %% trace symmetric periodic orbit with maximum equal to the threshold (theta=0.5) in two parameters
% find extrema of $x$  along orbits on po_symmetry branch and pick orbits that
% have two extrema
c_A=[1,0,0,0,0,0];
smaxval=0.5;
sympo_uA_extrema=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',asym_brtinyTau2_df.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
% for i=1:length(po_symmetry_tinyTau2.point)
sympomax_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(asym_brtinyTau2_df.point),sympo_uA_extrema);
[~,theta_cross]=min(abs(sympomax_ua-(smaxval)));
%
%theta_cross=theta_cross
second_max=4;
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
sympo0=setfield(asym_brtinyTau2_df,'point',asym_brtinyTau2_df.point(theta_cross));
sympo0.point.parameter([ine.t0,ine.val])=[sympo_uA_extrema{theta_cross}(second_max),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
%
pt=asym_brtinyTau2_df.point(theta_cross);
figure(77)
clf
hold on; grid on
plot(pt.mesh,pt.profile(1:2,:),'LineWidth',2)
%%
parbd={'min_bound',[in.PR,1;in.df,0.003],'max_bound',[in.PR,75; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.08; 0,0.08],'print_residual_info',1};
[mfuncs,mbranch_ini,suc_max]=ChangeBranchParameters(funcs_audi,sympo0,1,'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{max_cond},'outputfuncs',true,'print_residual_info',1,parbd{:});

%%
mbranch=mbranch_ini;
figure(8)
hold on
mbranch=br_contn(mfuncs,mbranch,500);


%%
% mbranch.parameter.max_step(4)=0.02;
% mbranch.parameter.max_step(5)=0.02;
% mbranch.parameter.max_step(end)=0.02;
figure(8)
mbranch=br_rvers(mbranch);
mbranch=br_contn(mfuncs,mbranch,1000);
%%
% clear mbranch1
% mbranch=br_reorder(mbranch,in.PR);
% rp_ss_tau=arrayfun(@(x)x.parameter(in.PR),branch.point);
% [~,it]=min(abs(rp_ss_tau-61.5092));
% mbranch1=br_remove_points(mbranch,in.PR,it);
% %%
% mbranch1.parameter.max_step(4)=0.05;
% mbranch1.parameter.max_step(5)=0.05;
% mbranch.parameter.max_step(end)=0.05;
% %%
% figure(70)
% mbranch1=br_contn(mfuncs,mbranch1,1000);
%%
rp_ss_tau=arrayfun(@(x)x.parameter(in.PR),mbranch.point);
df_ss_tau=arrayfun(@(x)x.parameter(in.df),mbranch.point);
%%
figure(708)
%clf; 
hold on
plot(rp_ss_tau,df_ss_tau,'k.','MarkerSize',5)
%%

mbranch=br_remove_extracolumns(mbranch);
mbranch1=br_remove_extracolumns(mbranch);
save('sym_breaking_small_tau_coh_and_fiss.mat')
