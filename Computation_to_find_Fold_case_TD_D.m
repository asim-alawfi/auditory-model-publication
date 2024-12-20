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
load('br_crossing_threshold_try2.mat')
%% We compute a branch of one-parameter in df, then we fixed df=0.47,
% and then we change the compuation in PR
parbds={'min_bound',[in.PR,2;in.df,0.35],'max_bound',[in.PR,5; in.df,1],...
    'max_step',[in.PR,0.01; in.df,0.01; 0,0.005],'print_residual_info',1};
[~,it_pr]=min(abs(rp2_x-5));
[br_in_df_fold,sucp]=ChangeBranchParameters(funcs_ss,br_symmetry_wbifs(2),it_pr,...
    'degree',6,'intervals',120,'contpar',in.df,'step',0.01,'correc',true,parbds{:},...
    'extra_condition',true,'phase_condition',0,...
    'eigmatrix', 'sparse','matrix', 'sparse');
%%
figure(1)
clf
hold on
br_in_df_fold=br_rvers(br_in_df_fold);
br_in_df_fold=br_contn(funcs_ss,br_in_df_fold,300);
% br_in_df_fold=br_rvers(br_in_df_fold);
% br_in_df_fold=br_contn(funcs_ss,br_in_df_fold,200);
%% Now, 3from the previous branch we fixed df=0.47 and then we continue in PR
df2_x2=arrayfun(@(x)x.parameter(in.df),br_in_df_fold.point);
[~,it_df]=min(abs(df2_x2-0.47));
[br_pr,sucp]=ChangeBranchParameters(funcs_ss,br_in_df_fold,it_df,...
    'degree',6,'intervals',120,'contpar',in.PR,'step',0.001,'correc',true,parbds{:},...
    'extra_condition',true,'phase_condition',0,...
    'eigmatrix', 'sparse','matrix', 'sparse');
%%
figure(2)
clf
hold on
br_pr=br_contn(funcs_ss,br_pr,200);
br_pr=br_rvers(br_pr);
br_pr=br_contn(funcs_ss,br_pr,600);
%%
for i=1:length(br_pr.point)
    pm=br_pr.point(i);
    figure(9)
    clf
    hold on
    plot(pm.mesh,pm.profile(1:2,:))
    yline(0.5,'k--')
    grid on
end
%%
[br_pr,unst_br_pr,~,~]=br_stabl(funcs_ss,br_pr,0,0,'exclude_trivial',true);


%%
parbds={'min_bound',[in.PR,1;in.df,0.35],'max_bound',[in.PR,5; in.df,1],...
    'max_step',[in.PR,0.01; in.df,0.01; 0,0.01],'print_residual_info',1};
bif_loc=find(diff(unst_br_pr));
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_asy2,br_pr_asym,suc_v]=SetupPsol(funcs_audi,br_pr, bif_loc(end),'print_residual_info',1,...
    'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,parbds{:},...
    nspoev1args{:});
%%
figure(334)
hold on
br_pr_asym=br_contn(funcs_asy2,br_pr_asym,10);
%%
[br_pr_asym,unst_non_sym,~,~]=br_stabl(funcs_asy2,br_pr_asym,0,0,'exclude_trivial',true);
pr2_x2=arrayfun(@(x)x.parameter(in.PR),br_pr.point);
pr2_nonsym=arrayfun(@(x)x.parameter(in.PR),br_pr_asym.point);

Sint_A1=dde_lincond_struct(size(br_pr.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[1,0,0,0,0,0],'condprojint',[0,0.5]);
Sint_B2=dde_lincond_struct(size(br_pr.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[0,1,0,0,0,0],'condprojint',[0.5,1]);
yax_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),br_pr.point);
yax_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),br_pr.point);
non_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),br_pr_asym.point);
non_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),br_pr_asym.point);
yax_sym=yax_Sint_A1-yax_Sint_B2;
yax_nonsym=non_sym_Sint_A1-non_sym_Sint_B2;
%%
set_one_bif={'FontWeight','bold','FontSize',12,'FontName','Aril'};
loc_bif_non_sym=find(diff(unst_non_sym));
clrs2=lines();
figure(1)
clf;
tiledlayout(1,2)
nexttile
hold on; grid on
plot(pr2_x2(unst_br_pr==0),yax_sym(unst_br_pr==0),'.','Color',clrs2(4,:),...
    'LineWidth',1,'MarkerSize',10)
plot(pr2_x2(unst_br_pr>=1),yax_sym(unst_br_pr>=1),'k--','LineWidth',3)%,'Color',grayColor)%,'LineWidth',1)
plot(pr2_nonsym(unst_non_sym==0),yax_nonsym(unst_non_sym==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10)
plot(pr2_nonsym(unst_non_sym>=1),yax_nonsym(unst_non_sym>=1),'.','Color',clrs2(3,:),'LineWidth',3,'MarkerSize',10)%,'g.')%,'MarkerS',5)
plot(pr2_x2(bif_loc),yax_sym(bif_loc),'k.','MarkerSize',35)
plot(pr2_nonsym(loc_bif_non_sym(1)),yax_nonsym(loc_bif_non_sym(1)),'.','Color',clrs2(2,:),'MarkerSize',35)
plot(pr2_nonsym(loc_bif_non_sym(end)),yax_nonsym(loc_bif_non_sym(end)),'.','Color',clrs2(2,:),'MarkerSize',35)
xlim([3.5,5])
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
ylabel('$\int_0^{1/2}u_B(t)-u_A(t+1/2)\mathrm{d}t$','interpreter','latex','FontSize',16);
legend('sym stable POs','sym untsable POs','non-sym stable POs','non-sym untsable POs',...
    'symmetry-breaking','','fold of POs','FontSize',12)
title(' $(t_\mathrm{d }, D) = (0.05,0.05)$ with $d_\mathrm{f}=0.47$','interpreter','latex','FontSize',16)
set(gca,set_one_bif{:})
nexttile
td=[0.022,0.05];
D=[0.015,0.05];
plot(td(1),D(1),'bo',td(1),D(2),'ro',td(2),D(1),'ko',td(2),D(2),'go','LineWidth',3)
xlim([0.0 0.07])
ylim([0.0 0.07])
xlabel('$t_\mathrm{d}$','FontName','Courier','Interpreter','latex','FontSize',12)
ylabel('$D$','interpreter','latex','FontSize',16)
xticks([0,0.022,0.05])
yticks([0,0.022,0.05])
set(gca,set_one_bif{:})
grid on
%%
br_in_df_fold=br_remove_extracolumns(br_in_df_fold);
br_pr=br_remove_extracolumns(br_pr);
br_pr_asym=br_remove_extracolumns(br_pr_asym);
%%
save('one_paramter_bif_fold_case.mat')
