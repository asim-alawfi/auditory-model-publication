clear;
format compact
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
    load('branch_of_sympos_original_case_and_threshold_crossing_coj.mat')
%% (1) Pick PO at PR=38 and then we continue in df
rpl_x=arrayfun(@(x)x.parameter(in.PR),po2_symmetry_wbifs.point);
[~,ind]=min(abs(rpl_x-38));
%% (2) continuoe in df
%%% 
parbd={'min_bound',[in.PR,1;in.df,0.02],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.002; 0,0.002],'print_residual_info',1};
%%%
[branch_prlarge,suc_r]=ChangeBranchParameters(funcs_s2,po2_symmetry_wbifs,ind,'degree',6,'intervals',120,'contpar',in.df,'step',1e-3,'correc',true,parbd{:});%,...
      %'extra_condition',true,'phase_condition',0,'eigmatrix', 'sparse','matrix', 'sparse');
%
branch_prlarge=br_rvers(branch_prlarge);
figure(4)
clf;
hold on
branch_prlarge=br_contn(funcs_s2,branch_prlarge,700);
%%
for i=300: length(branch_prlarge.point)
    pp=branch_prlarge.point(i);
    figure(17)
    clf
    hold on
    plot(pp.mesh,pp.profile(1:2,:),'Linewidth',2)
    yline(0.5,'k--','Linewidth',2)
    title('pr=',i)
    pause(0.1)
end
%%
branch_prlarge=br_remove_extracolumns(branch_prlarge)
save('results_PR_large.mat')
%%
pp=branch_prlarge.point(456);
figure(1)
clf
hold on
plot(pp.mesh,pp.profile(1:2,:),'Linewidth',2)
yline(0.5,'k--','Linewidth',2)
%title('pr=',i)
pause(0.1)
[it,vv]=min(pp.profile(1,:))
%%
%% trace symmetric periodic orbit with maximum equal to the threshold (theta=0.5) in two parameters
% find extrema of $x$  along orbits on po_symmetry branch and pick orbits that
% have two extrema
smaxval=par(in.theta); % 
c_A=[1,0,0,0,0,0]; % c for u_A: to only extract the values of u_A
sympo_uA_extrema=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',branch_prlarge.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(branch_prlarge.point),sympo_uA_extrema);
[~,theta_cross]=min(abs(sympomax_ua+(smaxval+1e-1)));
% pp=branch_prlarge.point(theta_cross)


%%
second_max=1;
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
sympo0=setfield(branch_prlarge,'point',branch_prlarge.point(456));
sympo0.point.parameter([ine.t0,ine.val])=[pp.mesh(vv),smaxval];%[sympo_uA_extrema{theta_cross}(second_max),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
[Sfuncs,br_stat_touch,suc_max]=ChangeBranchParameters(funcs_s2,sympo0,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{psolsym,max_cond},'outputfuncs',true,...
    'print_residual_info',1);
%%
br_stat_touch.parameter.max_step(:,2)=[0.05;0.05;0.05];
figure(3)
clf;
br_stat_touch=br_contn(Sfuncs,br_stat_touch,150);
%%
br_stat_touch=br_rvers(br_stat_touch);
%%
br_stat_touch.parameter.max_step(:,2)=[0.05;0.01;0.01];
br_stat_touch.parameter.min_bound(end)=0.001;
%
br_stat_touch=br_contn(Sfuncs,br_stat_touch,200);

rpp_x=arrayfun(@(x)x.parameter(in.PR),br_stat_touch.point);

rdff_x=arrayfun(@(x)x.parameter(in.df),br_stat_touch.point);
figure(6)
hold on
plot(rpp_x,rdff_x,'k-','LineWidth',3)
%%
br_stat_touch=br_remove_extracolumns(br_stat_touch)
save('PR_large_df_small_part2.mat')
%%
for i=1:length(br_stat_touch.point)
    pp=br_stat_touch.point(i);
    figure(1)
    clf
    hold on
    plot(pp.mesh,pp.profile(1:2,:),'Linewidth',2)
    yline(0.5,'k--','Linewidth',2)
    ylim([0,1])
    title('pr=',i)
    pause(0.1)
end
%%
crossing_tha=diff(sign(pp.profile(1,:)-0.5))==2;
crossing_thb=diff(sign(pp.profile(2,:)-0.5))==2;
crossing_inda=find(crossing_tha~= 0);
crossing_indb=find(crossing_thb~= 0);
pp=branch_prlarge.point(420);
figure(18)
clf
hold on
plot(pp.mesh,pp.profile(1:2,:),'Linewidth',2)
plot(pp.mesh(crossing_inda),pp.profile(1,crossing_inda),'k.','MarkerSize',30)
plot(pp.mesh(crossing_indb),pp.profile(2,crossing_indb),'k.','MarkerSize',30)
yline(0.5,'k--','Linewidth',2)
%title('pr=',i)
ylim([0,1])
figure(6) 
hold on
plot(pp.parameter(in.PR),pp.parameter(in.df),'x','MarkerSize',6,'LineWidth',2)
%%
% load('PR_large_df_small_part2.mat')
load('branching_towrds_asymmetric_sols_original_case.mat')
load('sym_breaking_original_case_try2.mat')
load('identifying_solutions_in_regions_try2.mat')
%
set_one_bif={'LineWidth',2,'Box','on','FontSize',14,'FontWeight','normal'};

clrs2=lines();
graycolor=0.7*[1 1 1];
[~,it]=min(abs(rp_thta-34));
pt=mbr_wbifs.point(it);
t5=find(diff(sign(pr_m-6.658)));
pt5=mbranch_df.point(t5(1));
figure(1140)
clf;
hold on; 
%%%%%%%%%%%%%%%%%% Bifurcation results
plt1=plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b-','color',clrs2(1,:),'LineWidth',4);
%plot(rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plt_thsy=plot(rp_thta(dum==0),df_thta(dum==0),'m-','Linewidth',3);
plt_thsy2=plot(rp_thta(dum>=1),df_thta(dum>=1),'','color',graycolor,'Linewidth',2);
plt_thunsym=plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'--','Color','k','LineWidth',3);
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',20,'LineWidth',3)
plot(pr_m(t5(1)),df_m(t5(1)),'x','Color',clrs2(3,:),'MarkerSize',20,'LineWidth',3)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'rx',...
    'Color',clrs2(7,:),'MarkerSize',20,'LineWidth',3)
plot(rpp_x,rdff_x,'k-','LineWidth',3)

%plot(non_symper.point(it4).parameter(in.PR),non_symper.point(it4).parameter(in.df),'gx','MarkerSize',25,'LineWidth',3);
%%%%%%%% Experimental range
ft=17;
plt_text={'sym-breaking','touching \theta: stable sym','touching \theta: unstable sym ','touching \theta: stable non-sym'};
plt_vect=[plt1;plt_thsy;plt_thsy2;plt_thunsym];
legend(plt_vect,plt_text,'FontSize',ft,'FontWeight','normal')
%%%%%%%%%%%%

set(gca, set_one_bif{:})
text(5.5,0.1 ,'1 (int)', 'FontSize', 20, 'Color', 0.4*[1 1 1],'FontWeight','bold')
text(7,0.25,'2 (bist)', 'FontSize', 20, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
text(15,0.5,'3 (seg)', 'FontSize', 20, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=20,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=20,FontWeight='bold')
title('A','FontSize',ft,'FontWeight','normal')
ylim([0,1] )
yticks([0,0.5,1])
%%
[br_stat_touch,nunst_pr,dom_pr,triv_defect_pr]=br_stabl(Sfuncs,br_stat_touch,0,0,'exclude_trivial',1);