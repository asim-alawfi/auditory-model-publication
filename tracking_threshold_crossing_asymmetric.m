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
 load('branch_of_sympos_original_case_and_threshold_crossing.mat')
%% One-parameter continuation with fixing PR=25, and varying df. 
% We n
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.01; in.df,0.01; 0,0.005],'print_residual_info',1};
[~,it]=min(abs(rp2_x-25));
[br_df,sucp]=ChangeBranchParameters(funcs_s2,po2_symmetry_wbifs,it,...
      'degree',6,'intervals',120,'contpar',in.df,'step',1e-3,'correc',true,parbds{:},...
      'extra_condition',true,'phase_condition',0);
br_df=br_rvers(br_df);
figure(4)
clf;
hold on
br_df=br_contn(funcs_s2,br_df,300);
br_df=br_rvers(br_df);
br_df=br_contn(funcs_s2,br_df,10000);
[nunst_df,dom_df,triv2_defect_df,br_df.point]=GetStability(br_df,'funcs',funcs_s2,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df1=find(diff(nunst_df));
%% plotting one-par bif for symmetric POs in (df,max(u_A))-plan with stability
df_x1=arrayfun(@(x)x.parameter(in.df),br_df.point);
df_pr_x1=arrayfun(@(x)x.parameter(in.PR),br_df.point);
ymx_df=arrayfun(@(x)max(x.profile(1,:)),br_df.point);
ymn_df=arrayfun(@(x)min(x.profile(1,:)),br_df.point);
figure(49)
clf;
plot(df_x1(nunst_df==0),ymx_df(nunst_df==0),'ob',...
    df_x1(nunst_df>=1),ymx_df(nunst_df>=1),'xk')
xlabel('df')
ylabel('max u_A')
grid on
%%  Branching-off twoards asymmetric solutions ( varying df)
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_df,nonsymper_df,suc_v]=SetupPsol(funcs_audi,br_df,chang_stb_df1(2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.df,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.005; 0,0.001]);
figure(99)
clf;
hold on
nonsymper_df=br_contn(funcs_df,nonsymper_df,447);
nonsymper_df=br_remove_extracolumns(nonsymper_df);
for i=1:1%length(nonsymper_df.point)
xp=nonsymper_df.point(141);
figure(100)
clf;
plot(xp.mesh*xp.period,xp.profile(1:2,:),'LineWidth',2)
grid on
drawnow
end
% Plloting stability for asymmetric POs in (df,max(u_A))-plan
[nunst_dfs,dom_dfs,triv2_defect_dfs,nonsymper_df.point]=GetStability(nonsymper_df,'funcs',funcs_df,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df=find(diff(nunst_dfs));
df_xs=arrayfun(@(x)x.parameter(in.df),nonsymper_df.point);
df_xs_pr=arrayfun(@(x)x.parameter(in.PR),nonsymper_df.point);
ymxs_df=arrayfun(@(x)max(x.profile(1,:)),nonsymper_df.point);
ymns_df=arrayfun(@(x)min(x.profile(1,:)),nonsymper_df.point);
%%
figure(40)
clf;
hold on
plot(df_xs(nunst_dfs==0),ymxs_df(nunst_dfs==0),'og',...
    df_xs(nunst_dfs>=1),ymxs_df(nunst_dfs>=1),'xk')
plot(df_x1(nunst_df==0),ymx_df(nunst_df==0),'ob',...
    df_x1(nunst_df>=1),ymx_df(nunst_df>=1),'xk')
xlabel('d_f')
ylabel('max u_A')
title('varying d_f at fixed PR=25')
grid on
%
figure(45)
clf;
hold on
plot(df_xs(nunst_dfs==0),ymxs_df(nunst_dfs==0),'og',...
    df_xs(nunst_dfs>=1),ymxs_df(nunst_dfs>=1),'xk')
%
figure(45)
hold on
plot(df_xs(141),ymxs_df(141),'sr','LineWidth',3)
plot(df_xs(425),ymxs_df(425),'sg','LineWidth',3)
%
figure(3)
clf;
hold on
plot(df_xs_pr(141),df_xs(141),'sr','LineWidth',3)
%plot(df_xs_pr(425),df_xs(425),'sg','LineWidth',3)
%%
br_df=br_remove_extracolumns(br_df);
nonsymper_df=br_remove_extracolumns(nonsymper_df);
%save('tracking_threshold_crossing_unsymmetric_part1.mat')

%%
c_A=[1,0,0,0,0,0];
smaxval=0.5;
uA_extrema=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',nonsymper_df.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(nonsymper_df.point),uA_extrema);
%%
figure(999)
clf;
plot(sympomax_ua,'o')
grid on
%% We pick the PO solution with touching the threshold and continiue in (pr,df,t), then we plot in (pr,df)-space
[~,it_cross]=min(abs(sympomax_ua-(smaxval+1e-4)));
second_ua_peak=3;
%icross4=141;
%
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
% symmetric conditios
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
sympo=setfield(nonsymper_df,'point',nonsymper_df.point(it_cross));
sympo.point.parameter([ine.t0,ine.val])=[uA_extrema{it_cross}(second_ua_peak),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
[mfuncs_df,mbranch_df,suc_max]=ChangeBranchParameters(funcs_df,sympo,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{max_cond},'outputfuncs',true,...
    'print_residual_info',1,'max_step',[in.PR,0.1; in.df,0.05; 0,0.1]);
mbranch_df.parameter.max_bound(3)=40;
%%
figure(12)
clf;
hold on
mbranch_df=br_contn(mfuncs_df,mbranch_df,400);
mbranch_df=br_rvers(mbranch_df);
mbranch_df=br_contn(mfuncs_df,mbranch_df,700);
%%
[mbranch_df_wbifs,df_tests,m5_bifs,df_bifind]=MonitorChange(mfuncs_df,mbranch_df,'print_residual_info',0);
nunst_m=GetStability(mbranch_df_wbifs,'funcs',mfuncs_df,...
    'exclude_trivial',true);
df_m=arrayfun(@(x)x.parameter(in.df),mbranch_df_wbifs.point);
pr_m=arrayfun(@(x)x.parameter(in.PR),mbranch_df_wbifs.point);
%%
figure(11)
clf;
hold on %clf;
plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
    pr_m(nunst_m>=1),df_m(nunst_m>=1),'yo','LineWidth',1)
grid on
%%
mbranch_df_wbifs=br_remove_extracolumns(mbranch_df_wbifs);
mbranch_df=br_remove_extracolumns(mbranch_df);
%%
save('tracking_threshold_crossing_asymmetric_try2.mat')
%%  simulation of asymmetric POs with touching threshold
g=@(x)0.5;
for i=20:20%200%length(mbranch_df_wbifs.point)
xp=mbranch_df_wbifs.point(i);
perTR=xp.period;
figure(1100)
clf;
hold on
plot(xp.mesh*xp.period,xp.profile(1:2,:),'LineWidth',2)
fplot(g,[0,perTR],'g--','LineWidth',2);
grid on
axis tight %drawnow
legend('u_A','u_B')
end
