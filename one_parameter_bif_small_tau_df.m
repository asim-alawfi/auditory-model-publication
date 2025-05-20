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
%%  
% upload from file 'Solve23_2df'
load('dde23_small_tau.mat')
%% SetupFuncs and parameter bounds.
funcs_audi=set_symfuncs(@symbolic_auditory_with_symmetry_version,'sys_tau',@()[in.D, in.TD],...
   'sys_cond',@sys_cond);
parbd={'min_bound',[in.PR,1;in.df,0.02],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,08; in.df,0.01; 0,0.01],'print_residual_info',1};
%% Set up symmetry condition for dde-biftool implementation
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
sproj=zeros(6); sproj(1:2,1:2)=eye(2);% make all zeros except u_A in the first column and u_B in the second column.

% symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.25,3)'*[1,1],'condprojmat',[1 0 0 0 0 0],'stateproj',sproj);
% psolsymfine=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
%     'condprojint',linspace(0,1,100)'*[1,1]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
%% picking up a periodic solution from dde23 and continuing computation in PR with fixed df=0.73.
[funcs_tau2_df,po_symmetry_tinyTau2_df,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df_tinytau2,in.df,par_tinytau,...
    'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',6,'intervals',120,'extracolumns','auto','newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',{psolsym},'outputfuncs',true);
%%
figure(2)
clf;
po_symmetry_tinyTau2_df=br_contn(funcs_tau2_df,po_symmetry_tinyTau2_df,2000);
po_symmetry_tinyTau2_df=br_rvers(po_symmetry_tinyTau2_df);
po_symmetry_tinyTau2_df=br_contn(funcs_tau2_df,po_symmetry_tinyTau2_df,2000);

%% Computing stability 
[po_symmetry_tinyTau2_df,nunst_sym_tinytau2_df,dom_df,triv_defect_df]=br_stabl(funcs_tau2_df,po_symmetry_tinyTau2_df,0,0,'exclude_trivial',1);
bifpoints_tinytau2_df=find(diff(nunst_sym_tinytau2_df));
%
df2_tinytau2=arrayfun(@(x)x.parameter(in.df),po_symmetry_tinyTau2_df.point);
max2_tinyytau2_df=arrayfun(@(x)max(x.profile(1,:)),po_symmetry_tinyTau2_df.point);
min2_tinyytau_df=arrayfun(@(x)min(x.profile(1,:)),po_symmetry_tinyTau2_df.point);
%% plotting stability
figure(2);clf;hold on;grid on
plot(df2_tinytau2(nunst_sym_tinytau2_df==0),max2_tinyytau2_df(nunst_sym_tinytau2_df==0),'bo',df2_tinytau2(nunst_sym_tinytau2_df>=1),max2_tinyytau2_df(nunst_sym_tinytau2_df>=1),'kx')%,'LineWidth',1)
%plot(df2_tinytau2(nunst_sym_tinytau2_df==0),min2_tinyytau_df(nunst_sym_tinytau2_df==0),'bo',df2_tinytau2(nunst_sym_tinytau2_df>=1),min2_tinyytau_df(nunst_sym_tinytau2_df>=1),'kx')%,'LineWidth',1)
%%
po_symmetry_tinyTau2_df=br_remove_extracolumns(po_symmetry_tinyTau2_df);


%% Now, branching off from the branch of symmetric solutions to a branch of non-symmetric solutions
% % (a note form me:period doubling bifurcation are expected to be detected).
 S_tinytau2=dde_lincond_struct(size(po_symmetry_tinyTau2_df.point(1).profile,1),'profile',...
     'shift',[1,2],'condprojmat',[1,0,0,0,0,0],'condprojint',[0,0.5],'stateproj',sproj);
 plm={@(p)p.parameter(po_symmetry_tinyTau2_df.parameter.free),@(p)dde_psol_lincond(p,S_tinytau2)};
 parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,40; in.df,1],...
     'max_step',[in.PR,0.01; in.df,0.01; 0,0.01],'print_residual_info',1};
 sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
     'condprojint',linspace(0.1,0.5,6)');
 poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
 nspoev1args=addprefix('SetupPOEV1',poev1args);
 [funcs_asy_tinytau2_df,asym_brtinyTau2_df,suc_v]=SetupPsol(funcs_audi,po_symmetry_tinyTau2_df,bifpoints_tinytau2_df(2),'print_residual_info',1,...
    'outputfuncs',true,'branch_off','POEV1','contpar',in.df,parbds{:},...
   nspoev1args{:},'plot_measure',plm);
%%
figure(3)
%clf
hold on
asym_brtinyTau2_df=br_contn(funcs_asy_tinytau2_df,asym_brtinyTau2_df,165);

%%
asym_brtinyTau2_df=br_remove_extracolumns(asym_brtinyTau2_df);
[asym_brtinyTau2_df,nunst_nonsym_tinytau2_df,domnons_df,triv_defect_df]=br_stabl(funcs_asy_tinytau2_df,asym_brtinyTau2_df,0,0,'exclude_trivial',1);
asdf2_tinytau2=arrayfun(@(x)x.parameter(in.df),asym_brtinyTau2_df.point);
Sint_A=dde_lincond_struct(size(asym_brtinyTau2_df.point(1).profile,1),'profile','trafo',Rsym,...
    'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'stateproj',sproj,'condprojint',[0,0.5]);
% Sint_B=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
%     'shift',[1,4],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.25,0.75]);
yax_sym_df=arrayfun(@(x)dde_psol_lincond(x,Sint_A),po_symmetry_tinyTau2_df.point);
yax_nonsym_df=arrayfun(@(x)dde_psol_lincond(x,Sint_A),asym_brtinyTau2_df.point);
%
figure(50)
clf; hold on 
grid on
plot(df2_tinytau2(nunst_sym_tinytau2_df==0),yax_sym_df(nunst_sym_tinytau2_df==0),'b.',...
    df2_tinytau2(nunst_sym_tinytau2_df>=1),yax_sym_df(nunst_sym_tinytau2_df>=1),'k.','LineWidth',1)
plot(asdf2_tinytau2(nunst_nonsym_tinytau2_df==0),yax_nonsym_df(nunst_nonsym_tinytau2_df==0),'g.','MarkerSize',10)
plot(asdf2_tinytau2(nunst_nonsym_tinytau2_df>=1),yax_nonsym_df(nunst_nonsym_tinytau2_df>=1),'k.','MarkerSize',30)
%%
save('one_parameter_bif_small_tau_df.mat')
