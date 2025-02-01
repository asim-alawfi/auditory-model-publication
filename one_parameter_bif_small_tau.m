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
parbd={'min_bound',[in.PR,1;in.df,0],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,08; in.df,0.01; 0,0.02],'print_residual_info',1};
%% Set up symmetry condition for dde-biftool implementation
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
sproj=zeros(6); sproj(1)=1; sproj(8)=1;

% symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)'*[1,1],'condprojmat',[1 0 0 0 0 0],'stateproj',sproj);
% psolsymfine=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
%     'condprojint',linspace(0,1,100)'*[1,1]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
%% picking up a periodic solution from dde23 and continuing computation in PR with fixed df=0.73.
[funcs_tau2,po_symmetry_tinyTau2,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df_tinytau2,in.PR,par_tinytau,...
    'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',{psolsym},'outputfuncs',true);
%%
figure(2)
clf;
po_symmetry_tinyTau2=br_contn(funcs_tau2,po_symmetry_tinyTau2,2000);
po_symmetry_tinyTau2=br_rvers(po_symmetry_tinyTau2);
po_symmetry_tinyTau2=br_contn(funcs_tau2,po_symmetry_tinyTau2,2000);

%% Computing stability and special points
% [po_symmetryTau_wbifs,po2_symmetry_nunst,po2_symmetry_bifs,po2_symmetry_bifind]=MonitorChange(funcs_tau,po_symmetryTau,...
%     'range',2:length(po_symmetryTau.point),'printlevel',1,'print_residual_info',0,...
%     'min_iterations',5);
[po_symmetry_tinyTau2,nunst_sym_tinytau2,dom,triv_defect]=br_stabl(funcs_tau2,po_symmetry_tinyTau2,0,0,'exclude_trivial',1);
bifpoints_tinytau2=find(diff(nunst_sym_tinytau2));
%
rp2_tinytau2=arrayfun(@(x)x.parameter(in.PR),po_symmetry_tinyTau2.point);
max2_tinyytau2=arrayfun(@(x)max(x.profile(1,:)),po_symmetry_tinyTau2.point);
min2_tinyytau=arrayfun(@(x)min(x.profile(1,:)),po_symmetry_tinyTau2.point);
%% plotting stability
figure(2);clf;hold on;grid on
plot(rp2_tinytau2(nunst_sym_tinytau2==0),max2_tinyytau2(nunst_sym_tinytau2==0),'bo',rp2_tinytau2(nunst_sym_tinytau2>=1),max2_tinyytau2(nunst_sym_tinytau2>=1),'kx')%,'LineWidth',1)
plot(rp2_tinytau2(nunst_sym_tinytau2==0),min2_tinyytau(nunst_sym_tinytau2==0),'bo',rp2_tinytau2(nunst_sym_tinytau2>=1),min2_tinyytau(nunst_sym_tinytau2>=1),'kx')%,'LineWidth',1)
%%
po_symmetry_tinyTau2=br_remove_extracolumns(po_symmetry_tinyTau2);


%% Now, branching off from the branch of symmetric solutions to a branch of non-symmetric solutions
% % (a note form me:period doubling bifurcation are expected to be detected).
 S_tinytau2=dde_lincond_struct(size(po_symmetry_tinyTau2.point(1).profile,1),'profile',...
     'shift',[1,2],'condprojmat',[1,0,0,0,0,0],'condprojint',[0,0.5],'stateproj',sproj);
 plm={@(p)p.parameter(po_symmetry_tinyTau2.parameter.free),@(p)dde_psol_lincond(p,S_tinytau2)};
 parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,40; in.df,1],...
     'max_step',[in.PR,0.01; in.df,0.01; 0,0.02],'print_residual_info',1};
 sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
     'condprojint',linspace(0.1,0.5,6)');
 poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
 nspoev1args=addprefix('SetupPOEV1',poev1args);
 [funcs_asy_tinytau2,asym_brtinyTau2,suc_v]=SetupPsol(funcs_audi,po_symmetry_tinyTau2,bifpoints_tinytau2(1),'print_residual_info',1,...
    'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,parbds{:},...
   nspoev1args{:},'plot_measure',plm);
%%
figure(3)
%clf
hold on
asym_brtinyTau2=br_contn(funcs_asy_tinytau2,asym_brtinyTau2,1500);
asym_brtinyTau2=br_rvers(asym_brtinyTau2);
asym_brtinyTau2=br_contn(funcs_asy_tinytau2,asym_brtinyTau2,10);
%%
asym_brtinyTau2=br_remove_extracolumns(asym_brtinyTau2);
%%
[asym_brtinyTau2,nunst_nonsym_tinytau2,domnons,triv_defect0]=br_stabl(funcs_asy_tinytau2,asym_brtinyTau2,0,0,'exclude_trivial',1);

asrp2_tinytau2=arrayfun(@(x)x.parameter(in.PR),asym_brtinyTau2.point);
%%
Sint_A=dde_lincond_struct(size(asym_brtinyTau2.point(1).profile,1),'profile','trafo',Rsym,...
    'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'stateproj',sproj,'condprojint',[0,0.5]);
% Sint_B=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
%     'shift',[1,4],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.25,0.75]);
yax_sym=arrayfun(@(x)dde_psol_lincond(x,Sint_A),po_symmetry_tinyTau2.point);
yax_nonsym=arrayfun(@(x)dde_psol_lincond(x,Sint_A),asym_brtinyTau2.point);
%%
figure(50)
clf; hold on 
grid on
plot(rp2_tinytau2(nunst_sym_tinytau2==0),yax_sym(nunst_sym_tinytau2==0),'b.',...
    rp2_tinytau2(nunst_sym_tinytau2>=1),yax_sym(nunst_sym_tinytau2>=1),'k.','LineWidth',1)
plot(asrp2_tinytau2(nunst_nonsym_tinytau2==0),yax_nonsym(nunst_nonsym_tinytau2==0),'g.','MarkerSize',10)
plot(asrp2_tinytau2(nunst_nonsym_tinytau2>=1),yax_nonsym(nunst_nonsym_tinytau2>=1),'k.','MarkerSize',30)
%%
save('one_parameter_bif_small_tau.mat')
