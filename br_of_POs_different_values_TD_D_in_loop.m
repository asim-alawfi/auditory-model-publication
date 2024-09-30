 
% Inevestigating the impact of the time duration and delay parameters for
% the auditory model. 
% we vary PR where we consider different values for TD and D: 
% Case (1) we increase the time duration to  TD=0.05 (which was TD=0.022)
% Case (2) we increase the time duration and the delay D to  TD=D=0.05 (which was TD=0.022, D=0.015)
% Case (3) we increase the delay to D=0.05 (which was D=0.015)
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
%% We pick initial POs for the continuation for each case using dde23
load('sol23_j.mat') 
load('sol23_j3.mat') 
%%  Set up initial date
parbd={'min_bound',[in.PR,1;in.df,0],'max_bound',[in.PR,20.2; in.df,1],...
'max_step',[in.PR,0.1; in.df,0.01; 0,0.1],'print_residual_info',1};
funcs_audi=set_symfuncs(@symbolic_auditory_with_symmetry_version,'sys_tau',@()[in.D, in.TD],...
'sys_cond',@sys_cond);
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
% symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
'condprojint',linspace(0.1,0.5,6)'*[1,1]);
psolsymfine=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
'condprojint',linspace(0,1,100)'*[1,1]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
'UniformOutput',false),args(2:2:end)),1,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  (1) branches of one-parameter bifurcation by varying PR with different values of TD and D.
dde23_sols={sol23_j1,sol23_j2,sol23_j4};
dde23_pars={par_j,par_j1,par_j3};
%%   One-parameter bifurcation of POs by varying PR
% br_symmetry(1): with TD=0.05, D=0.15 (case: increasing TD)
% br_symmetry(2): with TD=D=0.05       (case: TD=D) 
% br_symmetry(3): with TD=0.022, D=0.05 (case: increasing D)
for i=1:3
[funcs_ss,br_symmetry(i),suc(i)]=branch_from_sol_audi_sym(funcs_audi,dde23_sols{i},in.PR,dde23_pars{i},...
'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
'adapt_mesh_after_correct',1,...funcs_audi
'usercond',{psolsym},'outputfuncs',true);
figure(i)
clf;
br_symmetry(i)=br_contn(funcs_ss,br_symmetry(i),1000);
br_symmetry(i)=br_rvers(br_symmetry(i));
br_symmetry(i)=br_contn(funcs_ss,br_symmetry(i),1000);
br_symmetry(i)=br_stabl(funcs_ss,br_symmetry(i),0,0);
[br_symmetry_wbifs(i),nunst_po{i},bif(i,:),p_bif(:,i)]=MonitorChange(funcs_ss,br_symmetry(i),...
'range',2:length(br_symmetry(i).point),'printlevel',1,'print_residual_info',0,...
'min_iterations',5);
nunst_g{i}=GetStability(br_symmetry_wbifs(i),'exclude_trivial',true);
end
%% remove extra columns
for i=1:3
br_symmetry_wbifs(i)=br_remove_extracolumns(br_symmetry_wbifs(i));
br_symmetry(i)=br_remove_extracolumns(br_symmetry(i));
end
%%  Plotting stability
rp1_x=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(1).point);
max1_y=arrayfun(@(x)max(x.profile(1,:)),br_symmetry_wbifs(1).point);
max1_y2=arrayfun(@(x)min(x.profile(1,:)),br_symmetry_wbifs(1).point);
%
rp2_x=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(2).point);
max2_y=arrayfun(@(x)max(x.profile(1,:)),br_symmetry_wbifs(2).point);
max2_y2=arrayfun(@(x)min(x.profile(1,:)),br_symmetry_wbifs(2).point);
%
rp3_x=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(3).point);
max3_y=arrayfun(@(x)max(x.profile(1,:)),br_symmetry_wbifs(3).point);
max3_y2=arrayfun(@(x)min(x.profile(1,:)),br_symmetry_wbifs(3).point);
%% save the computation of one-parameter bifurcation in PR
%save('one-parameter_bif_for_all_cases_of_td_d.mat')
%%  Tracking symmetry branch with touching threshold  
c_A=[1,0,0,0,0,0];
smaxval=par_j(in.theta); % the value of the threshold (our maximum tracking value)
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
for i=1:3
sym_uA_extema=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',br_symmetry(i).point,'uniformoutput',false);
second_max_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(br_symmetry(i).point),sym_uA_extema); % we pick the second maximum of u_A
[~,theta_cross]=min(abs(second_max_ua-(smaxval+1e-4)));
itcross3=3; % It is supposed that we have 4 extrema, we pick the third one (i.e the second max of u_A) 
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
symmetry_po(i)=setfield(br_symmetry(i),'point',br_symmetry(i).point(theta_cross));
symmetry_po(i).point.parameter([ine.t0,ine.val])=[sym_uA_extema{theta_cross}(itcross3),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
[mfuncs_ss,br_crossing(i),suc_max(i)]=ChangeBranchParameters(funcs_ss,symmetry_po(i),1,...
'contpar',[ine.PR,ine.df,ine.t0],...
'usercond',{psolsym,max_cond},'outputfuncs',true,...
'print_residual_info',1);
end
br_crossing(1).parameter.max_bound(3)=32;
br_crossing(2).parameter.max_bound(3)=19.9;
br_crossing(3).parameter.max_bound(3)=40;
%
for i=1:3
figure(i)
hold on
clf;
br_crossing(i)=br_contn(mfuncs_ss,br_crossing(i),450);
br_crossing(i)=br_rvers(br_crossing(i));
br_crossing(i)=br_contn(mfuncs_ss,br_crossing(i),300);
end
%% Removing extra columns 
for i=1:3
br_crossing(i)=br_remove_extracolumns(br_crossing(i));
end

%% Computiong the stability
for i=1:3
[br_crossing_wbifs(i),nunst_pc{i},bif_px(i,:),pc_bif(:,i)]=MonitorChange(mfuncs_ss,br_crossing(i),...
    'range',2:length(br_crossing(i).point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
end
%%
save('br_crossing_threshold_try2.mat')
%%