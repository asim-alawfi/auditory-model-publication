%% Original Case: Computation of one-parameter computation in PR with fixing df=0.73, 
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
% upload from file 'Solve23_2df'
load('sol23_2df.mat') 
%% SetupFuncs and parameter bounds.
funcs_audi=set_symfuncs(@symbolic_auditory_with_symmetry_version,'sys_tau',@()[in.D, in.TD],...
   'sys_cond',@sys_cond);
parbd={'min_bound',[in.PR,0.5;in.df,0],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.01; 0,0.1],'print_residual_info',1};
%% Set up symmetry condition for dde-biftool implementation
% the system is symmetric with respect to the linear transformation:
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
% symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.25,2)'*[1,1],'condprojmat',[1,0,0,0,0,0]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
%% picking up a periodic solution from dde23 and continuing computation in PR (here df=0.73).
[funcs_s2,po2_symmetry,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df,in.PR,par1,...
    'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',{psolsym},'outputfuncs',true);
figure(2)
clf;
po2_symmetry=br_contn(funcs_s2,po2_symmetry,1000);
po2_symmetry=br_rvers(po2_symmetry);
po2_symmetry=br_contn(funcs_s2,po2_symmetry,1000);
po2_symmetry=br_rvers(po2_symmetry);
%% Computing stability and special points
[po2_symmetry_wbifs,nunst3_sym,po2_symmetry_bifs,po2_symmetry_bifind]=MonitorChange(funcs_s2,po2_symmetry,...
    'range',2:length(po2_symmetry.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
%% [nunst2_sym,dom,triv2_defect_sym,po2_symmetry.point]=GetStability(po2_symmetry,'funcs',funcs_s2,...
%     'exclude_trivial',true);%,'recompute',true);
chang2_stb_sym=find(diff(nunst3_sym));
rp2_x=arrayfun(@(x)x.parameter(in.PR),po2_symmetry_wbifs.point);
max2_y=arrayfun(@(x)max(x.profile(1,:)),po2_symmetry_wbifs.point);
min2_y2=arrayfun(@(x)min(x.profile(1,:)),po2_symmetry_wbifs.point);
% plotting stability
figure(2);clf;hold on;grid on
plot(rp2_x(nunst3_sym==0),max2_y(nunst3_sym==0),'bo',rp2_x(nunst3_sym>=1),max2_y(nunst3_sym>=1),'kx')%,'LineWidth',1)
plot(rp2_x(nunst3_sym==0),min2_y2(nunst3_sym==0),'bo',rp2_x(nunst3_sym>=1),min2_y2(nunst3_sym>=1),'kx')%,'LineWidth',1)
%% trace symmetric periodic orbit with maximum equal to the threshold (theta=0.5) in two parameters
% find extrema of $x$  along orbits on po_symmetry branch and pick orbits that
% have two extrema
smaxval=par(in.theta); % 
c_A=[1,0,0,0,0,0]; % c for u_A: to only extract the values of u_A
sympo_uA_extrema=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',po2_symmetry.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(po2_symmetry.point),sympo_uA_extrema);
[~,theta_cross]=min(abs(sympomax_ua-(smaxval+1e-4)));
%
second_max=3;
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
sympo0=setfield(po2_symmetry,'point',po2_symmetry.point(theta_cross));
sympo0.point.parameter([ine.t0,ine.val])=[sympo_uA_extrema{theta_cross}(second_max),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
[mfuncs,mbranch,suc_max]=ChangeBranchParameters(funcs_s2,sympo0,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{psolsym,max_cond},'outputfuncs',true,...
    'print_residual_info',1);
mbranch.parameter.max_bound(3)=40;
figure(3)
clf;
mbranch=br_contn(mfuncs,mbranch,3000);
mbranch=br_rvers(mbranch);
mbranch=br_contn(mfuncs,mbranch,3000);
[mbr_wbifs,dum,m_bifs,m_bifind]=MonitorChange(mfuncs,mbranch,'print_residual_info',0);
%[dum,m_dom,m_triv]=GetStability(mbr_wbifs,'exclude_trivial',true);
rp_thta=arrayfun(@(x)x.parameter(in.PR),mbr_wbifs.point);
df_thta=arrayfun(@(x)x.parameter(in.df),mbr_wbifs.point);
% plot stability for threshold crossing branch 
figure(3)
clf;hold on 
plot(rp_thta(dum==0),df_thta(dum==0),'g.',...
    rp_thta(dum>=1),df_thta(dum>=1),'r.','LineWidth',5)
grid on
xlabel('r_p')
ylabel('d_f')
%
%
%length(mbr_wbifs.point)
[~,it]=min(abs(rp_thta-34));
for i=1:1%length(mbr_wbifs.point)%1%m_bifind:-1:1
p=mbr_wbifs.point(it);
xm=p.mesh*p.period;
ym=p.profile(:,:);
figure(90)
clf;
plot(xm,ym(1:2,:),'LineWidth',2);
hold on
yline(0.5,'g--','LineWidth',2)
set(gca,'Fontweight','bold')
xlabel('t')
ylabel('Neurons activity')
legend('u_A','u_B')
grid on
%drawnow
end

%
po2_symmetry=br_remove_extracolumns(po2_symmetry);
po2_symmetry_wbifs=br_remove_extracolumns(po2_symmetry_wbifs);
mbranch=br_remove_extracolumns(mbranch);
mbr_wbifs=br_remove_extracolumns(mbr_wbifs);
%%
[po2_symmetry,nunst_pobr1,dom_pos,triv_defect_pos]=br_stabl(funcs_s2,po2_symmetry,0,1,'exclude_trivial',true);
%%
%
%save('branch_off_pos_original_case_and_threshold_crossing.mat')
%
save('branch_of_sympos_original_case_and_threshold_crossing_coj.mat')
%%
pt=po2_symmetry.point(120);

[r,J]=psolsym(pt,pt)

%%
 J1=J(1).profile(1,55:61)
J2=J(1).profile(2,289:301)
%%
s1=sum(J(1).profile,2)

%%
v=dde_coll_eva(pt.profile,pt.mesh,[0.1,0.6],pt.degree)
