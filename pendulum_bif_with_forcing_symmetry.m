clear;
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
 %% Bifurcation computation with enforcing symmetry in the pendulum model represented by  J.Sieber et al (2004)

 %% Set parameter bounds and initial condition
parnames={'a','b','tau'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
par([in.a, in.b, in.tau])=...
   [-1,   1.6,   1];
parbd={'min_bound',[in.a,-1;in.b,-4],'max_bound',[in.a,2; in.b,4.5],...
        'max_step',[in.a,0.05;in.b,0.1;0,0.05],'print_residual_info',1};
 funcs=set_symfuncs(@demo_pendulum_sym,'sys_tau',@()[in.tau]);
%% (1)  One-parameter branch for steady-state solutions with enforing symmetry
xini=[0;0];
xdim=length(xini);
Msym=-eye(length(xini));
pini=dde_stst_create('x',xini,'parameter',par);
symcond_stst=@(p,pref)dde_stst_lincond(p,xdim,'x','trafo',Msym,'shift',[1,2],'condprojint',linspace(0.1,0.25,2)'*[1,1]); 
[funcs_sym,stat_br,suceq]=SetupStst(funcs,'point',pini,'contpar',in.a,'step',0.01,parbd{:},...
                         'usercond',{symcond_stst},'outputfuncs',true,'extra_condition',true,...
                          'extracolumns','auto');
figure(1)
clf
ax1=gca;xlabel(ax1,'a');
stat_br=br_contn(funcs_sym,stat_br,300,'plotaxis',ax1);
%[stat_br,nunst_eq,demoeq,triv_defecteq]=br_stabl(funcs_sym,stat_br,0,0);
%%
[stat_br_wbifs,stattestfuncs,ind_bifeqs,bifeqs_typeseq]=LocateSpecialPoints(funcs_sym,stat_br);
[stat_br_wbifs,nunst_eq,demoeq,triv_defecteq]=br_stabl(funcs_sym,stat_br_wbifs,0,0);
chgsteq=find(diff(nunst_eq));
%% Plotting stability of steady-state solutions
%% One-parameter continuation for periodid orbits with enforcing symmetry
% braching off at a Hopf bifurcation, and continuing POs in one-paramete with enforcing symmetry
hopf_ind=ind_bifeqs(strcmp(bifeqs_typeseq,'hopf'));
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Msym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)'*[1,1]);

[pfuncs_sym,po_symbr,sucpo]=SetupPsol(funcs,stat_br_wbifs,hopf_ind,...
    'degree',6,'intervals',50,'print_residual_info',1,...
    'point.matrix','sparse','point.extra_condition',1,...
    'point.extracolumns','auto','max_step',[0,1;in.a,5e-3],...
    'continuation.plot_measure',{@(p)p.parameter(in.a),@(p)p.period},...
    'continuation.stops',{@(p)p(end).period>40},...
    'usercond',{psolsym},'outputfuncs',true);
figure(1)
hold on
po_symbr=br_contn(pfuncs_sym,po_symbr,200);
[po_symbr,nunst_po,demopo,triv_defectpo]=br_stabl(pfuncs_sym,po_symbr,0,0,'exclude_trivial',true);
%%
 [po_symbr_wbifs,po_nunst,sympobifs,sympobifind]=MonitorChange(pfuncs_sym,po_symbr,...
    'range',2:length(po_symbr.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
 %%
 clarray=colormap('lines');
 a_par2=arrayfun(@(x)x.parameter(in.a),stat_br_wbifs.point);
a_par=arrayfun(@(x)x.parameter(in.a),po_symbr_wbifs.point);
x_symstst=arrayfun(@(x)x.x(1),stat_br_wbifs.point);
x_sympo=arrayfun(@(x)max(x.profile(1,:)),po_symbr_wbifs.point);
y_sympo=arrayfun(@(x)min(x.profile(1,:)),po_symbr_wbifs.point);
figure(55)
clf
hold on; grid on
plot(a_par2(nunst_eq==0),x_symstst(nunst_eq==0),'-','MarkerSize',10,'LineWidth',2.5,'color',clarray(2,:));
plot(a_par2(nunst_eq>=1),x_symstst(nunst_eq>=1),'.','MarkerSize',8,'color',clarray(1,:));
%
plot(a_par(po_nunst==0),x_sympo(po_nunst==0),'.','MarkerSize',7,'color',clarray(5,:));
plot(a_par(po_nunst>=1),x_sympo(po_nunst>=1),'.','MarkerSize',7,'color',clarray(6,:));
plot(a_par(po_nunst==0),y_sympo(po_nunst==0),'.','MarkerSize',7,'color',clarray(5,:));
plot(a_par(po_nunst>=1),y_sympo(po_nunst>=1),'.','MarkerSize',7,'color',clarray(6,:));
plot(a_par2(chgsteq(1)),x_symstst(chgsteq(1)), 'p','Marker', 'p', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'black')
plot(a_par2(chgsteq(2)),x_symstst(chgsteq(2)),'s', 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'black')
plot(a_par(sympobifind(1)),x_sympo(sympobifind(1)), '^','Marker', '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')
plot(a_par(sympobifind(1)),y_sympo(sympobifind(1)), '^','Marker', '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black')
legend('stable sts','unstable sts','stable psol','unstable psol','','','fold','hopf',...
    'posl symbk')
xlabel('a','FontSize',16,'Interpreter','latex')
ylabel('$x_{1}$','FontSize',16,'Interpreter','latex')
grid on
set(gca,'FontWeight','bold')
   %   'MarkerS',15,'color',clarray(2,:));
%% (3) Continue symmetry-breaking in (a,b)-parameter space

pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Msym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.4,2)'*[1,1],'stateproj',eye(xdim),'condprojmat',eye(xdim));
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
[funcs_sym_bk,symbk_br,succ]=SetupPOEV1(funcs,po_symbr_wbifs,sympobifind(1),...
    'contpar',[in.a,in.b],'dir',in.a,parbd{:},...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
%
figure(8)
clf;
hold on
symbk_br=br_contn(funcs_sym_bk,symbk_br,300);
%%
symbk_br=br_remove_extracolumns(symbk_br);
[symbk_br,nunst_symbk,dombk,triv_defect_bk]=br_stabl(funcs_sym_bk,symbk_br,0,0,'exclude_trivial',1);
%%
 a_par=arrayfun(@(x)x.parameter(in.a),symbk_br.point);
 b_par=arrayfun(@(x)x.parameter(in.b),symbk_br.point);
figure(10)
clf; 
 hold on 
plot(a_par(nunst_symbk==0),b_par(nunst_symbk==0),'b-',...
    a_par(nunst_symbk>=1),b_par(nunst_symbk>=1),'k.','LineWidth',3)
grid on
legend('sym-breaking','')

xlabel('a','FontSize',16,'Interpreter','latex')
ylabel('b','FontSize',16,'Interpreter','latex')
grid on
set(gca,'FontWeight','bold')
%
stat_br_wbifs=br_remove_extracolumns(stat_br_wbifs);
stat_br=br_remove_extracolumns(stat_br);
po_symbr_wbifs=br_remove_extracolumns(po_symbr_wbifs);
po_symbr=br_remove_extracolumns(po_symbr);
%%
save('pundulum_bif_with_symmetry.mat')
%%
