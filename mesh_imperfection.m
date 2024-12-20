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

%% We pick a periodic solution using dde23, then we continue one-parameter continuation for POs in $PR$
  load('sol23_DF.mat')
%%  (1) Computation without using symmetry condtion.
%funcs_audi=set_symfuncs(@symbolic_auditory_with_symmetry_version,'sys_tau',@()[in.D, in.TD],'sys_cond',@sys_cond);
funcs_audi=set_symfuncs(@symbolic_auditory_with_symmetry_version,'sys_tau',@()[in.D, in.TD],'sys_cond',@sys_cond);
parbd={'min_bound',[in.PR,15;in.df,0],'max_bound',[in.PR,18; in.df,1],...
    'max_step',[in.PR,0.1; in.df,0.01; 0,0.1],'print_residual_info',1};
[funcs_symb,per_br_nosymcond,suc]=branch_from_sol_audi_sym(funcs_audi,dde23_sols,in.PR,par,...
    'extra_condition',true,'phase_condition',0,parbd{:},...
    'degree',4,'intervals',100,'outputfuncs',true);
per_br_nosymcond.method.point.newton_max_iterations=10;
per_br_nosymcond.method.point.adapt_mesh_after_correct=1;
per_br_nosymcond.method.point.matrix='sparse';
per_br_nosymcond.method.stability.matrix='sparse';
figure(100)
clf;
per_br_nosymcond=br_contn(funcs_symb,per_br_nosymcond,100);
%%
[per_br_nosymcond,unst_stab_prs,domm,triv_defectt]=br_stabl(funcs_symb,per_br_nosymcond,0,0,'exclude_trivial',1);
% %chang_stabs=find(diff(unst_stab_prs));
 %% plot stability
PRsy_ax=arrayfun(@(x)x.parameter(in.PR),per_br_nosymcond.point);
ysy_ax=arrayfun(@(x)max(x.profile(1,:)),per_br_nosymcond.point);
y2sy_ax=arrayfun(@(x)min(x.profile(1,:)),per_br_nosymcond.point);
 %% (2) Computation when  symmetry condition is implemented.
 Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
% symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.45,2)'*[1,1]);
% psolsymfine=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
%     'condprojint',linspace(0,1,100)'*[1,1]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
[funcs_symb_sm,per_br_withsymcond,suc2]=branch_from_sol_audi_sym(funcs_audi,dde23_sols,in.PR,par,...
    'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',4,'intervals',100,'extracolumns','auto','newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',{psolsym},'outputfuncs',true);
figure(100)
hold
per_br_withsymcond=br_contn(funcs_symb_sm,per_br_withsymcond,100);
%%
[per_br_withsymcond,nunst_sym3,dom,triv_defect]=br_stabl(funcs_symb_sm,per_br_withsymcond,0,0,'exclude_trivial',1);
PRsycond_ax=arrayfun(@(x)x.parameter(in.PR),per_br_withsymcond.point);
ysycond_ax=arrayfun(@(x)max(x.profile(1,:)),per_br_withsymcond.point);
y2sycond_ax=arrayfun(@(x)min(x.profile(1,:)),per_br_withsymcond.point);
per_br_withsymcond=br_remove_extracolumns(per_br_withsymcond)
%%
clrs=lines();
bif_loc=find(diff(nunst_sym3));
%%
pt1=per_br_nosymcond.point(1);
pt2=per_br_nosymcond.point(bif_loc-1);
%
branch=per_br_nosymcond;
figure(1)
clf;
tiledlayout(6,5)
nexttile([6,3])
hold on; grid on
plt1=plot(PRsy_ax(1:bif_loc),ysy_ax(1:bif_loc),'.-','Color',clrs(1,:),'MarkerSize',30);
plt2=plot(PRsy_ax(bif_loc:end),ysy_ax(bif_loc:end),'.-','Color',clrs(5,:),'MarkerSize',30);
plt3=plot(PRsycond_ax(nunst_sym3>=1),ysycond_ax(nunst_sym3>=1),'ko-','LineWidth',2);
plt4=plot(PRsy_ax(bif_loc),ysy_ax(bif_loc),'.','Color','red','MarkerSize',50);
plt_text={'symmetric POs: stable','non-symmetric POs: stable','symmetric POs: unstable','symmetry-breaking'};
plt_vec=[plt1(1);plt2(1);plt3(1);plt4(1)];
legend(plt_vec,plt_text,'Interpreter','latex','EdgeColor',0.5*[1,1,1],'FontSize',24)
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier', 'LineWidth',2)
nexttile([2,2])
err_plot_half(pt1,gca);
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier', 'LineWidth',2)
nexttile([2,2])
err_plot_half(pt2,gca);
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier', 'LineWidth',2)
nexttile([2,2])
[esym,msym]=arrayfun(@(p)err_plot_half(p),branch.point);
periods=arrayfun(@(p)p.period,branch.point);
pe=semilogy(gca,1:length(periods),esym,'LineWidth',3,...
    'DisplayName','$\Delta_{\hat{e}}$');
hold(gca,'on');
pm=semilogy(gca,1:length(periods),msym,'LineWidth',2,...
    'DisplayName','$\Delta_{\hat{m}}$');
ylim(gca,[1e-5,1e-1]);
xlabel(gca,'point number along branch','Interpreter','latex');
legend(gca,[pe,pm],'Interpreter','latex')
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold',...
    'FontName','Courier')
%%
mbrown= [0.65, 0.33, 0.13];
pclr=[0.75, 0, 0.75];
figure(505);clf;
tiledlayout(6,7,"TileSpacing","compact");
nexttile([6,3]);
hold on;grid on
plt1=plot(PRsy_ax(1:bif_loc),ysy_ax(1:bif_loc),'.-','Color',clrs(1,:),'MarkerSize',30);
plt2=plot(PRsy_ax(bif_loc:end),ysy_ax(bif_loc:end),'.-','Color',clrs(5,:),'MarkerSize',30);
plt3=plot(PRsycond_ax(nunst_sym3>=1),ysycond_ax(nunst_sym3>=1),'ko-','LineWidth',2);
plt4=plot(PRsy_ax(bif_loc),ysy_ax(bif_loc),'^','Color','k','MarkerFaceColor', clrs(3,:),'MarkerSize',20);
plt5=plot(PRsy_ax(bif_loc-1),ysy_ax(bif_loc-1),'^','Color','k','MarkerFaceColor', clrs(7,:),'MarkerSize',20);
plt6=plot(PRsy_ax(1),ysy_ax(1),'^','Color','k','MarkerFaceColor','red','MarkerSize',20);
plt_text={'symmetric POs: stable','non-symmetric POs: stable','symmetric POs: unstable','symmetry-breaking'};
plt_vec=[plt1(1);plt2(1);plt3(1);plt4(1)];
legend(plt_vec,plt_text,'Interpreter','latex','EdgeColor',0.5*[1,1,1],'FontSize',24)
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier', 'LineWidth',2)
title('(a)')
%%%%%%%%%%%%%
nexttile([2,2])
hold on
grid on
plot(pt1.mesh,pt1.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{A}$','Interpreter','latex','FontSize',16,'Location','Best')
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b): solution of the red rectangle')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier')
%%%%%%%%%%%%%%%
nexttile([2,2])
hold on
grid on
plot(pt2.mesh,pt2.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{A}$','Interpreter','latex','FontSize',16,'Location','Best')
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b): solution of the red rectangle')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier')
% %%%%%%%%%%%
nexttile([2,2])
ntfine=100001;
tcoarse=pt1.mesh(1:pt1.degree:end);
nt=length(tcoarse);
m_int=interp1(tcoarse,linspace(0,1,nt),'linear','pp');
tfine=linspace(pt1.mesh(1),pt1.mesh(end),ntfine);
t1=tfine(tfine<0.5);
t2=tfine(tfine>=0.5);
m=fnder(m_int);
hold on; grid on 
plot(t1,fnval(m,t1),'-','Color',clrs(5,:),'LineWidth',3)
plot(t1,fnval(m,t2(1:end-1)),'-','Color',clrs(7,:),'LineWidth',2)
yline(1,'k--','LineWidth',2)
xlabel('period','FontSize',12,'FontName','Cambria')
text(0.4,1.35,'$y=1$','Interpreter','latex','FontSize',16,'FontWeight','bold')
title('(c)')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier','box','off')
legend('$m(t)$','$m(t+0.5)$','Interpreter','latex')
xlim([0,0.5])
% %%%%%%%
nexttile([2,2])
tcoarsept2=pt2.mesh(1:pt2.degree:end);
ntpt2=length(tcoarsept2);
m_intpt2=interp1(tcoarsept2,linspace(0,1,ntpt2),'linear','pp');
tfinept2=linspace(pt2.mesh(1),pt2.mesh(end),ntfine);
t1pt2=tfinept2(tfinept2<0.5);
t2pt2=tfinept2(tfinept2>=0.5);
mpt2=fnder(m_intpt2);
hold on; grid on 
plot(t1pt2,fnval(mpt2,t1pt2),'-','Color',clrs(5,:),'LineWidth',3)
plot(t1pt2,fnval(mpt2,t2pt2(1:end-1)),'-','Color',clrs(7,:),'LineWidth',2)
yline(1,'k--','LineWidth',2)
xlabel('period','FontSize',12,'FontName','Cambria')
text(0.4,1.35,'$y=1$','Interpreter','latex','FontSize',16,'FontWeight','bold')
title('(c)')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier','box','off')
legend('$m(t)$','$m(t+0.5)$','Interpreter','latex')
xlim([0,0.5])
%
nexttile([2,2])
hold on; grid on
yline(0,'k--','LineWidth',3)
plot(t1,fnval(m,t1)-fnval(m,t2(1:end-1)),'LineWidth',2)
xlim([0,0.5])
xlabel('period','FontSize',16,'FontName','Cambria')
legend('','$m(t)-m(t+0.5)$','Interpreter','latex')
title('(d)')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier')
%%
nexttile([2,2])
hold on; grid on
yline(0,'k--','LineWidth',3)
plot(t1pt2,fnval(mpt2,t1pt2)-fnval(mpt2,t2pt2(1:end-1)),'LineWidth',2)
xlim([0,0.5])
xlabel('period','FontSize',16,'FontName','Cambria')
legend('','$m(t)-m(t+0.5)$','Interpreter','latex')
title('(d)')
set(gca,'FontSize',16,'FontWeight','bold','FontName','Courier')
%%%%%%%%

%%

save('imperfection_in_mesh.mat')
