clear
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
load('asymmetric_branches_case12_try2.mat')
%%
p_bifD=p_bif_3; %location of the bifurcation points  
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.05; in.df,0.01; 0,0.01],'print_residual_info',1};
[perdb_branch,suc]=DoublePsol(funcs_audi,asym_brs3_wbifs,p_bifD(2),'contpar',in.PR,'print_residual_info',true,parbds{:});
%%
figure(1)
clf
[perdb_branch,succ]=br_contn(funcs_audi,perdb_branch,300);
%%
[perdb_branch,nunst_db,dom,triv_defect]=br_stabl(funcs_audi,perdb_branch,0,0,'exclude_trivial',1);
chang_db=find(diff(nunst_db));
%% plotting one-par bif for symmetric POs in (df,max(u_A))-plan with stability
pr_db=arrayfun(@(x)x.parameter(in.PR),perdb_branch.point);
ymx_db=arrayfun(@(x)max(x.profile(1,:)),perdb_branch.point);
%%
[~,its]=min(abs(pr_db-5));
for i=1:1
    pi=perdb_branch.point(i);
    figure(27)
    clf
    hold on
    plot(pi.mesh,pi.profile(1:2,:),'LineWidth',2)
    yline(0.5,'k--','LineWidth',2)
    grid on
   ylim([0,1])
   xlim([0,1])
end
%%
figure(1)
clf;
hold on;grid on
plot(pr_db(nunst_db==0),ymx_db(nunst_db==0),'.',pr_db(nunst_db>=1),ymx_db(nunst_db>=1),'.','MarkerSize',10,'LineWidth',2)
figure(44)
hold on
plot(pr_db(nunst_db==0),'.','MarkerSize',10,'LineWidth',2)
plot(pr_db(nunst_db>=1),'.','MarkerSize',10,'LineWidth',2)
%%
save('period_doubling_bifPart1.mat')
%%
parbds={'min_bound',[in.PR,1;in.df,0.02],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.05; in.df,0.02; 0,0.02],'print_residual_info',1};
plm={@(p)p.parameter(asym_brs3_wbifs.parameter.free(1)),@(p)p.parameter(asym_brs3_wbifs.parameter.free(2))};
[pdfuncs,perid_branch,suc]=SetupPeriodDoubling(funcs_audi,asym_brs3_wbifs,p_bifD(3),'contpar',[in.PR,in.df],'dir',in.PR,...
    'print_residual_info',1,'print_residual_info',1,'step',0.02,'TorusInit.closest',[],parbds{:},plm);
%%
figure(3)
%clf;
perid_branch=br_contn(pdfuncs,perid_branch,200);
%%
perid_branch=br_rvers(perid_branch);
perid_branch=br_contn(pdfuncs,perid_branch,50);
%%
[perid_branch,nunst_dbbif,dombif,triv_defectbif]=br_stabl(pdfuncs,perid_branch,0,0,'exclude_trivial',1);
chang_dbbif=find(diff(nunst_dbbif));
pr_dbbif=arrayfun(@(x)x.parameter(in.PR),perid_branch.point);
df_dbbif=arrayfun(@(x)x.parameter(in.df),perid_branch.point);
save('period_doubling_bifPart2.mat')
%%
plm={@(p)p.parameter(asym_brs3_wbifs.parameter.free),@(p)p.parameter(asym_brs3_wbifs.parameter.free)};
[pdfuncs2,perid_branch2,sucd2]=SetupPeriodDoubling(funcs_audi,asym_brs3_wbifs,p_bifD(7),'contpar',[in.PR,in.df],'dir',in.PR,...
    'print_residual_info',1,'print_residual_info',1,'step',0.02,'TorusInit.closest',[],parbds{:},plm);
figure(3)
%clf;
perid_branch2=br_contn(pdfuncs2,perid_branch2,150);
%%
perid_branch2=br_rvers(perid_branch2);
perid_branch2=br_contn(pdfuncs2,perid_branch2,100);
%%
[perid_branch2,nunst_dbbif2,dombif2,triv_defectbif2]=br_stabl(pdfuncs2,perid_branch2,0,0,'exclude_trivial',1);

pr_dbbif2=arrayfun(@(x)x.parameter(in.PR),perid_branch2.point);
df_dbbif2=arrayfun(@(x)x.parameter(in.df),perid_branch2.point);
%%
figure(3)
clf;
hold on; grid on
plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'.','LineWidth',2)
plot(pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'.','LineWidth',2)
%%
rp3_Dx=arrayfun(@(x)x.parameter(in.PR),asym_brs3_wbifs.point);
df3_Dx=arrayfun(@(x)x.parameter(in.df),asym_brs3_wbifs.point);

figure(4)
clf;
hold on; grid on
plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'.','LineWidth',2)
plot(pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'.','LineWidth',2)
plot(rp3_Dx(nunst_3==0),df3_Dx(nunst_3==0),'.',rp3_Dx(nunst_3>=1),df3_Dx(nunst_3>=1),'.','MarkerSize',10,'LineWidth',2)
%
load('symmetry_breaking_D_b2_try2.mat')

load('touching_theta_case_increased_D_try2.mat')
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
%
figure(34)
clf
hold on 
plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'r.',...
     rp_b(nunst_sb>=1),df_b(nunst_sb>=1),'k.','LineWidth',1.5)
plot(rp_b2(uns==0),df_b2(uns==0),'r.',...
     rp_b2(uns>=1),df_b2(uns>=1),'k.','LineWidth',1.5)
plot(rp_bt2(nunst_pc3==0),df_bt2(nunst_pc3==0),'m--','LineWidth',2)
plot(rp_bt2(nunst_pc3>=1),df_bt2(nunst_pc3>=1),'.','Color',[0.7 0.7 0.7],'LineWidth',2)
plot(pr3_m(nunst3_mm==0),df3_m(nunst3_mm==0),'k.',...
    pr3_m(nunst3_mm>=1),df3_m(nunst3_mm>=1),'y.','LineWidth',2)
plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'b.',...
     pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'b.','LineWidth',2)
grid on
%% Computation of Torus bifrucation
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,40; in.df,1],...
    'max_step',[in.PR,0.05; in.df,0.02; 0,0.05],'print_residual_info',1};
[trfuncs,tr_branch,succ]=SetupTorusBifurcation(funcs_audi,br_crossing_wbifs(3),pc_bif(1,3), 'contpar',[in.PR,in.df],'dir',in.PR,...
    'step',0.01,parbds{:});%,'adapt_mesh_after_correct',0,'remesh',false);
%%
figure(34)
tr_branch=br_contn(trfuncs,tr_branch,300);
%%
tr_branch=br_rvers(tr_branch);
tr_branch=br_contn(trfuncs,tr_branch,190);


%% Compute stability 
[tr_branch,nunst_tr,dom_tr,triv_defect_tr]=br_stabl(trfuncs,tr_branch,0,0,'exclude_trivial',1);
save('period_doubling_bifPart3.mat')
%%

rp_btr=arrayfun(@(x)x.parameter(in.PR),tr_branch.point);
df_btr=arrayfun(@(x)x.parameter(in.df),tr_branch.point);
%

%
[~,it]=min(abs(rp_btr-34));
%find(diff(sign(rp_tr-34)))
px3=tr_branch.point(it);
px3_sol=trfuncs.get_comp(px3,'solution');
%px3.parameter(in.df)=px3.parameter(in.df)+0.2;
%px3.parameter(in.PR)=px3.parameter(in.PR)+0.2;
%
y_ini=px3_sol.profile(1:end,1);
his=@(t)dde_coll_eva(px3_sol.profile,px3_sol.mesh,1+t/px3_sol.period,px3_sol.degree);
tic
sol23_test=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,px3.parameter,in),lag,y_ini+0.05,[0,5],...
    ddeset('RelTol',1e-7,'AbsTol',1e-7,'Events',@(t,y,z)event(y)));
toc
%%
pr3_m0=pr3_m(nunst3_mm==0);
df3_m0=df3_m(nunst3_mm==0);
pr3_m1=pr3_m(nunst3_mm>=1);
df3_m1=df3_m(nunst3_mm>=1);
tr_loc=find(diff(nunst_pc3));
clrs=lines();
figure(99)
clf;
tiledlayout(8,6)
nexttile([8 3])
hold on 
d1=plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'-',rp_b2,df_b2,'-','Color',clrs(1,:),'LineWidth',3);
d4=plot(rp_bt2(nunst_pc3==0),df_bt2(nunst_pc3==0),'--','Color',clrs(7,:),'LineWidth',4);
d5=plot(rp_bt2(nunst_pc3>=1),df_bt2(nunst_pc3>=1),'.','Color',[0.5 0.5 0.5],'LineWidth',2);
d6=plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'-',...
    pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'-','Color',clrs(3,:),'LineWidth',4);
d7=plot(rp_btr(nunst_tr==0),df_btr(nunst_tr==0),'-','Color',clrs(5,:),'LineWidth',4);
d8=plot(rp_bt2(tr_loc(2)),df_bt2(tr_loc(2)),'o','MarkerSize',20,'Color',clrs(7,:),'LineWidth',2,'MarkerFaceColor', 'k') ;
d2=plot(pr3_m0(1:2:end),df3_m0(1:2:end),'k.','LineWidth',4,'MarkerSize',15);
d3=plot(pr3_m1(1:1:end),df3_m1(1:1:end),'.','Color',[0.7,0.7,0.7],'LineWidth',4,'MarkerSize',15);
%xlabel('$r_\mathrm{p}$','FontSize',18,'interpreter','latex')
title(' (a)','FontName','Courier',...
    'FontSize',26,'FontWeight','bold'),%'interpreter','latex',
legend([d1(1),d2,d3,d4,d5,d6(1),d7(1),d8(1)],...
    {'symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: unstable non-symmetric POs','touching-threshold: stable symmetric POs',...
    'touching-threshold: unstable symmetric POs','period-doubling ','torus bifurcation',' torus bifurcation'...
    },FontSize=18,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Courier',FontSize=26,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Courier',FontSize=26,FontWeight='bold')

grid on
ylim([0,1])
set(gca,'FontWeight','bold','LineWidth',2,'Box','on')
nexttile([4 3])
[~,its]=min(abs(pr_db-4.2));
pi=perdb_branch.point(its);
hold on
plot(pi.mesh*pi.period,pi.profile(1:2,:),'LineWidth',4)
plot(pi.mesh*pi.period,pi.profile(3:4,:),'--','LineWidth',2)
yline(0.5,'k--','Threshold','LineWidth',3)
xlabel('period','FontName','Courier',FontSize=22,FontWeight='bold')
%ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
title('(b)','FontName','Courier','FontSize',26,'FontWeight','bold')
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Courier',FontSize=28,FontWeight='bold')
set(gca, 'FontWeight','bold','LineWidth',2,'Box','on')
grid on
ylim([0,1])
xlim([0,pi.period])
grid on
nexttile([4 3])
plot(sol23_test.x,sol23_test.y(1:2,:),...
    sol23_test.x,sol23_test.y(3:4,:),'LineWidth',2);
yline(0.5,'k--','LineWidth',1)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','','interpreter','latex','FontName','Courier',FontSize=28,FontWeight='bold')
xlabel('time(s)','FontName','Courier',FontSize=22,FontWeight='bold')
grid on
title('(c)','FontName','Courier','FontSize',26,'FontWeight','bold')
set(gca,'FontWeight','bold','LineWidth',2,'Box','on')
%%
save('period_doubling_bifPart3_include_dde23simulation.mat')
%%
clear
load('period_doubling_bifPart3.mat')
rp_btr=arrayfun(@(x)x.parameter(in.PR),tr_branch.point);
df_btr=arrayfun(@(x)x.parameter(in.df),tr_branch.point);
[~,it]=min(abs(rp_btr-34));
%find(diff(sign(rp_tr-34)))
px3=tr_branch.point(it);
px3_sol=trfuncs.get_comp(px3,'solution');
px3.parameter(in.df)=px3.parameter(in.df)-0.003;
%px3.parameter(in.PR)=px3.parameter(in.PR)+0.2;

y_ini=px3_sol.profile(1:end,1);
his=@(t)dde_coll_eva(px3_sol.profile,px3_sol.mesh,1+t/px3_sol.period,px3_sol.degree);
tic
sol23_test2=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,px3.parameter,in),lag,y_ini+0.02,[0,10],...
    ddeset('RelTol',1e-7,'AbsTol',1e-7,'Events',@(t,y,z)event(y)));
figure(111)
plot(sol23_test2.x,sol23_test2.y(1:2,:),...
    sol23_test2.x,sol23_test2.y(3:4,:),'LineWidth',1);
yline(0.5,'k--','LineWidth',1)
xlabel('time(s)','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
grid on
title('(c)','interpreter','latex','FontSize',16)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
stab_location=find(diff(nunst3_mm));
sab_loc_sym=find(diff(nunst_pc3));
figure(77)
clf;
tiledlayout(8,3,'TileSpacing','compact')
nexttile([8 2])
hold on 
d1=plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'-',rp_b2,df_b2,'-','Color',clrs(1,:),'LineWidth',7);
d4=plot(rp_bt2(sab_loc_sym(2):sab_loc_sym(3)),df_bt2(sab_loc_sym(2):sab_loc_sym(3)),'-','Color','m','LineWidth',6);
d5=plot(rp_bt2(1:sab_loc_sym(2)),df_bt2(1:sab_loc_sym(2)),'-','Color',[0.5 0.5 0.5],'LineWidth',6);
plot(rp_bt2(sab_loc_sym(3):end),df_bt2(sab_loc_sym(3):end),'-','Color',[0.5 0.5 0.5],'LineWidth',6);
d6=plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'-',...
    pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'-','Color',clrs(3,:),'LineWidth',6);
d7=plot(rp_btr(nunst_tr==0),df_btr(nunst_tr==0),'-','Color',clrs(5,:),'LineWidth',5);
d8=plot(rp_bt2(tr_loc(2)),df_bt2(tr_loc(2)),'o','MarkerSize',20,'Color',clrs(7,:),'LineWidth',2,'MarkerFaceColor', 'k') ;
d2=plot(pr3_m(1:stab_location(1)),df3_m(1:stab_location(1)),'k-','LineWidth',6)
plot(pr3_m(stab_location(3):end),df3_m(stab_location(3):end),'k-','LineWidth',6)
plot(pr3_m(stab_location(2):stab_location(3)),df3_m(stab_location(2):stab_location(3)),'k-','LineWidth',6)

d3=plot(pr3_m(stab_location(1):stab_location(2)),df3_m(stab_location(1):stab_location(2)),'-','Color',[0.5,0.5,0.5],'LineWidth',4,'MarkerSize',15);
%xlabel('$r_\mathrm{p}$','FontSize',18,'interpreter','latex')
%title(' (a)','FontName','Cambria','FontSize',26,'FontWeight','normal')%'interpreter','latex',
legend([d1(1),d2,d3,d4,d5,d6(1),d7(1),d8(1)],...
    {'symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: unstable non-symmetric POs','touching-threshold: stable symmetric POs',...
    'touching-threshold: unstable symmetric POs','period-doubling ','torus bifurcation',' torus bifurcation'...
    },FontSize=28,FontWeight='normal',FontName='Cambria')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Cambria',FontSize=26,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Cambria',FontSize=26,FontWeight='bold')
xlim([0,40])
grid on
ylim([0,1])
set(gca,'FontWeight','bold','LineWidth',2,'Box','on')
nexttile([2 1])
[~,its]=min(abs(pr_db-4.2));
pi=perdb_branch.point(its);
hold on
plot(pi.mesh*pi.period,pi.profile(1:2,:),'LineWidth',4)
plot(pi.mesh*pi.period,pi.profile(3:4,:),'--','LineWidth',2)
yline(0.5,'k--','Threshold','LineWidth',3)
xlabel('period','FontName','Cambria',FontSize=22,FontWeight='bold')
%title('(b)','FontName','Cambria','FontSize',26,'FontWeight','normal'),
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Courier',FontSize=28,FontWeight='bold')
set(gca, 'FontWeight','bold','LineWidth',2,'Box','on')
grid on
ylim([0,1])
xlim([0,pi.period])
grid on
nexttile([2 1])
plot(sol23_test.x,sol23_test.y(1:2,:),...
    sol23_test.x,sol23_test.y(3:4,:),'LineWidth',2);
yline(0.5,'k--','LineWidth',1)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','',....
    'interpreter','latex','FontName','Courier',FontSize=28,FontWeight='bold')
xlabel('time(s)','FontName','Cambria',FontSize=22,FontWeight='bold')
grid on
%title('(c)','FontName','Cambria','FontSize',26,'FontWeight','normal')
set(gca,'FontWeight','bold','LineWidth',2,'Box','on')
%%%%%%%%%%
nexttile([4 1])
hold on
plot(pr_m,df_m,'-','Color',clrs2(6,:),'LineWidth',6)
plot(pr3_m,df3_m,'k-','LineWidth',5)
%legend({sprintf('PB: (t_{d},D)=(0.022,0.015)'),....
  %  sprintf('PB: (t_{d},D)=(0.022,0.05)'),},...
    %'Fontname','Courier','FontSize',26)
xlabel('$r_\mathrm{p}$','FontSize',30,'interpreter','latex','FontName','Courier')
ylabel('$d_\mathrm{f}$','FontSize',30,'interpreter','latex','FontName','Courier') 
text(10, 0.85, {'Comparison of', 'perceptual boundary'}, 'FontSize', 20, 'Fontname','Cambria','BackgroundColor', 'non', 'EdgeColor', 'non','FontWeight','bold');
%title('(d)','FontName','Cambria','FontSize',26,'FontWeight','normal')

grid on 
set(gca,'FontWeight','bold', 'LineWidth',2,'Box','on')
