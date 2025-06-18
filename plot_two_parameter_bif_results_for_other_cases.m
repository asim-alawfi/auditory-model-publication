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
load('symmetry_breaking_TD_try2.mat')
load('touching_theta_case_increased_TD_try2.mat')
%%
graycolor=[0.7 0.7 0.7];
mbrown= [0.65, 0.33, 0.13];
clrs=lines();
nunst_pct1=nunst_pc{1};
rp_theta1=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(1).point);
df_theta1=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(1).point);
%
set_bif={'LineWidth',2,'Box','on','FontSize',10,'FontWeight','bold'};
figure(13)
clf;
hold on
plot(rp_btd(nunst_td==0),df_btd(nunst_td==0),'-','Color',clrs(1,:),'LineWidth',3)

plot(rp_theta1(nunst_pct1==0),df_theta1(nunst_pct1==0),'m.','LineWidth',1.5)
plot(rp_theta1(nunst_pct1>=1),df_theta1(nunst_pct1>=1),'.','Color',[0.7 0.7 0.7],'LineWidth',1.5)
%plot(pr1_m,df1_m,'k--','LineWidth',2)
plot(pr1_m(nunst1_mm==0),df1_m(nunst1_mm==0),'k.',...
    pr1_m(nunst1_mm>=1),df1_m(nunst1_mm>=1),'r.','LineWidth',2)
set(gca, set_bif{:})
yticks([0.2,0.4,0.6,0.8,1])
title('(a)','FontSize',14)
grid on
ylim([0,1])
%%
load('symmetry_breaking_TD_and_D_try2.mat')
load('touching_theta_case_increased_TD_and_D_try2.mat')
load('foldbifrucation_case_TD_D_results.mat')
%
nunst_pct2=nunst_pc{2};
rp_theta2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(2).point);
df_theta2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(2).point);
rp_fld=arrayfun(@(x)x.parameter(in.PR),pfbranch.point);
df_fld=arrayfun(@(x)x.parameter(in.df),pfbranch.point);
%%
figure(23)
clf;
hold on
plot(rp_btdd(nunst_tdd==0),df_btdd(nunst_tdd==0),'-','Color',clrs(1,:),'LineWidth',2)
%plot(rp_btdd(nunst_tdd>=1),df_btdd(nunst_tdd>=1),'r.','LineWidth',2)
plot(rp_theta2(nunst_pct2==0),df_theta2(nunst_pct2==0),'m-','LineWidth',2)
plot(rp_theta2(nunst_pct2>=1),df_theta2(nunst_pct2>=1),'-','Color',[0.7 0.7 0.7],'LineWidth',1.5)
plot(pr2_m(nunst2_mm==0),df2_m(nunst2_mm==0),'k--','LineWidth',2)
%pr2_m(nunst2_mm>=1),df2_m(nunst2_mm>=1),'y.','LineWidth',2)
plot(rp_fld,df_fld,'-','Color',clrs(3,:),'LineWidth',2)
set(gca, set_bif{:})
legend('symmetry-breaking','POs with touching \theta =0.5','','','fold of PO','FontSize',16)
title('(b)','FontSize',14)
ylim([0,1])
yticks([0.2,0.4,0.6,0.8,1])
grid on
%%
load('symmetry_breaking_D_b2_try2.mat')

load('touching_theta_case_increased_D_try2.mat')
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
changes_in_stab=find(diff(nunst3_mm));
figure(34)
clf
hold on
plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'r-',rp_b2(uns==0),df_b2(uns==0),'r-','LineWidth',2)
plot(rp_bt2(nunst_pc3==0),df_bt2(nunst_pc3==0),'m--','LineWidth',2)
plot(rp_bt2(nunst_pc3>=1),df_bt2(nunst_pc3>=1),'.','Color',[0.7 0.7 0.7],'LineWidth',2)
%plot(pr3_m(nunst3_mm==0),df3_m(nunst3_mm==0),'k.',...
plot(pr3_m(1:changes_in_stab(1)),df3_m(1:changes_in_stab(1)),'k-','LineWidth',2)
plot(pr3_m(changes_in_stab(2):end),df3_m(changes_in_stab(2):end),'k-','LineWidth',2)
plot(pr3_m(changes_in_stab(1):changes_in_stab(2)),df3_m(changes_in_stab(1):changes_in_stab(2)),'b-','LineWidth',2)
grid on
ylim([0,1])
set(gca, set_bif{:})
%%

%
load('period_doubling_bifPart3_include_dde23simulation.mat')
%%
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
changes_in_stab=find(diff(nunst3_mm));
pr3_m0=pr3_m(nunst3_mm==0);
df3_m0=df3_m(nunst3_mm==0);
pr3_m1=pr3_m(nunst3_mm>=1);
df3_m1=df3_m(nunst3_mm>=1);
tr_loc=find(diff(nunst_pc3));
chang_stab_touch_sym=find(diff(nunst_pc3));
clrs=lines();
set_bif={'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold'};
figure(99)
clf;
tiledlayout(8,6)
nexttile([8 3])
hold on
d1=plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'-',rp_b2,df_b2,'-','Color',clrs(1,:),'LineWidth',3);
d4=plot(rp_bt2(chang_stab_touch_sym(2):chang_stab_touch_sym(3)),df_bt2(chang_stab_touch_sym(2):chang_stab_touch_sym(3)),'-','Color','m','LineWidth',3);
d5=plot(rp_bt2(1:chang_stab_touch_sym(2)),df_bt2(1:chang_stab_touch_sym(2)),'-','Color',[0.5 0.5 0.5],'LineWidth',3)
plot(rp_bt2(chang_stab_touch_sym(3):end),df_bt2(chang_stab_touch_sym(3):end),'-','Color',[0.5 0.5 0.5],'LineWidth',3);
d6=plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'-',...
    pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'-','Color',clrs(3,:),'LineWidth',3);
d7=plot(rp_btr(nunst_tr==0),df_btr(nunst_tr==0),'-','Color',clrs(5,:),'LineWidth',4);
d8=plot(rp_bt2(tr_loc(2)),df_bt2(tr_loc(2)),'o','MarkerSize',20,'Color',clrs(7,:),'LineWidth',2,'MarkerFaceColor', 'k') ;
d2=plot(pr3_m(1:changes_in_stab(1)),df3_m(1:changes_in_stab(1)),'k--','LineWidth',3)
plot(pr3_m(changes_in_stab(2):end),df3_m(changes_in_stab(2):end),'k--','LineWidth',3);%plot(pr3_m0(1:2:end),df3_m0(1:2:end),'k.','LineWidth',4,'MarkerSize',15);
d3=plot(pr3_m(changes_in_stab(1):changes_in_stab(2)),df3_m(changes_in_stab(1):changes_in_stab(2)),'--','Color',[0.5 0.5 0.5],'LineWidth',2)%plot(pr3_m1(1:1:end),df3_m1(1:1:end),'.','Color',[0.7,0.7,0.7],'LineWidth',4,'MarkerSize',15);
set(gca, set_bif{:})
legend([d1(1),d2,d3,d4,d5,d6(1),d7(1),d8(1)],...
    {'symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: unstable non-symmetric POs','touching-threshold: stable symmetric POs',...
    'touching-threshold: unstable symmetric POs','period-doubling ','torus bifurcation',' torus bifurcation',...
    },'FontSize',18,'FontWeight','bold','Interpreter','latex')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=32,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=32,FontWeight='bold')
grid on
xlim([0,40])
ylim([0,1])
title(' (a)', 'FontSize',16,'FontWeight','bold'),%'interpreter','latex',
legend([d1(1),d2,d3,d4,d5,d6(1),d7(1),d8(1)],...
    {'symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: unstable non-symmetric POs','touching-threshold: stable symmetric POs',...
    'touching-threshold: unstable symmetric POs','period-doubling ','torus bifurcation',' torus bifurcation'...
    },FontSize=18,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=26,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=26,FontWeight='bold')
grid on
yticks([0.2,0.4,0.6,0.8,1])
nexttile([4 3])
[~,its]=min(abs(pr_db-4.2));
pi=perdb_branch.point(its);
hold on
plot(pi.mesh*pi.period,pi.profile(1:2,:),'LineWidth',4)
plot(pi.mesh*pi.period,pi.profile(3:4,:),'--','LineWidth',2)
yline(0.5,'k--','Threshold','LineWidth',3)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Courier',FontSize=28,FontWeight='bold')
set(gca, set_bif{:})
xlabel('period',FontSize=18,FontWeight='bold')
%ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
title('(b)','FontSize',16)
grid on
ylim([0,1])
xlim([0,pi.period])
yticks([0.2,0.4,0.6,0.8,1])
grid on
nexttile([4 3])
plot(sol23_test.x,sol23_test.y(1:2,:),...
    sol23_test.x,sol23_test.y(3:4,:),'LineWidth',2);
yline(0.5,'k--','LineWidth',1)
grid on
set(gca,set_bif{:})
title('(c)','FontSize',16)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','','interpreter','latex',FontSize=40,FontWeight='bold')
xlabel('time(s)',FontSize=18,FontWeight='bold')
yticks([0.2,0.4,0.6,0.8,1])
%%
set_bif={'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold'};
chang_stab_touch_asym1=find(diff(nunst1_mm));
chang_stab_touch_sym1=find(diff(nunst_pct1));
figure(18)
clf;
tiledlayout(1,2,"TileSpacing","compact")
nexttile
hold on
plot(rp_btd,df_btd,'-','Color',clrs(1,:),'LineWidth',3)
plot(rp_theta1(chang_stab_touch_sym1(2):chang_stab_touch_sym1(3)),df_theta1(chang_stab_touch_sym1(2):chang_stab_touch_sym1(3)),'m-','LineWidth',3)
plot(rp_theta1(1:chang_stab_touch_sym1(2)),df_theta1(1:chang_stab_touch_sym1(2)),'.','Color',0.5*[1 1 1],'LineWidth',2)
plot(rp_theta1(chang_stab_touch_sym1(3):end),df_theta1(chang_stab_touch_sym1(3):end),'-','Color',0.5*[1 1 1],'LineWidth',3)

plot(pr1_m(1:chang_stab_touch_asym1(3)),df1_m(1:chang_stab_touch_asym1(3)),'k--','LineWidth',3)
plot(pr1_m(chang_stab_touch_asym1(end):end),df1_m(chang_stab_touch_asym1(end):end),'k--','LineWidth',3)
plot(pr1_m(chang_stab_touch_asym1(3):chang_stab_touch_asym1(4)),df1_m(chang_stab_touch_asym1(3):chang_stab_touch_asym1(4)),'-','Color',0.5*[1 1 1],'LineWidth',3)
% plot(pr1_m(nunst1_mm==0),df1_m(nunst1_mm==0),'k--',...
%     pr1_m(nunst1_mm>=1),df1_m(nunst1_mm>=1),'k--','LineWidth',2)
set(gca, set_bif{:})
yticks([0.2,0.4,0.6,0.8,1])
title('(a)','FontSize',20)
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=26,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=26,FontWeight='bold')
grid on
ylim([0,1])
nexttile
hold on
plot(rp_btdd(nunst_tdd==0),df_btdd(nunst_tdd==0),'-','Color',clrs(1,:),'LineWidth',3)
%plot(rp_btdd(nunst_tdd>=1),df_btdd(nunst_tdd>=1),'r.','LineWidth',2)
plot(pr2_m(nunst2_mm==0),df2_m(nunst2_mm==0),'k--','LineWidth',3)
plot(rp_theta2(nunst_pct2==0),df_theta2(nunst_pct2==0),'m-','LineWidth',3)
plot(rp_theta2(nunst_pct2>=1),df_theta2(nunst_pct2>=1),'-','Color',0.5*[1 1 1],'LineWidth',1.5)

%pr2_m(nunst2_mm>=1),df2_m(nunst2_mm>=1),'y.','LineWidth',2)
plot(rp_fld,df_fld,'-','Color',clrs(7,:),'LineWidth',3)

% ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=26,FontWeight='bold')
set(gca, set_bif{:})
legend('symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: stable symmetric POs','touching-threshold: unstable POs',...
    'fold of PO','FontSize',20)
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=26,FontWeight='bold')
title('(b)','FontSize',20)
ylim([0,1])
yticks([0.2,0.4,0.6,0.8,1])
grid on
%%
load('period_doubling_bifPart3_include_dde23simulation.mat')
%%
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
changes_in_stab=find(diff(nunst3_mm));
pr3_m0=pr3_m(nunst3_mm==0);
df3_m0=df3_m(nunst3_mm==0);
pr3_m1=pr3_m(nunst3_mm>=1);
df3_m1=df3_m(nunst3_mm>=1);
tr_loc=find(diff(nunst_pc3));
chang_stab_touch_sym=find(diff(nunst_pc3));
ft=13;
set_bif={'LineWidth',2,'Box','on','FontSize',12}%,'FontWeight','bold'};
clrs=lines();
figure(99)
clf;
tiledlayout(2,3,'TileSpacing','compact')
nexttile([2 2])
hold on
d1=plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'-',rp_b2,df_b2,'-','Color',clrs(1,:),'LineWidth',3);
d4=plot(rp_bt2(chang_stab_touch_sym(2):chang_stab_touch_sym(3)),df_bt2(chang_stab_touch_sym(2):chang_stab_touch_sym(3)),'-','Color','m','LineWidth',3);
d5=plot(rp_bt2(1:chang_stab_touch_sym(2)),df_bt2(1:chang_stab_touch_sym(2)),'-','Color',[0.5 0.5 0.5],'LineWidth',3)
plot(rp_bt2(chang_stab_touch_sym(3):end),df_bt2(chang_stab_touch_sym(3):end),'-','Color',[0.5 0.5 0.5],'LineWidth',3);
d6=plot(pr_dbbif(nunst_dbbif==0),df_dbbif(nunst_dbbif==0),'-',...
    pr_dbbif2(nunst_dbbif2==0),df_dbbif2(nunst_dbbif2==0),'-','Color',clrs(3,:),'LineWidth',3);
d7=plot(rp_btr(nunst_tr==0),df_btr(nunst_tr==0),'-','Color',clrs(5,:),'LineWidth',4);
%d8=plot(rp_bt2(tr_loc(2)),df_bt2(tr_loc(2)),'o','MarkerSize',20,'Color',clrs(7,:),'LineWidth',2,'MarkerFaceColor', 'k') ;
d2=plot(pr3_m(1:changes_in_stab(1)),df3_m(1:changes_in_stab(1)),'k--','LineWidth',3)
plot(pr3_m(changes_in_stab(2):end),df3_m(changes_in_stab(2):end),'k--','LineWidth',3);%plot(pr3_m0(1:2:end),df3_m0(1:2:end),'k.','LineWidth',4,'MarkerSize',15);
d3=plot(pr3_m(changes_in_stab(1):changes_in_stab(2)),df3_m(changes_in_stab(1):changes_in_stab(2)),'--','Color',[0.5 0.5 0.5],'LineWidth',2)%plot(pr3_m1(1:1:end),df3_m1(1:1:end),'.','Color',[0.7,0.7,0.7],'LineWidth',4,'MarkerSize',15);
[~,its]=min(abs(pr_db-4.2));
pi=perdb_branch.point(its);
plot(pi.parameter(in.PR),pi.parameter(in.df),'kx','MarkerSize',8,'LineWidth',2)
text(pi.parameter(in.PR)-1.25,pi.parameter(in.df)+0.01 ,'B', 'FontSize', 11, 'Color', 0.0*[1 1 1],'FontWeight','bold')
plot(px3.parameter(in.PR),px3.parameter(in.df),'kx','MarkerSize',8,'LineWidth',2)
text(px3.parameter(in.PR)-1,px3.parameter(in.df)+0.03 ,'c', 'FontSize', 13, 'Color', 0.0*[1 1 1],'FontWeight','bold')
set(gca, set_bif{:})
% legend([d1(1),d2,d3,d4,d5,d6(1),d7(1),d8(1)],...
%     {'symmetry-breaking','touching-threshold: stable non-symmetric POs','touching-threshold: unstable non-symmetric POs','touching-threshold: stable symmetric POs',...
%     'touching-threshold: unstable symmetric POs','period-doubling ','torus bifurcation',' torus bifurcation'...
%     },FontSize=30,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=15,FontWeight='normal')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=15,FontWeight='normal')

xlim([0,40])
ylim([0,1])
title(' A', 'FontSize',ft,'FontWeight','normal'),%'interpreter','latex',
legend([d1(1),d2,d3,d4,d5,d6(1),d7(1)],...
    {'sym-breaking','touching \theta: stable non-sym','touching \theta: unstable non-sym','touching \theta: stable sym',...
    'touching \theta: unstable sym','period-doubling ',' quasi-periodic bif'},FontSize=ft,FontWeight='normal')
xticks([5,10,15,20,25,30,35,40])
yticks([0,0.2,0.4,0.6,0.8,1])
nexttile
hold on
plot(pi.mesh*pi.period,pi.profile(1:2,:),'LineWidth',2)
plot(pi.mesh*pi.period,pi.profile(3:4,:),'LineWidth',2)
yline(0.5,'k--','LineWidth',3)
set(gca, set_bif{:})
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','','interpreter','latex',FontSize=ft,FontWeight='normal')
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('B','FontSize',ft,'FontWeight','normal')

ylim([0,1])
xlim([0,pi.period])
%xticks([0.2,0.5,1])
yticks([0,0.5,1])
nexttile
plot(sol23_test.x,sol23_test.y(1:2,:),...
    sol23_test.x,sol23_test.y(3:4,:),'LineWidth',2);
yline(0.5,'k--','LineWidth',2)
set(gca,set_bif{:})
title('C','FontSize',ft,'FontWeight','normal')
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','','interpreter','latex',FontSize=ft,FontWeight='bold')
xlabel('time (s)',FontSize=ft,FontWeight='normal')
yticks([0,0.5,1])
xticks([0.5,1,1.5,2])
xlim([0,2])
%%
load('branching_towrds_asymmetric_sols_original_case.mat')
%%
load('tracking_threshold_crossing_asymmetric_try2.mat')
%  rp_symbk1=arrayfun(@(x)x.parameter(in.PR),Symbk_org_br_with_stab.point);
%  df_symbk1=arrayfun(@(x)x.parameter(in.df),Symbk_org_br_with_stab.point);
% df_m=arrayfun(@(x)x.parameter(in.df),mbranch_df_wbifs.point);
% pr_m=arrayfun(@(x)x.parameter(in.PR),mbranch_df_wbifs.point);
%%
ft=14;
ft2=17;
set_bif={'LineWidth',2,'Box','on','FontSize',13,'FontWeight','normal'};

fiss=readtable('fission.csv'); fiss=fiss{:,:}; [~,idx]=sort(fiss(:,2)); fiss=fiss(idx,:);
coher=readtable('coherence.csv'); coher=coher{:,:}; [~,idx]=sort(coher(:,2)); coher=coher(idx,:);
figure(100)
clf;
tiledlayout(1,3,'TileSpacing','compact')
nexttile
hold on;%grid on
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
pltbif1=plot(pr1_m,df1_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
xline(20,'--','Color',0.4*[1 1 1],'LineWidth',2)
text(10.5, 0.65, 'valid region','Color', 0.4*[1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Rotation',0,'FontSize',14,'FontWeight','normal');

yticks([0,0.5,1])
ylim([0,1])
set(gca,set_bif{:})
lgtext={'case: $(t_\mathrm{d},D)=(0.05,0.015)$',' case: $ (t_\mathrm{d},D)=(0.022,0.015)$', 'exp data'};%,...
   % 'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1)];%;plt3(1)];
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=ft2,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=ft2,FontWeight='bold')
legend(lgvect,lgtext,...
    'interpreter','latex',FontSize=ft,FontWeight='normal')
title('A','FontSize',ft,'FontWeight','normal')
nexttile
hold on; 
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
pltbif1=plot(pr2_m,df2_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
xline(20,'--','Color',0.5*[1 1 1],'LineWidth',2)
text(10.5, 0.65, 'valid region','Color', 0.4*[1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Rotation',0,'FontSize',14,'FontWeight','normal');
yticks([0,0.5,1])
ylim([0,1])
set(gca,set_bif{:})
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=ft2,FontWeight='bold')
lgtext={'case: $(t_\mathrm{d},D)=(0.05,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$','exp data'};%...
   % 'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1)]%;plt3(1)];
legend(lgvect,lgtext,'interpreter','latex',FontSize=ft,FontWeight='normal')
title('B','FontSize',ft,'FontWeight','normal')
nexttile
hold on; 
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',8);
pltbif1=plot(pr3_m,df3_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
yticks([0,0.5,1])
ylim([0,1])
set(gca,set_bif{:})
title('C','FontSize',ft,'FontWeight','normal')
%legend('case: $(t_\mathrm{d},D)=(0.022,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$','interpreter','latex',FontSize=17,FontWeight='normal')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=ft2,FontWeight='normal')
lgtext={'case: $(t_\mathrm{d},D)=(0.022,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$','exp data'};%...
   % 'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1)];%plt3(1)];
legend(lgvect,lgtext,'interpreter','latex',FontSize=ft,FontWeight='normal')
%%
%%%%
ft=14;
set_bif={'LineWidth',2,'Box','on','FontSize',13,'FontWeight','normal'};

fiss=readtable('fission.csv'); fiss=fiss{:,:}; [~,idx]=sort(fiss(:,2)); fiss=fiss(idx,:);
coher=readtable('coherence.csv'); coher=coher{:,:}; [~,idx]=sort(coher(:,2)); coher=coher(idx,:);
figure(2000)
clf;
%tiledlayout(1,3)%,'TileSpacing','compact')
rr=1; cc=1;
%nexttile([rr cc])
hold on; %grid on
pltbif1=plot(pr3_m,df3_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',clrs(7,:),'LineWidth',2,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',clrs(5,:),'LineWidth',2,'MarkerSize',8);
yticks([0,0.5,1])
ylim([0,1])
%xlim([0,40])
set(gca,set_bif{:})
title('(d)','FontSize',ft,'FontWeight','normal','interpreter','latex')
%legend('case: $(t_\mathrm{d},D)=(0.022,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$','interpreter','latex',FontSize=17,FontWeight='normal')
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='normal')
lgtext={'case: $(t_\mathrm{d},D)=(0.022,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$',...
    'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1);plt3(1)];
legend(lgvect,lgtext,'interpreter','latex',FontSize=ft,FontWeight='normal')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')
set(figure(2000), 'Position', [10, 300, 700, 515])
figure(2001)
clf
hold on; %grid on
pltbif1=plot(pr2_m,df2_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
%xline(20,'--','Color',0.5*[1 1 1],'LineWidth',2)
%text(10.5, 0.65, 'valid region','Color', 0.4*[1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Rotation',0,'FontSize',14,'FontWeight','normal');
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',clrs(7,:),'LineWidth',2,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',clrs(5,:),'LineWidth',2,'MarkerSize',8);
yticks([0,0.5,1])
ylim([0,1])
%xlim([0,30])
set(gca,set_bif{:})
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='bold')
lgtext={'case: $(t_\mathrm{d},D)=(0.05,0.05)$','case: $(t_\mathrm{d},D)=(0.022,0.015)$',...
    'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1);plt3(1)];
legend(lgvect,lgtext,'interpreter','latex',FontSize=ft,FontWeight='normal')
title('(e)','FontSize',ft,'FontWeight','normal','interpreter','latex')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')

set(figure(2001), 'Position', [10, 300, 700, 515])
figure(2003)
clf
hold on
pltbif1=plot(pr1_m,df1_m,'k-','LineWidth',3);
pltbif2=plot(pr_m,df_m,'color',clrs(1,:),'LineWidth',3);
% xline(20,'--','Color',0.4*[1 1 1],'LineWidth',2)
% text(10.5, 0.65, 'valid region','Color', 0.4*[1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Rotation',0,'FontSize',14,'FontWeight','normal');
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color',clrs(7,:),'LineWidth',2,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'x','Color',clrs(5,:),'LineWidth',2,'MarkerSize',8);
yticks([0,0.5,1])
ylim([0,1])
%xlim([0,35])
set(gca,set_bif{:})
lgtext={'case: $(t_\mathrm{d},D)=(0.05,0.015)$',' case: $ (t_\mathrm{d},D)=(0.022,0.015)$',...
    'coherence curve','fission curve'};
lgvect=[pltbif1;pltbif2;plt2(1);plt3(1)];
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='bold')
%ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')
legend(lgvect,lgtext,...
    'interpreter','latex',FontSize=ft,FontWeight='normal')
% legend('case: $(t_\mathrm{d},D)=(0.05,0.015)$',' case: $ (t_\mathrm{d},D)=(0.022,0.015)$',...
%     'co','interpreter','latex',FontSize=17,FontWeight='normal')
title('(f)','FontSize',ft,'FontWeight','normal','interpreter','latex')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')

set(figure(2003), 'Position', [10, 300, 700, 515])
