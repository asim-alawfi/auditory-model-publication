%% Plotting one-parameter bifurcation (symmetric and non-symmetric POs) t_d, D cases
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
%%
load('asymmetric_branches_case12_try2.mat')
%%
set_one_bif={'FontWeight','bold','FontSize',12,'FontName','Aril'};
codprojmat=[1,0,0,0,0,0];
%%  Plotting for one-parameter branch for $(t_\mathrm{d }, D) = (0.05,0.015)$
Scond=dde_lincond_struct(size(asym_brs1_wbifs.point(1).profile,1),'profile',...
    'shift',[1,2],'condprojmat',codprojmat,'condprojint',[0,0.5],'trafo',Rsym);
%%
Sint_A1=dde_lincond_struct(size(asym_brs1_wbifs.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[1,0,0,0,0,0],'condprojint',[0,0.5]);
Sint_B2=dde_lincond_struct(size(asym_brs1_wbifs.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[0,1,0,0,0,0],'condprojint',[0.5,1]);
%
rp1_sx=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(1).point);
rp_asy1=arrayfun(@(x)x.parameter(in.PR),asym_brs1_wbifs.point);
%
yax1_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),br_symmetry_wbifs(1).point);
yax1_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),br_symmetry_wbifs(1).point);
non1_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),asym_brs1_wbifs.point);
non1_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),asym_brs1_wbifs.point);
%
max1_sy= yax1_Sint_A1-yax1_Sint_B2;%1arrayfun(@(x)dde_psol_lincond(x,Scond),br_symmetry_wbifs(1).point);
yaxis_asy1=non1_sym_Sint_A1-non1_sym_Sint_B2;%arrayfun(@(x)dde_psol_lincond(x,Scond),asym_brs1_wbifs.point);
%%
clrs2=lines();
figure(1)
clf;
hold on
plot(rp1_sx(nunst_po{1}'==0),max1_sy(nunst_po{1}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10),
plot(rp1_sx(nunst_po{1}'>=1),max1_sy(nunst_po{1}'>=1),'k--','LineWidth',3)%'MarkerSize',7)
plot(rp_asy1(nunst_1==0),yaxis_asy1(nunst_1==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10)
plot(rp_asy1(nunst_1>=1),yaxis_asy1(nunst_1>=1),'r.','MarkerSize',7)
plot(rp1_x(p_bif(:,1)),max1_sy(p_bif(:,1)),'k.','MarkerSize',25)
title(' $(t_\mathrm{d }, D) = (0.05,0.015)$','interpreter','latex','FontSize',16)
legend('stable symmetric POs','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif','')
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
set(gca,set_one_bif{:})
grid on
%% Plotting for one-parameter branch for $(t_\mathrm{d }, D) = (0.05,0.05)$
yax2_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),br_symmetry_wbifs(2).point);
yax2_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),br_symmetry_wbifs(2).point);
non2_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),asym_brs2_wbifs.point);
non2_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),asym_brs2_wbifs.point);
max2_sy= yax2_Sint_A1-yax2_Sint_B2;
yaxis_asy2=non2_sym_Sint_A1-non2_sym_Sint_B2;
rp2_sx=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(2).point);
rp_asy2=arrayfun(@(x)x.parameter(in.PR),asym_brs2_wbifs.point);
%
figure(2)
clf;
hold on
plot(rp2_sx(nunst_po{2}'==0),max2_sy(nunst_po{2}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10)
plot(rp2_sx(nunst_po{2}'>=1),max2_sy(nunst_po{2}'>=1),'k--','LineWidth',3)
plot(rp_asy2(nunst_td==0),yaxis_asy2(nunst_td'==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10)
    plot(rp_asy2(nunst_td>=1),yaxis_asy2(nunst_td>=1),'r.','MarkerSize',7)
plot(rp2_x(p_bif(:,2)),max2_sy(p_bif(:,2)),'k.','MarkerSize',25)
title('$(t_\mathrm{d }, D) = (0.05,0.05)$','interpreter','latex','FontSize',16)
legend('stable symmetric POs','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif','')
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
set(gca,set_one_bif{:})
grid on
%%  Plotting for one-parameter branch for $(t_\mathrm{d }, D) = (0.022,0.05)$
%load('asymm_brD.mat')
%%
rp3_Dx=arrayfun(@(x)x.parameter(in.PR),asym_brs3_wbifs.point);
rp3_sx=arrayfun(@(x)x.parameter(in.PR),br_symmetry_wbifs(3).point);
yax3_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),br_symmetry_wbifs(3).point);
yax3_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),br_symmetry_wbifs(3).point);
non3_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),asym_brs3_wbifs.point);
non3_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),asym_brs3_wbifs.point);
max3_sy= yax3_Sint_A1-yax3_Sint_B2;
max3_Dy=non3_sym_Sint_A1-non3_sym_Sint_B2;
%%
p_bifD=p_bif_3;
figure(3)
clf
hold on
plot(rp3_sx(nunst_po{3}'==0),max3_sy(nunst_po{3}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10)
plot(rp3_sx(nunst_po{3}'>=1),max3_sy(nunst_po{3}'>=1),'k--','LineWidth',3)
plot(rp3_Dx(nunst_3==0),max3_Dy(nunst_3==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10) %nunst_poD
plot(rp3_Dx(nunst_3>=1),max3_Dy(nunst_3>=1),'.','Color','[0.5 0.5 0.5]','LineWidth',3,'MarkerSize',8)
plot(rp3_x(p_bif(:,3)),max3_sy(p_bif(:,3)),'k.','MarkerSize',35);
plot(rp3_Dx(p_bifD(1)),max3_Dy(p_bifD(1)),'.','Color',clrs2(2,:),'MarkerSize',35) %p_bifD
plot(rp3_Dx(p_bifD(2)),max3_Dy(p_bifD(2)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(3)),max3_Dy(p_bifD(3)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(6)),max3_Dy(p_bifD(6)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(7)),max3_Dy(p_bifD(7)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(8)),max3_Dy(p_bifD(8)),'.','Color',clrs2(2,:),'MarkerSize',35)
legend('sym stable POs','sym untsable POs','non-sym stable POs','non-sym untsable POs',...
    'symmetry-breaking','fold of POs','period doubling bifurcation')
title(' $(t_\mathrm{d }, D) = (0.022,0.05)$','interpreter','latex','FontSize',16)
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
ylabel('$\int_0^{1/2}u_B(t)-u_A(t+1/2)\mathrm{d}t$','interpreter','latex','FontSize',16);
set(gca,set_one_bif{:})%set(gca,"FontWeight",'bold')
grid on
%%
figure(888)
clf
tiledlayout(1,3)
nexttile
hold on
plot(rp3_sx(nunst_po{3}'==0),max3_sy(nunst_po{3}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10)
plot(rp3_sx(nunst_po{3}'>=1),max3_sy(nunst_po{3}'>=1),'k--','LineWidth',3)
plot(rp3_Dx(nunst_3==0),max3_Dy(nunst_3==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10) %nunst_poD
plot(rp3_Dx(nunst_3>=1),max3_Dy(nunst_3>=1),'.','Color','[0.5 0.5 0.5]','LineWidth',3,'MarkerSize',8)
plot(rp3_x(p_bif(:,3)),max3_sy(p_bif(:,3)),'k.','MarkerSize',35);
plot(rp3_Dx(p_bifD(1)),max3_Dy(p_bifD(1)),'.','Color',clrs2(2,:),'MarkerSize',35) %p_bifD
plot(rp3_Dx(p_bifD(2)),max3_Dy(p_bifD(2)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(3)),max3_Dy(p_bifD(3)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(6)),max3_Dy(p_bifD(6)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(7)),max3_Dy(p_bifD(7)),'.','Color',clrs2(1,:),'MarkerSize',35)
plot(rp3_Dx(p_bifD(8)),max3_Dy(p_bifD(8)),'.','Color',clrs2(2,:),'MarkerSize',35)
legend('sym stable POs','sym untsable POs','non-sym stable POs','non-sym untsable POs',...
    'symmetry-breaking','period doubling bifurcation','fold of POs','FontSize',12)
title(' $(t_\mathrm{d }, D) = (0.022,0.05)$','interpreter','latex','FontSize',16)
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
ylabel('$\int_0^{1/2}u_B(t)-u_A(t+1/2)\mathrm{d}t$','interpreter','latex','FontSize',18);
set(gca,set_one_bif{:})%set(gca,"FontWeight",'bold')
grid on
nexttile
hold on
plot(rp2_sx(nunst_po{2}'==0),max2_sy(nunst_po{2}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10)
plot(rp2_sx(nunst_po{2}'>=1),max2_sy(nunst_po{2}'>=1),'k--','LineWidth',3)
plot(rp_asy2(nunst_td==0),yaxis_asy2(nunst_td'==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10)
    plot(rp_asy2(nunst_td>=1),yaxis_asy2(nunst_td>=1),'r.','MarkerSize',7)
plot(rp2_x(p_bif(:,2)),max2_sy(p_bif(:,2)),'k.','MarkerSize',35)
title('$(t_\mathrm{d }, D) = (0.05,0.05)$','interpreter','latex','FontSize',16)
legend('stable symmetric POs','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif','','FontSize',12)
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
set(gca,set_one_bif{:})
grid on
nexttile
hold on; grid on
plot(rp1_sx(nunst_po{1}'==0),max1_sy(nunst_po{1}'==0),'.','Color',clrs2(4,:),'LineWidth',3,'MarkerSize',10),
plot(rp1_sx(nunst_po{1}'>=1),max1_sy(nunst_po{1}'>=1),'k--','LineWidth',3)%'MarkerSize',7)
plot(rp_asy1(nunst_1==0),yaxis_asy1(nunst_1==0),'.','Color',clrs2(5,:),'LineWidth',3,'MarkerSize',10)
plot(rp_asy1(nunst_1>=1),yaxis_asy1(nunst_1>=1),'r.','MarkerSize',7)
plot(rp1_x(p_bif(:,1)),max1_sy(p_bif(:,1)),'k.','MarkerSize',35)
title('$(t_\mathrm{d }, D) = (0.05,0.015)$','interpreter','latex','FontSize',16)
legend('stable symmetric POs','unstable symmetric POs','stable non-symmetric POs','','symmetry-breaking Bif','FontSize',12)
xlabel('$r_\mathrm{p}$','FontName','Courier','Interpreter','latex','FontSize',12)
set(gca,set_one_bif{:})

