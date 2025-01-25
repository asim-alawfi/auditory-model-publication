
%% Plot the all results of the orginial case (one-parameter bifurcation, symmetry-breaking, touching threshold ...)
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
%load('branch_of_sympos_original_case_and_threshold_crossing_coj.mat')
load('branching_towrds_asymmetric_sols_original_case.mat')
load('sym_breaking_original_case_try2.mat')
load('identifying_solutions_in_regions_try2.mat')
%load('tracking_threshold_crossing_asymmetric_try2.mat') 
%%
graycolor=[0.7 0.7 0.7];
mbrown= [0.65, 0.33, 0.13];

%% plot one-parameter bifurcation in PR and df=0.73 symmetric and asymmetric solutions
rp2_x=arrayfun(@(x)x.parameter(in.PR),po2_symmetry.point);
max2_y=arrayfun(@(x)max(x.profile(1,:)),po2_symmetry.point);
min2_y2=arrayfun(@(x)min(x.profile(1,:)),po2_symmetry.point);
rp_unsym=arrayfun(@(x)x.parameter(in.PR),nonsymper_wbifs.point);
y_unsym=arrayfun(@(x)max(x.profile(1,:)),nonsymper_wbifs.point);
y2_unsym=arrayfun(@(x)min(x.profile(1,:)),nonsymper_wbifs.point);
[~,it_intg]=min(abs(rp2_x-3));
[~,it_seg]=min(abs(rp2_x-8.15));
it_bis=find(diff(sign(rp_unsym-4.0)));
%%
p_intg=po2_symmetry.point(it_intg);
p_seg=po2_symmetry.point(it_seg);
p_bis=nonsymper_wbifs.point(it_bis(1));
p_bis2=nonsymper_wbifs.point(it_bis(2));
%% We plot int_{0}^{1/2} u_A(t)-u_B(t+1/2) dt by setting:
Sint_A1=dde_lincond_struct(size(po2_symmetry.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[1,0,0,0,0,0],'condprojint',[0,0.5]);
Sint_B2=dde_lincond_struct(size(po2_symmetry.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[0,1,0,0,0,0],'condprojint',[0.5,1]);
yax_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),po2_symmetry.point);
yax_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),po2_symmetry.point);
non_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),nonsymper_wbifs.point);
non_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),nonsymper_wbifs.point);
yax_sym=yax_Sint_A1-yax_Sint_B2;
yax_nonsym=non_sym_Sint_A1-non_sym_Sint_B2;
%%
bif_loc=find(diff(nunst2_sym));
grayColor=[0.4 0.4 0.4];
pColor=[0 0.5 0];%[0.3010 0.7450 0.9330]%[0.9290 0.6940 0.1250]%[0.4940 0.1840 0.5560];
clr=colorcube();%lines()
p2c=[0.4940 0.1840 0.5560];
%%
clrs2=lines();
figure(1)
clf;
tiledlayout(4,4)
nexttile([4 2])
hold on; grid on
plot(rp2_x(1:bif_loc(1)),yax_sym(1:bif_loc(1)),'-','Color',clrs2(4,:),'LineWidth',7)
plot(rp2_x(bif_loc(2):end),yax_sym(bif_loc(end):end),'-','Color',clrs2(4,:),'LineWidth',7)
plot(rp2_x(nunst2_sym>=1),yax_sym(nunst2_sym>=1),'k--','LineWidth',4)%,'Color',grayColor)%,'LineWidth',1)
plot(rp_unsym(uns_um==0),yax_nonsym(uns_um==0),'.-','Color',clrs2(5,:),'LineWidth',7)%,'MarkerSize',10)%,'g.')%,'MarkerS',5)
plot(rp2_x(bif_loc),yax_sym(bif_loc),'k.','MarkerSize',55)
plot(rp2_x(it_intg+8),yax_sym(it_intg+8),'.','Color','red','MarkerSize',55)
plot(rp2_x(it_seg-8),yax_sym(it_seg-8),'.','Color',mbrown,'MarkerS',55)
plot(rp_unsym(it_bis),yax_nonsym(it_bis),'x','Color','red','LineWidth',5,'MarkerS',25)
%yline(0.5,'k--','Threshold')
legend('stable symmetric POs','','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif',...
    'solution in a1','solution in a2','solution in a3 and a4','interpreter','latex','FontName','Cambria',FontSize=32)%,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Cambria',FontSize=26,FontWeight='bold')
ylabel('$\int_0^{1/2}u_A(t)-u_B(t+1/2)\mathrm{d}t$','interpreter','latex','FontName','Cambria',FontSize=26,FontWeight='bold');
xlim([0,rp2_x(1)])
title('(A)','FontSize',26,'FontName','Cambria')
%ylabel('Second peak of $u_\mathrm{A}$ per period','FontName','Courier','Interpreter','latex','FontSize',12)
set_one_bif={'LineWidth',2,'Box','on'};
set(gca,set_one_bif{:})
xlim([1,30])

% %%%%%%%%%%%%%%%%%%%%%
nexttile([2 1])
hold on
plot(p_intg.mesh*p_intg.period,p_intg.profile(1:2,:),'LineWidth',7)
plot(0.6,0.9,'.','Color','red','MarkerSize',90)
yline(0.5,'k--','Threshold','LineWidth',5)
ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Cambria',FontSize=20,FontWeight='bold')
title('(a1)','FontSize',26,'FontName','Cambria')
set(gca,set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_intg.period])
%%%%%%%%%%%%%%%%%%
nexttile([2 1])
hold on 
plot(p_seg.mesh*p_seg.period,p_seg.profile(1:2,:),'LineWidth',7)
plot(0.2,0.9,'.','Color',mbrown,'MarkerSize',90)
yline(0.5,'k--','Threshold','LineWidth',5)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Cambria',FontSize=28,FontWeight='bold')
title('(a2)','FontSize',26,'FontName','Cambria')
set(gca, set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_seg.period])
%
nexttile([2 1])
hold on 
plot(p_bis.mesh*p_bis.period,p_bis.profile(1:2,:),'LineWidth',7)
plot(0.4,0.9,'x','Color','red','MarkerSize',50,'LineWidth',7)

yline(0.5,'k--','Threshold','LineWidth',5)
ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Cambria',FontSize=20,FontWeight='bold')
xlabel('period','FontName','Cambria',FontSize=26,FontWeight='bold')
title('(a3)','FontSize',26,'FontName','Cambria')
set(gca,  set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_bis.period])
nexttile([2 1])
hold on 
plot(p_bis2.mesh*p_bis2.period,p_bis2.profile(1:2,:),'LineWidth',7)
plot(0.4,0.9,'x','Color','red','MarkerSize',50,'LineWidth',7)
yline(0.5,'k--','Threshold','LineWidth',5)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','','interpreter','latex','FontName','Courier',FontSize=12,FontWeight='bold')
xlabel('period','FontName','Cambria',FontSize=26,FontWeight='bold')
title('(a4)','FontSize',26,'FontName','Cambria')
set(gca,  set_one_bif{:})
grid on
ylim([0,1])%axis tight
xlim([0,p_bis2.period])
%%
[~,it]=min(abs(rp_thta-34));
pt=mbr_wbifs.point(it);
t5=find(diff(sign(pr_m-6.658)));
pt5=mbranch_df.point(t5(1));
figure(2)
clf;
tiledlayout(4,4)
nexttile([4 2])
hold on; grid on
plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b.-',...
rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plot(rp_thta(dum==0),df_thta(dum==0),'m.-','MarkerSize',8)
plot(rp_thta(dum>=1),df_thta(dum>=1),'.-','color',graycolor,'MarkerSize',8)
plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'k--','LineWidth',1.2)
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',12,'LineWidth',2)
plot(pr_m(t5(1)),df_m(t5(1)),'kx','MarkerSize',12,'LineWidth',2)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'rx',...
non_symper.point(it4).parameter(in.PR),non_symper.point(it4).parameter(in.df),'gx','MarkerSize',12,'LineWidth',2);
legend(gca,'symmetry-breaking','touching \theta: stable symmetric POs','touching \theta: unstable symmetric POs',...
'touching \theta:stable non-symmetric POs','', 'solution in b1','solution in b2','solution in c1','solution in c2',...
'','FontSize',12,'FontName','Cambria')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Cambria',FontSize=16,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Cambria',FontSize=16,FontWeight='bold')
title('(a)','FontSize',12,'FontName','Cambria')
ylim([0,1] )
set(gca, 'FontWeight','bold')
nexttile([2 1])
plot(pt.mesh*pt.period,pt.profile(1:2,:),'LineWidth',2);
yline(0.5,'k--','Threshold','LineWidth',3)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
ylabel('Neural activity','FontName','Cambria',FontSize=12)
%xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b1)','FontSize',12,'FontName','Cambria')
xlim([0,pt.period])
ylim([0,0.8])
set(gca, 'FontWeight','bold')
grid on
nexttile([2 1])
plot(pt5.mesh*pt5.period,pt5.profile(1:2,:),'LineWidth',2);
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
%xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b2)','FontSize',12,'FontName','Cambria')
xlim([0,pt5.period])
ylim([0,0.8])
set(gca, 'FontWeight','bold')
grid on
nexttile([2 1])
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
ylabel('Neural activity','FontName','Cambria',FontSize=12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(c1)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
xlim([0,pt3.period])
ylim([0,0.8])
grid on
nexttile([2 1])
plot(pt4.mesh*pt4.period,pt4.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(c2)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
xlim([0,pt4.period])
ylim([0,0.8])
grid on
%% Ploting the maximum of the second peak of u_A over the interval I_{B}
c_A2=zeros(1,12);
c_A2(1)=1;
ua_eval_asymbk=@(p)max(c_A2*dde_coll_eva(p.profile,p.mesh,linspace(0.5,1,1e4),p.degree));
second_peak=zeros(1,length(Symbk_org_br_with_stab.point));
for i=1:length(Symbk_org_br_with_stab.point)
    pi=Symbk_org_br_with_stab.point(i);
    second_peak(i)=ua_eval_asymbk(pi);
end
figure(3)
clf; hold on
plot(rp_symbk1,second_peak,'b.','MarkerSize',15)
% clrs=colormap(jet);
% tmp=@(x)round(x*255+1);
% mscal=tmp((second_peak-min(second_peak))/(max(second_peak)-min(second_peak)));
% figure(3)
% clf; hold on
% for i=1:length(second_peak)
% plot(rp_symbk1(i),second_peak(i),'.','Color',clrs(mscal(i),:),'MarkerSize',15)
% end
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Cambria',FontSize=16,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Cambria',FontSize=16,FontWeight='bold')
set(gca, 'FontWeight','bold')
grid on
%%
p_sym=mbr_wbifs.point(it);
inp=auditory_forcing(mbr_wbifs,sig,it,in);
figure(4)
clf;
subplot(2,1,1)
hold on
plot(p_sym.mesh*p_sym.period,p_sym.profile(1:2,:),'LineWidth',2)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','','Interpreter','latex','FontSize',12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
grid on
ylim([0,1])%axis tight
xlim([0,p_sym.period])
grid on
hold on
subplot(2,1,2)
hold on
plot(p_sym.mesh*p_sym.period,inp(1,:),'b',p_sym.mesh*p_sym.period,inp(2,:),'r','LineWidth',2)
grid on
legend('$i_\mathrm{A}(t)$','$i_\mathrm{B}(t)$','Interpreter','latex','FontSize',12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
ylim([0,8])
xlim([0,p_sym.period])
%%
my_branch=mbr_wbifs;
point_ind=30;
point=mbr_wbifs.point(point_ind);
forces=auditory_forcing(mbr_wbifs,sig,point_ind,in);
figure(5)
clf;
subplot(2,1,1);
hold on
plot(point.mesh*point.period,point.profile(1:2,:),'LineWidth',2)
yline(0.5,'k--','LineWidth',2)
subplot(2,1,2)
hold on
plot(point.mesh*point.period,forces(1,:),'b.',...
    point.mesh*point.period,forces(2,:),'r.','LineWidth',2)
title('PR=',point.parameter(in.PR))
%%
figure(111)
clf;
tiledlayout(12,8)
nexttile([6 4])
hold on; grid on
plot(rp2_x(1:bif_loc(1)),yax_sym(1:bif_loc(1)),'-','Color',clrs2(4,:),'LineWidth',7)
plot(rp2_x(bif_loc(2):end),yax_sym(bif_loc(end):end),'-','Color',clrs2(4,:),'LineWidth',7)
plot(rp2_x(nunst2_sym>=1),yax_sym(nunst2_sym>=1),'k--','LineWidth',4)%,'Color',grayColor)%,'LineWidth',1)
plot(rp_unsym(uns_um==0),yax_nonsym(uns_um==0),'.-','Color',clrs2(5,:),'LineWidth',7)%,'MarkerSize',10)%,'g.')%,'MarkerS',5)
plot(rp2_x(bif_loc),yax_sym(bif_loc),'k.','MarkerSize',55)
plot(rp2_x(it_intg+8),yax_sym(it_intg+8),'.','Color','red','MarkerSize',55)
plot(rp2_x(it_seg-8),yax_sym(it_seg-8),'.','Color',mbrown,'MarkerS',55)
plot(rp_unsym(it_bis),yax_nonsym(it_bis),'x','Color','red','LineWidth',5,'MarkerS',25)
legend('stable symmetric POs','','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif',...
    'solution in a1','solution in a2','solution in a3 and a4','interpreter','latex','FontName','Courier',FontSize=26,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Courier',FontSize=22,FontWeight='bold')
ylabel('$\int_0^{1/2}u_A(t)-u_B(t+1/2)\mathrm{d}t$','interpreter','latex','FontName','Courier',FontSize=22,FontWeight='bold');
xlim([0,rp2_x(1)])
title('(A)','FontSize',22,'FontName','Courier')
set_one_bif={'LineWidth',2,'Box','on'};
set(gca,set_one_bif{:})
xlim([1,30])
% %%%%%%%%%%%%%%%%%%%%%
nexttile([6 4])
hold on; grid on
plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b-','color',clrs2(1,:),'LineWidth',4)
%plot(rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plot(rp_thta(dum==0),df_thta(dum==0),'m-','Linewidth',3)
plot(rp_thta(dum>=1),df_thta(dum>=1),'','color',graycolor,'Linewidth',3)
plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'--','Color','k','LineWidth',3)
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',20,'LineWidth',4)
plot(pr_m(t5(1)),df_m(t5(1)),'kx','MarkerSize',20,'LineWidth',4)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'rx',...
non_symper.point(it4).parameter(in.PR),non_symper.point(it4).parameter(in.df),'gx','MarkerSize',20,'LineWidth',4);
legend(gca,'symmetry-breaking','touching $\theta$: stable symmetric POs','touching $\theta$: unstable symmetric POs',...
'touching $\theta$: stable non-symmetric POs','', 'solution in b1','solution in b2','solution in b3','solution in b4',...
'','interpreter','latex','FontName','Courier','FontSize',18,'FontWeight','bold')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Courier',FontSize=22,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Courier',FontSize=22,FontWeight='bold')
title('(B)','FontSize',22,'FontName','Courier')
ylim([0,1] )
set(gca, set_one_bif{:})
%%%%%
nexttile([3 2])
plot(p_intg.mesh*p_intg.period,p_intg.profile(1:2,:),'LineWidth',5)
yline(0.5,'k--','Threshold','LineWidth',4)
%ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=20,FontWeight='bold')
title('(a1)','FontSize',22,'FontName','Courier')
set(gca,set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_intg.period])
xticks([p_intg.period/2,p_intg.period])
%
nexttile([3 2])
hold on 
plot(p_seg.mesh*p_seg.period,p_seg.profile(1:2,:),'LineWidth',7)
yline(0.5,'k--','Threshold','LineWidth',4)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Courier',FontSize=26,FontWeight='bold')
title('(a2)','FontSize',22,'FontName','Courier')
set(gca, set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_seg.period])
xticks([p_seg.period/2,p_seg.period])
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile([3 2])
plot(pt.mesh*pt.period,pt.profile(1:2,:),'LineWidth',7);
yline(0.5,'k--','Threshold','LineWidth',5)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
title('(b1)','FontSize',22,'FontName','Courier')
xlim([0,pt.period])
ylim([0,1])
xticks([pt.period/2,pt.period])
set(gca, 'FontWeight','bold')
grid on
nexttile([3 2])
plot(pt5.mesh*pt5.period,pt5.profile(1:2,:),'LineWidth',5);
yline(0.5,'k--','Threshold','LineWidth',4)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
%xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b2)','FontSize',22,'FontName','Courier')
xlim([0,pt5.period])
xticks([pt5.period/2,pt5.period])

ylim([0,1])
set(gca, set_one_bif{:})
grid on


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile([3 2])
hold on 
plot(p_bis.mesh*p_bis.period,p_bis.profile(1:2,:),'LineWidth',5)
yline(0.5,'k--','Threshold','LineWidth',4)
%ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=20,FontWeight='bold')
xlabel('time(s)','interpreter','latex','FontName','Courier',FontSize=20,FontWeight='bold')
title('(a3)','FontSize',22,'FontName','Courier')
set(gca,  set_one_bif{:})
grid on
ylim([0,1])
xlim([0,p_bis.period])
xticks([p_bis.period/2,p_bis.period])

nexttile([3 2])
hold on 
plot(p_bis2.mesh*p_bis2.period,p_bis2.profile(1:2,:),'LineWidth',5)
yline(0.5,'k--','Threshold','LineWidth',4)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','','interpreter','latex','FontName','Courier',FontSize=12,FontWeight='bold')
xlabel('time(s)','interpreter','latex','FontName','Courier',FontSize=20,FontWeight='bold')
title('(a4)','FontSize',22,'FontName','Courier')
set(gca, set_one_bif{:})
grid on
ylim([0,1])%axis tight
xlim([0,p_bis2.period])
xticks([p_bis2.period/2,p_bis2.period])

nexttile([3 2])
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',5)
yline(0.5,'k--','Threshold','LineWidth',4)
%legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b3)','FontSize',22,'FontName','Courier')
set(gca,  set_one_bif{:})
xlim([0,pt3.period])
xticks([pt3.period/2,pt3.period])

ylim([0,1])
grid on
nexttile([3 2])
plot(pt4.mesh*pt4.period,pt4.profile(1:2,:),'LineWidth',5)
yline(0.5,'k--','Threshold','LineWidth',4)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',12)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
xlabel('time(s)','FontSize',12,'FontName','Cambria')
title('(b4)','FontSize',22,'FontName','Courier')
set(gca, set_one_bif{:})
xlim([0,pt4.period])
xticks([pt4.period/2,pt4.period])

ylim([0,1])
grid on
% %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%
figure(114)
clf;
tiledlayout(4,4)
nexttile([4 2])
% %%%%%%%%%%%%%%%%%%%%%
hold on; grid on
plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b-','color',clrs2(1,:),'LineWidth',7)
%plot(rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plot(rp_thta(dum==0),df_thta(dum==0),'m-','Linewidth',5)
plot(rp_thta(dum>=1),df_thta(dum>=1),'','color',graycolor,'Linewidth',4)
plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'--','Color','k','LineWidth',5)
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',20,'LineWidth',4)
plot(pr_m(t5(1)),df_m(t5(1)),'x','Color',clrs2(3,:),'MarkerSize',20,'LineWidth',4)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'rx',...
    'Color',clrs2(7,:),'MarkerSize',20,'LineWidth',4)
plot(non_symper.point(it4).parameter(in.PR),non_symper.point(it4).parameter(in.df),'gx','MarkerSize',20,'LineWidth',4);
legend(gca,'symmetry-breaking','touching $\theta$: stable symmetric POs','touching $\theta$: unstable symmetric POs',...
'touching $\theta$: stable non-symmetric POs','', 'solution in b1','solution in b2','solution in b3','solution in b4',...
'','interpreter','latex','FontName','Cambria','FontSize',24)
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Cambria',FontSize=28,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex','FontName','Cambria',FontSize=28,FontWeight='bold')
title('(B)','FontSize',26,'FontName','Cambria')
ylim([0,1] )
set(gca, set_one_bif{:})
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile([2 1])
hold on
plot(pt.mesh*pt.period,pt.profile(1:2,:),'LineWidth',7);
yline(0.5,'k--','Threshold','LineWidth',5)
plot(0.05,0.9,'x','Color','m','MarkerSize',50,'LineWidth',7)
title('(b1)','FontSize',26,'FontName','Cambria')
xlim([0,pt.period])
ylim([0,1])
set(gca,  set_one_bif{:})
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile([2 1])
hold on
plot(pt5.mesh*pt5.period,pt5.profile(1:2,:),'LineWidth',7);
yline(0.5,'k--','Threshold','LineWidth',5)
plot(0.25,0.9,'x','Color',clrs2(3,:),'MarkerSize',50,'LineWidth',7)
title('(b2)','FontSize',26,'FontName','Cambria')
xlim([0,pt5.period])
%xticks([pt5.period/2,pt5.period])
ylim([0,1])
set(gca, set_one_bif{:})
grid on
nexttile([2 1])
hold on
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',7)
plot(0.275,0.9,'x','Color',clrs2(7,:),'MarkerSize',50,'LineWidth',7)
yline(0.5,'k--','Threshold','LineWidth',5)
xlabel('period','FontName','Cambria',FontSize=26,FontWeight='bold')
title('(b3)','FontSize',26,'FontName','Cambria')
set(gca,  set_one_bif{:})
xlim([0,pt3.period])
ylim([0,1])
grid on
nexttile([2 1])
hold on
plot(pt4.mesh*pt4.period,pt4.profile(1:2,:),'LineWidth',7)
plot(0.2,0.9,'x','Color','g','MarkerSize',50,'LineWidth',7)
yline(0.5,'k--','Threshold','LineWidth',5)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',30)
%ylabel('Neural activity','FontName','Cambria',FontSize=12)
xlabel('period','FontName','Cambria',FontSize=26,FontWeight='bold')
title('(b4)','FontSize',26,'FontName','Cambria')
set(gca, set_one_bif{:})
xlim([0,pt4.period])

ylim([0,1])
grid on
%% 
load('audi_simulations.mat')
%%
grid on
figure(444)
clf
tiledlayout(2,7,"TileSpacing","compact")%(12,8)
nexttile([2,3])
hold on
fiss=readtable('fission.csv'); fiss=fiss{:,:}; [~,idx]=sort(fiss(:,2)); fiss=fiss(idx,:);
coher=readtable('coherence.csv'); coher=coher{:,:}; [~,idx]=sort(coher(:,2)); coher=coher(idx,:);
plt2=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'--','Color',[0 0 0],'LineWidth',3,'MarkerSize',8);
plt3=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'--','Color',[0 0 0],'LineWidth',3,'MarkerSize',8);
x_text_fiss = 16.7; 
y_text_fiss = 0.2; 
text(x_text_fiss, y_text_fiss, 'fission', 'FontSize', 20, 'Color', 'black','Rotation', 0)
x_text_coh = 14.2; % X-coordinate for the text
y_text_coh= 0.3; 
text(x_text_coh,y_text_coh,'coherence', 'FontSize', 20, 'Color', 'black','Rotation',0)
plot(par(in.PR),par(in.df),'.','Color',clrs2(5,:),'MarkerSize',50)
plot(par(in.PR),par2(in.df),'.','Color',clrs2(6,:),'MarkerSize',50)
plot(par(in.PR),par1(in.df),'.','Color',clrs2(7,:),'MarkerSize',50)
text(10.3,0.1 ,'integration (1)', 'FontSize', 28, 'Color', clrs2(5,:),'Rotation',0,'FontWeight','bold')
text(10,0.23,'bistability (2)', 'FontSize', 28, 'Color', clrs2(6,:),'Rotation',0,'FontWeight','bold')
text(10.3,par1(in.df),'segregation (3)', 'FontSize', 28, 'Color', clrs2(7,:),'Rotation',0,'FontWeight','bold')
% plt_legend={sprintf('Mathematical approx-\nmation formula'),sprintf('experimental data \n(fission & coherence)')};
plt_legend={sprintf('experimental data\n(fission & coherence)')};
plt_vec=[plt2(1)];
legend(plt_vec,plt_legend,'FontSize',34,'Location','best','FontName','Cambria','EdgeColor',0.7*[1,1,1])
set(gca,'Box','on','LineWidth',2)%'PlotBoxAspectRatio',[1.6,1,1],'LineWidth',2)
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=28,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=28,FontWeight='bold')
title('(a)','FontSize',24,'FontName','Cambria')
ylim([0,1])
xlim([7,20])
grid on
nexttile([1,2])
hold on
grid on;
plot(sol23.x,sol23.y(1:2,:),'LineWidth',6);
plot(sol23.x,sol23.y(3:4,:),'-','LineWidth',3);
yline(0.5,'k--','Threshold','LineWidth',2,'FontWeight','bold')
plot(3.875,0.75,'.','Color',clrs2(5,:),'MarkerSize',90)
text(3.72,0.725 ,'$\mathrm{I_{A}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
text(3.82,0.725 ,'$\mathrm{I_{B}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
xlim([3.7,3.9])
ylim([0,0.8])
set(gca,  set_one_bif{:},'XTickLabel', [], 'YTickLabel', [])
title('(b)','FontSize',24,'FontName','Cambria')
xlabel('period','FontName','Cambria',FontSize=20,FontWeight='bold')
nexttile([1,2])
hold on;
grid on
plot(sol23_df2.x,sol23_df2.y(1:2,:),'LineWidth',6);
plot(sol23_df2.x,sol23_df2.y(3:4,:),'LineWidth',3);
yline(0.5,'k--','Threshold','LineWidth',2,'FontWeight','bold')
plot(3.875,0.75,'.','Color',clrs2(6,:),'MarkerSize',90)
text(3.72,0.725 ,'$\mathrm{I_{A}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
text(3.82,0.725 ,'$\mathrm{I_{B}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
ylim([0,0.8])
set(gca,  set_one_bif{:},'XTickLabel', [], 'YTickLabel', [])
xlim([3.7,3.9])
title('(c)','FontSize',24,'FontName','Cambria')
xlabel('period','FontName','Cambria',FontSize=20,FontWeight='bold')
nexttile([1,2])
hold on
grid on
plot(sol23_df.x,sol23_df.y(1:2,:),'LineWidth',6);
plot(sol23_df.x,sol23_df.y(3:4,:),'LineWidth',3);
text(3.72,0.725 ,'$\mathrm{I_{A}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
text(3.82,0.725 ,'$\mathrm{I_{B}}$', 'FontSize', 20, 'Color', 'k','interpreter','latex','FontWeight','bold')
yline(0.5,'k--','Threshold','LineWidth',2,'FontWeight','bold')
plot(3.875,0.75,'.','Color',clrs2(7,:),'MarkerSize',90)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','interpreter','latex','FontName', ...
    'Cambria',FontSize=40,FontWeight='bold')
xlim([3.7,3.9])
ylim([0,0.8])
set(gca,  set_one_bif{:},'XTickLabel', [], 'YTickLabel', [])
title('(d)','FontSize',24,'FontName','Cambria')
xlabel('period','FontName','Cambria',FontSize=20,FontWeight='bold')
