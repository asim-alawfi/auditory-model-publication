
%% Plot the all results of the orginial case (one-parameter bifurcation, symmetry-breaking, touching threshold ...)
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
load('branch_of_sympos_original_case_and_threshold_crossing.mat')
load('branching_towrds_asymmetric_sols_original_case.mat')
load('sym_breaking_original_case_try2.mat')
load('identifying_solutions_in_regions_try2.mat')
load('tracking_threshold_crossing_asymmetric_try2.mat')
%%
graycolor=[0.7 0.7 0.7];
mbrown= [0.65, 0.33, 0.13];
set_one_bif={'FontWeight','bold','FontSize',12,'FontName','Aril'};
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
Ssym=dde_lincond_struct(size(po2_symmetry.point(1).profile,1),'profile',...
    'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'condprojint',[0,0.5]);
yax_sym=arrayfun(@(x)dde_psol_lincond(x,Ssym),po2_symmetry.point);

Snonsym=dde_lincond_struct(size(nonsymper_wbifs.point(1).profile,1),'profile',...
    'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'condprojint',[0,0.5]);
yax_nonsym=arrayfun(@(x)dde_psol_lincond(x,Snonsym),nonsymper_wbifs.point);
%%
bif_loc=find(diff(nunst2_sym));
grayColor=[0.4 0.4 0.4];
pColor=[0 0.5 0];%[0.3010 0.7450 0.9330]%[0.9290 0.6940 0.1250]%[0.4940 0.1840 0.5560];
clr=colorcube();%lines()
p2c=[0.4940 0.1840 0.5560];

%%
%%
figure(1)
clf;
tiledlayout(4,4)
nexttile([4 2])
hold on; grid on
plot(rp2_x(nunst2_sym==0),yax_sym(nunst2_sym==0),'b.')
plot(rp2_x(nunst2_sym>=1),yax_sym(nunst2_sym>=1),'k--','LineWidth',2)%,'Color',grayColor)%,'LineWidth',1)
plot(rp_unsym(uns_um==0),yax_nonsym(uns_um==0),'g.')%,'MarkerS',5)
plot(rp2_x(bif_loc),yax_sym(bif_loc),'k.','MarkerSize',35)
plot(rp2_x(it_intg+8),yax_sym(it_intg+8),'.','Color','red','MarkerSize',35)
plot(rp2_x(it_seg-8),yax_sym(it_seg-8),'.','Color',mbrown,'MarkerS',35)
plot(rp_unsym(it_bis),yax_nonsym(it_bis),'x','Color','red','LineWidth',3,'MarkerS',10)
%yline(0.5,'k--','Threshold')
legend('stable symmetric POs','unstable symmetric POs','stable non-symmetric POs','symmetry-breaking Bif',...
    'solution in b1','solution in b2','solution in c1 and c2','interpreter','latex','FontName','Courier',FontSize=14,FontWeight='bold')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
ylabel('$\int_0^{1/2}u_A(t)-u_B(t+1/2)\mathrm{d}t$','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold');
xlim([0,rp2_x(1)])
title('(a)','FontSize',12,'FontName','Cambria')
%ylabel('Second peak of $u_\mathrm{A}$ per period','FontName','Courier','Interpreter','latex','FontSize',12)
set(gca,set_one_bif{:})
xlim([1,30])

% %%%%%%%%%%%%%%%%%%%%%
nexttile([2 1])
plot(p_intg.mesh*p_intg.period,p_intg.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
title('(b1)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
grid on
ylim([0,1])
xlim([0,p_intg.period])
%
nexttile([2 1])
hold on 
plot(p_seg.mesh*p_seg.period,p_seg.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex','FontName','Courier',FontSize=14,FontWeight='bold')
title('(b2)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
grid on
ylim([0,1])
xlim([0,p_seg.period])
%
nexttile([2 1])
hold on 
plot(p_bis.mesh*p_bis.period,p_bis.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
ylabel('$u_\mathrm{A}$ and $u_\mathrm{B}$ activity','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
xlabel('time(s)','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
title('(c1)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
grid on
ylim([0,1])
xlim([0,p_bis.period])
nexttile([2 1])
hold on 
plot(p_bis2.mesh*p_bis2.period,p_bis2.profile(1:2,:),'LineWidth',3)
yline(0.5,'k--','Threshold','LineWidth',3)
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','','interpreter','latex','FontName','Courier',FontSize=12,FontWeight='bold')
xlabel('time(s)','interpreter','latex','FontName','Courier',FontSize=16,FontWeight='bold')
title('(c2)','FontSize',12,'FontName','Cambria')
set(gca, 'FontWeight','bold')
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