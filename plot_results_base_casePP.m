
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
load('branching_towrds_asymmetric_sols_original_case.mat')
load('sym_breaking_original_case_try2.mat')
load('identifying_solutions_in_regions_try2.mat')
load('PR_large_df_small_part2.mat')


%%
graycolor=[0.7 0.7 0.7];
mbrown= [0.65, 0.33, 0.13];

%% plot one-parameter bifurcation in PR and df=0.73 symmetric and asymmetric solutions
rp2_x=arrayfun(@(x)x.parameter(in.PR),po2_symmetry_wbifs.point);
max2_y=arrayfun(@(x)max(x.profile(1,:)),po2_symmetry_wbifs.point);
min2_y2=arrayfun(@(x)min(x.profile(1,:)),po2_symmetry_wbifs.point);
rp_unsym=arrayfun(@(x)x.parameter(in.PR),nonsymper_wbifs.point);
y_unsym=arrayfun(@(x)max(x.profile(1,:)),nonsymper_wbifs.point);
y2_unsym=arrayfun(@(x)min(x.profile(1,:)),nonsymper_wbifs.point);
[~,it_intg]=min(abs(rp2_x-3));
[~,it_seg]=min(abs(rp2_x-8.15));
it_bis=find(diff(sign(rp_unsym-4.0)));
%%
p_intg=po2_symmetry_wbifs.point(it_intg);
p_seg=po2_symmetry_wbifs.point(it_seg);
p_bis=nonsymper_wbifs.point(it_bis(1));
p_bis2=nonsymper_wbifs.point(it_bis(2));
%% We plot int_{0}^{1/2} u_A(t)-u_B(t+1/2) dt by setting:
Sint_A1=dde_lincond_struct(size(po2_symmetry.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[1,0,0,0,0,0],'condprojint',[0,0.5]);
Sint_B2=dde_lincond_struct(size(po2_symmetry.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,2],'condprojmat',-1,'stateproj',[0,1,0,0,0,0],'condprojint',[0.5,1]);
yax_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),po2_symmetry_wbifs.point);
yax_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),po2_symmetry_wbifs.point);
non_sym_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A1),nonsymper_wbifs.point);
non_sym_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B2),nonsymper_wbifs.point);
yax_sym=yax_Sint_A1-yax_Sint_B2;
yax_nonsym=non_sym_Sint_A1-non_sym_Sint_B2;
%%
bif_loc=find(diff(nunst3_sym));
grayColor=[0.4 0.4 0.4];
pColor=[0 0.5 0];%[0.3010 0.7450 0.9330]%[0.9290 0.6940 0.1250]%[0.4940 0.1840 0.5560];
clr=colorcube();%lines()
p2c=[0.4940 0.1840 0.5560];
%%
ft=13;
clrs2=lines();
figure(600)
clf;
tiledlayout(4,2,'TileSpacing','compact')
nexttile([2,2])
hold on
plot(rp2_x(1:bif_loc(2)),yax_sym(1:bif_loc(2)),'-','Color',clrs2(4,:),'LineWidth',3)
plot(rp2_x(bif_loc(3):end),yax_sym(bif_loc(3):end),'-','Color',clrs2(4,:),'LineWidth',3)
plot(rp2_x(nunst3_sym>=1),yax_sym(nunst3_sym>=1),'k--','LineWidth',3)%,'Color',grayColor)%,'LineWidth',1)
plot(rp_unsym(uns_um==0),yax_nonsym(uns_um==0),'.-','Color',clrs2(5,:),'LineWidth',3)%,'MarkerSize',10)%,'g.')%,'MarkerS',5)
plot(rp2_x(bif_loc(2:end)),yax_sym(bif_loc(2:end)),'k.','MarkerSize',40)
plot(rp2_x(it_intg-12),yax_sym(it_intg-12),'.','Color','red','MarkerSize',40)
plot(rp2_x(it_seg+8),yax_sym(it_seg+8),'.','Color',mbrown,'MarkerS',40)
plot(rp_unsym(it_bis),yax_nonsym(it_bis),'x','Color','red','LineWidth',3,'MarkerSize',15)
set_one_bif={'LineWidth',2,'Box','on','FontSize',10,'FontWeight','normal'};
set(gca,set_one_bif{:})
legend('stable sym','','unstable sym','stable non-sym','sym-breaking','FontSize',ft,'FontWeight','normal')%,'FontName','Times New Roman')%,FontWeight='bold')
%\mu_{\rotatebox[origin=c]{180}{T}}
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='bold')
ylabel('$\mu_{\bot}(t)$', ...
       'Interpreter','latex', ...
       'FontSize',17, ...
       'FontWeight','bold');
xlim([0,rp2_x(1)])
title('A','FontSize',ft,'FontWeight','normal')
yticks([-0.02,0,0.02])
xlim([1,20])

% %%%%%%%%%%%%%%%%%%%%%

nexttile
hold on
plot(p_intg.mesh*p_intg.period,p_intg.profile(1:2,:),'LineWidth',3)
plot(0.6,0.9,'.','Color','red','MarkerSize',40)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
title('B','FontSize',ft,'FontWeight','normal')
ylim([0,1])
xlim([0,p_intg.period])

xticks([0.33,0.66])
yticks([0,0.5,1])
%%%%%%%%%%%%%%%%%%
nexttile
hold on 
plot(p_seg.mesh*p_seg.period,p_seg.profile(1:2,:),'LineWidth',3)
plot(0.2,0.9,'.','Color',mbrown,'MarkerSize',40)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','interpreter','latex',FontSize=ft,FontWeight='bold')
title('C','FontSize',ft,'FontWeight','normal')
ylim([0,1])
xlim([0,p_seg.period])
xticks([0.12,0.24])
yticks([0,0.5,1])
nexttile
hold on 
plot(p_bis.mesh*p_bis.period,p_bis.profile(1:2,:),'LineWidth',3)
plot(0.4,0.9,'x','Color','red','MarkerSize',15,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca,  set_one_bif{:})
xlabel('period (s)','FontSize',ft,'FontWeight','normal')
title('D','FontSize',ft,'FontWeight','normal')

ylim([0,1])
xlim([0,0.5])

xticks([0.25,0.50])
yticks([0,0.5,1])
nexttile
hold on 
plot(p_bis2.mesh*p_bis2.period,p_bis2.profile(1:2,:),'LineWidth',3)
plot(0.4,0.9,'x','Color','red','MarkerSize',15,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca,  set_one_bif{:})
xlabel('period (s)','FontSize',ft,'FontWeight','normal')
title('E','FontSize',ft,'FontWeight','normal')
ylim([0,1])
xlim([0,0.5])

xticks([0.25,0.50])
yticks([0,0.5,1])
%% Ploting the maximum of the second peak of u_A over the interval I_{B}
c_A2=zeros(1,12);
c_A2(1)=1;
ua_eval_asymbk=@(p)max(c_A2*dde_coll_eva(p.profile,p.mesh,linspace(0.5,1,1e4),p.degree));
second_peak=zeros(1,length(Symbk_org_br_with_stab.point));
for i=1:length(Symbk_org_br_with_stab.point)
    pi=Symbk_org_br_with_stab.point(i);
    second_peak(i)=ua_eval_asymbk(pi);
end
%%
clrs=colormap(jet);
tmp=@(x)round(x*255+1);
mscal=tmp((second_peak-min(second_peak))/(max(second_peak)-min(second_peak)));
set_one_bif={'LineWidth',1.5,'Box','on','FontSize',8};
ft=11;
figure(320)
clf; 
tiledlayout(2,1)
nexttile
hold on
for i=1:length(second_peak)
    hold on
plot(rp_symbk1(i),df_symbk1(i),'.','Color',clrs(mscal(i),:),'MarkerSize',15)
end
colormap(jet);
clim([min(second_peak),max(second_peak)])
cb=colorbar;
cb.Limits=[min(second_peak),max(second_peak)];
cb.Label.FontSize = 12;
set(gca,  set_one_bif{:})
xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=13,FontWeight='normal')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=13,FontWeight='normal')
po=Symbk_org_br_with_stab.point(550);
plot(rp_symbk1(550),df_symbk1(550),'k.','MarkerSize',35)
ylim([0,1])
yticks([0,0.5,1])
title('A','FontSize',ft,'FontWeight','normal')
c_A=[1,0,0,0,0,0]; % c for u_A: to only extract the values of u_A
sympo_uA_extrema=arrayfun(@(p)dde_coll_roots(p,c_A2,'diff',1)',Symbk_org_br_with_stab.point,'uniformoutput',false);
ua_eval=@(p,t)c_A2*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(Symbk_org_br_with_stab.point),sympo_uA_extrema);
[~,theta_cross]=min(abs(sympomax_ua-(0.5+1e-4)));
%grid on
% figure(310)
% clf
nexttile
hold on
plot(po.mesh*po.period,po.profile(1:2,:),'LineWidth',2)
tv=sympo_uA_extrema{550};
plot(tv(3)*po.period,second_peak(550),'k.','MarkerSize',35)
yline(0.50,'k--','LineWidth',2,'FontSize',14)
yline(second_peak(550),'--','Color',0.5*[1 1 1],'LineWidth',1.5)
yticks([0,0.45,0.50,1])
xticks([0,0.14,0.28])
xlim([0,0.28])
%ylim([0,1])
set(gca,  set_one_bif{:})
xlabel('time (s)','FontSize',ft,'FontWeight','normal')
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','','interpreter','latex',FontSize=ft,FontWeight='normal')
title('B','FontSize',ft,'FontWeight','normal')
%%

load('audi_simulations.mat')


%%
set_one_bif={'LineWidth',2,'Box','on','FontSize',11,'FontWeight','bold'};
fiss=readtable('fission.csv'); fiss=fiss{:,:}; [~,idx]=sort(fiss(:,2)); fiss=fiss(idx,:);
coher=readtable('coherence.csv'); coher=coher{:,:}; [~,idx]=sort(coher(:,2)); coher=coher(idx,:);
%%
figure(9)
clf
hold on
set(gca, set_one_bif{:})
pltcoh=plot3(1000./coher(:,2),2.^(coher(:,1)/12)-1,6*ones(length(coher)),'x','Color','k','LineWidth',1,'MarkerSize',7);
pltfiss=plot3(1000./fiss(:,2),2.^(fiss(:,1)/12)-1,6*ones(length(fiss)),'.','Color','k','LineWidth',1,'MarkerSize',10);
x_text_fiss =7; 
y_text_fiss =0.31; 
%text(x_text_fiss, y_text_fiss, 'fission', 'FontSize', 14, 'Color', 'black','Rotation', 0)
x_text_coh = 15; 
y_text_coh= 0.327;
%text(x_text_coh,y_text_coh,'coherence', 'FontSize', 14, 'Color', 'black','Rotation',0)
text(x_text_fiss, 0.4, '2 (bist)', 'FontSize', 12, 'Color', 'black','Rotation', 0)
text(9, 0.1, '1 (int)', 'FontSize', 12, 'Color', 'black','Rotation', 0)
text(12, 0.6, '3 (seg)', 'FontSize', 12, 'Color', 'black','Rotation', 0)

vtext={'coherence','fission'};
text_vec=[pltcoh(1);pltfiss(1)];
legend(text_vec,vtext,'FontSize',12)
ylim([0,1] )
yticks([0,1])
xlim([min(1000./coher(:,2)),max(1000./coher(:,2))])
set(figure(9), 'Position',[10, 300, 450, 250])
%%
ft=13;
xspace=linspace(-1,1,1e6);
figure(700)
clf 
hold on
plot(xspace,sig(xspace),'k','Linewidth',2)
plot(xspace,SF(xspace),'-','Color',clrs2(1,:),'Linewidth',2)
xline(0.5,'--','Color',0.7*[1 1 1],'LineWidth',2)
xline(0.0,'--','Color',0.7*[1 1 1],'LineWidth',2)
xticks([-1,-0.5,0,0.5,1])
yticks([0,1])
set(gca, set_one_bif{:})
legend('$\theta=0$','$\theta=0.5$','interpreter','latex',FontSize=ft,FontWeight='normal')
%ylabel('$S_{\theta}(x)$','interpreter','latex',FontSize=ft,FontWeight='normal')
%xlabel('$x$','interpreter','latex',FontSize=ft,FontWeight='normal')
set(figure(700), 'Position', [10, 300, 450, 250])
%%
ft=13;
set_one_bif={'LineWidth',2,'Box','on','FontSize',10,'FontWeight','normal'};

figure(1144)
clf;
tiledlayout(3,4,'TileSpacing','compact')
nexttile([3 2])
% %%%%%%%%%%%%%%%%%%%%% Exp range
hold on;
%%%%%%%%%%%%%%%%%% Bifurcation results
plt1=plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b-','color',clrs2(1,:),'LineWidth',4);
%plot(rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plt_thsy=plot(rp_thta(dum==0),df_thta(dum==0),'m-','Linewidth',3);
plt_thsy2=plot(rp_thta(dum>=1),df_thta(dum>=1),'','color',graycolor,'Linewidth',2);
plt_thunsym=plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'--','Color','k','LineWidth',3);
plt_thup=plot(rpp_x,rdff_x,'-.','Color',clrs2(5,:),'LineWidth',2.5);
%%%%%%%% 
df_ax_pr=arrayfun(@(x)x.parameter(in.df),branch_prlarge.point);  %df2_y=arrayfun(@(x)x.parameter(in.df),branch_prlarge2.point);

[~,inds1]=min(abs(rpp_x-38));
p_pr=br_stat_touch.point(inds1);
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',15,'LineWidth',3)
plot(pr_m(t5(1)),df_m(t5(1)),'x','Color',clrs2(3,:),'MarkerSize',15,'LineWidth',3)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'.',...
    'Color',clrs2(7,:),'MarkerSize',35)%,'LineWidth',2)
plot(p_pr.parameter(in.PR),p_pr.parameter(in.df),'x','color','k','MarkerSize',15,'LineWidth',3)
[~,inds2]=min(abs(df_ax_pr-0.05));
p_pr2=branch_prlarge.point(inds2);
plot(p_pr2.parameter(in.PR),p_pr2.parameter(in.df),'.','color',0.45*[1 1 1],'MarkerSize',35)
plt_text={'sym-breaking','touching \theta: stable sym','touching \theta: unstable sym ','touching \theta: stable non-sym','touching \theta from below: stable sym'};
plt_vect=[plt1;plt_thsy;plt_thsy2;plt_thunsym;plt_thup];
legend(plt_vect,plt_text,'FontSize',ft,'FontWeight','normal')
%%%%%%%%%%%%

set(gca, set_one_bif{:})
text(5.5,0.1 ,'1 (int)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'FontWeight','bold')
text(7,0.25,'2 (bist)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
text(15,0.5,'3 (seg)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
text(32,0.05,'4 (sat)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')

xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')
title('A','FontSize',ft,'FontWeight','normal')
ylim([0,1] )
yticks([0,0.5,1])
xlim([0,40])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile%([2 1])
hold on
plot(pt.mesh*pt.period,pt.profile(1:2,:),'LineWidth',3);
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.05,0.9,'x','Color','m','MarkerSize',15,'LineWidth',3)
set(gca, set_one_bif{:})
title('B','FontSize',ft,'FontWeight','normal')
xlabel('time (s)',FontSize=ft,FontWeight='normal')

xlim([0,0.058])
%xticks([0,pt.period/2,pt.period])
xticks([0.029,0.058])
yticks([0,0.5,1])
ylim([0,1])
%grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile%([2 1])
hold on
plot(pt5.mesh*pt5.period,pt5.profile(1:2,:),'LineWidth',3);
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.25,0.9,'x','Color',clrs2(3,:),'MarkerSize',15,'LineWidth',3)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('C','FontSize',ft,'FontWeight','normal')
xlim([0,0.30])
xticks([0.15,0.30])
yticks([0,0.5,1])
ylim([0,1])
%grid on
nexttile%([2 1])
hold on
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',3)
plot(0.275,0.9,'.','Color',clrs2(7,:),'MarkerSize',35)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('D','FontSize',ft,'FontWeight','normal')
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',ft)
xlim([0,0.34])
xticks([0.17,0.34])
ylim([0,1])
yticks([0,0.5,1])%grid on
nexttile
hold on
plot(p_pr.mesh*p_pr.period,p_pr.profile(1:2,:),'LineWidth',3)
plot(0.04,0.9,'x','Color','k','MarkerSize',15,'LineWidth',3)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('E','FontSize',ft,'FontWeight','normal')
xlim([0,0.053])
xticks([0.026,0.053])
ylim([0,1])
yticks([0,0.5,1])
%%%
df_ax_pr=arrayfun(@(x)x.parameter(in.df),branch_prlarge.point);  %df2_y=arrayfun(@(x)x.parameter(in.df),branch_prlarge2.point);

[~,inds2]=min(abs(df_ax_pr-0.05));
p_pr2=branch_prlarge.point(inds2);
nexttile
hold on
plot(p_pr2.mesh*p_pr2.period,p_pr2.profile(1:2,:),'LineWidth',3)
plot(0.04,0.7,'.','Color',0.45*[1 1 1],'MarkerSize',35)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('F','FontSize',ft,'FontWeight','normal')
xlim([0,0.053])
xticks([0.026,0.053])
ylim([0,1])
yticks([0,0.5,1])
%%
%mbranch_df_wbifs

% % tiledlayout(1,2)
% % nexttile
% % plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--','LineWidth',1);
% figure(555)
% clf;
% for i=1:length(mbranch_df_wbifs.point) 
%     subplot(1,2,1)
%     hold on
%     plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--','LineWidth',1);
%     plot(pr_m(i),df_m(i),'k.','MarkerSize',15)
%     hold on
%     subplot(1,2,2)
%     po=mbranch_df_wbifs.point(i);
%     plot(po.mesh,po.profile(1:2,:),'LineWidth',2)
%     yline(0.5,'k--','LineWidth',2)
%     drawnow
% end
%%
inds=find(diff(sign(pr_m-7)));
figure(999)
clf;
subplot(1,3,1)
    hold on
    xline(7,'--','Color',0.5*[1 1 1],'LineWidth',2)

    plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k-','LineWidth',1.5);
    plot(pr_m(inds(1)),df_m(inds(1)),'bx','MarkerSize',12,'LineWidth',1.5)
        plot(pr_m(inds(2)),df_m(inds(2)),'rx','MarkerSize',12,'LineWidth',1.5)
subplot(1,3,2)
    hold on
    po=mbranch_df_wbifs.point(inds(1));
    plot(po.mesh,po.profile(1:2,:),'LineWidth',2)
    plot(0.75,0.6,'bx','MarkerSize',12,'LineWidth',1.5)
    yline(0.5,'k--','LineWidth',2)
    subplot(1,3,3)
    hold on
    po2=mbranch_df_wbifs.point(inds(2));
    plot(po2.mesh,po2.profile(1:2,:),'LineWidth',2)
     plot(0.75,0.6,'rx','MarkerSize',12,'LineWidth',1.5)

    yline(0.5,'k--','LineWidth',2)
 %%
 ft=13;
 ftm=12;
set_one_bif={'LineWidth',2,'Box','on','FontSize',10,'FontWeight','normal'};
 figure(11)
clf;
tiledlayout(3,4,'TileSpacing','compact')
nexttile([3 2])
% %%%%%%%%%%%%%%%%%%%%% Exp range
hold on;
%%%%%%%%%%%%%%%%%% Bifurcation results
plt1=plot(rp_symbk1(nunst_bk==0),df_symbk1(nunst_bk==0),'b-','color',clrs2(1,:),'LineWidth',4);
%plot(rp_symbk1(nunst_bk>=1),df_symbk1(nunst_bk>=1),'k.','MarkerSize',8)
plt_thsy=plot(rp_thta(dum==0),df_thta(dum==0),'m-','Linewidth',3);
plt_thsy2=plot(rp_thta(dum>=1),df_thta(dum>=1),'','color',graycolor,'Linewidth',2);
plt_thunsym=plot(pr_m(nunst_m==0),df_m(nunst_m==0),'k--',...
pr_m(nunst_m>=1),df_m(nunst_m>=1),'--','Color','k','LineWidth',3);
plt_thup=plot(rpp_x,rdff_x,'-.','Color',clrs2(5,:),'LineWidth',2.5);
%%%%%%%% 
df_ax_pr=arrayfun(@(x)x.parameter(in.df),branch_prlarge.point);  %df2_y=arrayfun(@(x)x.parameter(in.df),branch_prlarge2.point);

[~,inds1]=min(abs(rpp_x-38));
p_pr=br_stat_touch.point(inds1);
plot(rp_thta(it),df_thta(it),'mx','MarkerSize',ftm,'LineWidth',2.5)
plot(pr_m(inds(1)),df_m(inds(1)),'x','Color',clrs2(3,:),'MarkerSize',ftm,'LineWidth',2.5)
plot(pr_m(inds(2)),df_m(inds(2)),'x','Color',clrs2(5,:),'MarkerSize',ftm,'LineWidth',2.5)
plot(non_symper.point(it3).parameter(in.PR),non_symper.point(it3).parameter(in.df),'.',...
    'Color',clrs2(7,:),'MarkerSize',25)%,'LineWidth',2)
plot(p_pr.parameter(in.PR),p_pr.parameter(in.df),'x','color','k','MarkerSize',ftm,'LineWidth',2.5)
[~,inds2]=min(abs(df_ax_pr-0.05));
p_pr2=branch_prlarge.point(inds2);
plot(p_pr2.parameter(in.PR),p_pr2.parameter(in.df),'.','color',0.45*[1 1 1],'MarkerSize',25)
plt_text={'sym-breaking','touching \theta: stable sym','touching \theta: unstable sym ','touching \theta: stable non-sym','touching \theta from below: stable sym'};
plt_vect=[plt1;plt_thsy;plt_thsy2;plt_thunsym;plt_thup];
legend(plt_vect,plt_text,'FontSize',ft,'FontWeight','normal')
%%%%%%%%%%%%

set(gca, set_one_bif{:})
text(5.5,0.1 ,'1 (int)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'FontWeight','bold')
text(7,0.25,'2 (bist)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
text(15,0.5,'3 (seg)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')
text(32,0.05,'4 (sat)', 'FontSize', 13, 'Color', 0.4*[1 1 1],'Rotation',0,'FontWeight','bold')

xlabel('$r_\mathrm{p}$','interpreter','latex',FontSize=17,FontWeight='bold')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=17,FontWeight='bold')
title('A','FontSize',ft,'FontWeight','normal')
ylim([0,1] )
yticks([0,0.5,1])
xlim([0,40])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile%([2 1])
hold on
plot(pt.mesh*pt.period,pt.profile(1:2,:),'LineWidth',3);
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.05,0.9,'x','Color','m','MarkerSize',ftm,'LineWidth',2)
set(gca, set_one_bif{:})
title('B','FontSize',ft,'FontWeight','normal')
%xlabel('time (s)',FontSize=ft,FontWeight='normal')

xticks([0,pt.period/2,pt.period])
xlim([0,pt.period])
xtickformat('%3.2g')
yticks([0,0.5,1])
ylim([0,1])
%grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile%([2 1])
hold on
plot(po.mesh*po.period,po.profile(1:2,:),'LineWidth',3);
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.25,0.9,'x','Color',clrs2(3,:),'MarkerSize',ftm,'LineWidth',2.5)
set(gca, set_one_bif{:})
%xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('C','FontSize',ft,'FontWeight','normal')
%xlim([0,0.30])
xlim([0,po.period])
xticks([0,po.period/2,po.period])
xtickformat('%3.2g')
yticks([])
ylim([0,1])
nexttile%([2 1])
hold on
plot(po2.mesh*po2.period,po2.profile(1:2,:),'LineWidth',3);
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.25,0.9,'x','Color',clrs2(5,:),'MarkerSize',ftm,'LineWidth',2.5)
set(gca, set_one_bif{:})
%xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('D','FontSize',ft,'FontWeight','normal')
%xlim([0,0.30])
xlim([0,po2.period])
xticks([0,po2.period/2,po2.period])
xtickformat('%3.2g')
%xticks([0.15,0.30])
yticks([0,0.5,1])
ylim([0,1])
nexttile
hold on
plot(p_pr.mesh*p_pr.period,p_pr.profile(1:2,:),'LineWidth',3)
plot(0.04,0.9,'x','Color','k','MarkerSize',ftm,'LineWidth',2.5)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
%xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('E','FontSize',ft,'FontWeight','normal')
xlim([0,p_pr.period])
xticks([0,p_pr.period/2,p_pr.period])
xtickformat('%3.2g')
ylim([0,1])
yticks([])
nexttile%([2 1])
hold on
plot(pt3.mesh*pt3.period,pt3.profile(1:2,:),'LineWidth',3)
plot(0.275,0.9,'.','Color',clrs2(7,:),'MarkerSize',25)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('F','FontSize',ft,'FontWeight','normal')
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','Interpreter','latex','FontSize',ft)
xlim([0,pt3.period])
xticks([0,pt3.period/2,pt3.period])
xtickformat('%3.2g')
ylim([0,1])
yticks([0,0.5,1])%grid on

%%%
df_ax_pr=arrayfun(@(x)x.parameter(in.df),branch_prlarge.point);  %df2_y=arrayfun(@(x)x.parameter(in.df),branch_prlarge2.point);

[~,inds2]=min(abs(df_ax_pr-0.05));
p_pr2=branch_prlarge.point(inds2);
nexttile
hold on
plot(p_pr2.mesh*p_pr2.period,p_pr2.profile(1:2,:),'LineWidth',3)
plot(0.04,0.7,'.','Color',0.45*[1 1 1],'MarkerSize',25)%,'LineWidth',3)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
set(gca, set_one_bif{:})
xlabel('time (s)',FontSize=ft,FontWeight='normal')
title('G','FontSize',ft,'FontWeight','normal')
xlim([0,p_pr2.period])
xticks([0,p_pr2.period/2,p_pr2.period])
xtickformat('%3.2g')
ylim([0,1])
yticks([])