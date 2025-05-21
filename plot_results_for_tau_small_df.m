clear
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
  %%
  load('sym_breaking_small_tau_coh_and_fiss.mat')
  %load('sym_breaking_small_tau_part3.mat')
 % load('one_parameter_bif_small_tau_df.mat')
  % load('one_br_bif_in_df_small_tau.mat')
  % load('one_parbranch_tau_eqls02_includes_initialcomputation3.mat')
  % load('branch_touching_theta_s2.mat')
  %% We plot int_{0}^{1/2} u_A(t)-u_B(t+1/2) dt by setting:
  % Syms=dde_lincond_struct(size(asym_brtinyTau2.point(1).profile,1),'profile',...
  %   'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'condprojint',[0,0.5]);
%%
asdf2_tinytau2=arrayfun(@(x)x.parameter(in.df),asym_brtinyTau2_df.point);
Sint_A=dde_lincond_struct(size(asym_brtinyTau2_df.point(1).profile,1),'profile','trafo',Rsym,...
    'shift',[1,2],'condprojmat',[1,1,0,0,0,0],'stateproj',sproj,'condprojint',[0,0.5]);
% Sint_B=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
%     'shift',[1,4],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.25,0.75]);
yax_sym_df=arrayfun(@(x)dde_psol_lincond(x,Sint_A),po_symmetry_tinyTau2_df.point);
yax_nonsym_df=arrayfun(@(x)dde_psol_lincond(x,Sint_A),asym_brtinyTau2_df.point);
%% plotting results for one-parameter bifurcation in $d_f$
% figure(50)
% clf; hold on 
% grid on
% plot(df2_tinytau2(nunst_sym_tinytau2_df==0),yax_sym_df(nunst_sym_tinytau2_df==0),'b.',...
%     df2_tinytau2(nunst_sym_tinytau2_df>=1),yax_sym_df(nunst_sym_tinytau2_df>=1),'k.','LineWidth',1)
% plot(asdf2_tinytau2(nunst_nonsym_tinytau2_df==0),yax_nonsym_df(nunst_nonsym_tinytau2_df==0),'g.','MarkerSize',10)
% plot(asdf2_tinytau2(nunst_nonsym_tinytau2_df>=1),yax_nonsym_df(nunst_nonsym_tinytau2_df>=1),'k.','MarkerSize',30)
%% plotting results for two-parameter bifurcation in $(r_p,d_f)$

rp_s_tau=arrayfun(@(x)x.parameter(in.PR),symbk_brtu.point);
df_s_tau=arrayfun(@(x)x.parameter(in.df),symbk_brtu.point);
rp_ss_tau=arrayfun(@(x)x.parameter(in.PR),symbk_brtu2.point);
df_ss_tau=arrayfun(@(x)x.parameter(in.df),symbk_brtu2.point);
%
rptt_tau=arrayfun(@(x)x.parameter(in.PR),mbranch.point);
dftt_tau=arrayfun(@(x)x.parameter(in.df),mbranch.point);
%
load('tracking_threshold_crossing_asymmetric_try2.mat')
% Pick points for time profile
clrs2=lines();
[~,i1_intg]=min(abs(df2_tinytau2-0.2));
[~,i2_seg]=min(abs(df2_tinytau2-0.6));
  %i1_intg=pick_po(1);
pick_bispo=find(diff(sign(asdf2_tinytau2-0.37)));
  i1_bis=pick_bispo(1);
  i2_bis=pick_bispo(2);
ft=13;
mft=35;
set_bif={'LineWidth',2,'Box','on','FontSize',12};
p_i1_intg=po_symmetry_tinyTau2_df.point(i1_intg);
p_i1_bis=asym_brtinyTau2_df.point(i1_bis);
p_i2_seg=po_symmetry_tinyTau2_df.point(i2_seg);
figure(1001)
clf
tiledlayout(2,4)%,'TileSpacing','loose')
nexttile([2,2])
hold on
plot(rp_ss_tau,df_ss_tau,'-',rp_s_tau,df_s_tau,'-','LineWidth',3,'Color',clrs2(1,:))%,'MarkerSize',5)
plot(rptt_tau,dftt_tau,'k--','LineWidth',3)%,'MarkerSize',8)
%
plot(pr_m,df_m,'Color',0.8*[1 1 1],'LineWidth',3)
plot(p_i1_intg.parameter(in.PR),p_i1_intg.parameter(in.df),'k.','MarkerSize',mft)
plot(p_i2_seg.parameter(in.PR),p_i2_seg.parameter(in.df),'.','Color',clrs2(2,:),'MarkerSize',mft)
plot(p_i1_bis.parameter(in.PR),p_i1_bis.parameter(in.df),'.','Color',clrs2(5,:),'MarkerSize',mft)
set(gca, set_bif{:})
legend('sym-breaking','touching \theta: stable non-sym','base case','','FontSize',ft,'FontWeight','normal')
xlabel('$r_\mathrm{p}$','interpreter','latex','FontSize',15,FontWeight='normal')
ylabel('$d_\mathrm{f}$','interpreter','latex',FontSize=15,FontWeight='normal')
x_text_intg = 17.5%p_i1_intg.parameter(in.PR)+1; 
y_text_intg = p_i1_intg.parameter(in.df)+0.03; 
text(x_text_intg, y_text_intg, '1 (int)', 'FontSize', ft, 'Color', 0.4*[1 1 1],'FontWeight','bold')
x_text_seg =23% p_i2_seg.parameter(in.PR)+1; 
y_text_seg = 0.7%p_i2_seg.parameter(in.df); 
text(x_text_seg, y_text_seg,'3 (seg)', 'FontSize', ft, 'Color', 0.4*[1 1 1],'FontWeight','bold')
x_text_bis= 17;%p_i1_bis.parameter(in.PR)-4.3; 
y_text_bis = 0.4%p_i1_bis.parameter(in.df)+0.015; 
text(x_text_bis, y_text_bis,'2 (bist)', 'FontSize', ft, 'Color', 0.4*[1 1 1],'FontWeight','bold')

ylim([0,1])
yticks([0,0.5,1])
%grid on
%title('b','FontSize',18,'FontWeight','bold')
xlim([0,40])
title('A','FontSize',ft,'FontWeight','normal')
%
%%%%%%%%%%%%%
nexttile
hold on
plot(p_i1_intg.mesh*p_i1_intg.period,p_i1_intg.profile(1:2,:),'LineWidth',3)
plot(p_i1_intg.mesh*p_i1_intg.period,p_i1_intg.profile(3:4,:),'LineWidth',2)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.07,0.85,'k.','MarkerSize',mft)
set(gca, set_bif{:})
xlim([0,p_i1_intg.period])
ylim([0,1])
xticks([p_i1_intg.period/2,p_i1_intg.period])
yticks([0,0.5,1])
title('B','FontSize',ft,'FontWeight','normal')
xlabel('time (s)','FontSize',ft,'FontWeight','normal')
%%%%%%%%%%%%%
nexttile
hold on
plot(p_i1_bis.mesh*p_i1_bis.period,p_i1_bis.profile(1:2,:),'LineWidth',3)
plot(p_i1_bis.mesh*p_i1_bis.period,p_i1_bis.profile(3:4,:),'LineWidth',2)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.07,0.85,'.','Color',clrs2(5,:),'MarkerSize',mft)
set(gca, set_bif{:})
ylim([0,1])
xlim([0,p_i1_bis.period])
xticks([p_i1_bis.period/2,p_i1_bis.period])
yticks([0,0.5,1])
%xlabel('period','FontSize',20,'FontWeight','bold')
title('C','FontSize',ft,'FontWeight','normal')
xlabel('time (s)','FontSize',ft,'FontWeight','normal')

%%%%%%%%%%
nexttile
hold on
plot(p_i2_seg.mesh*p_i2_seg.period,p_i2_seg.profile(1:2,:),'LineWidth',3)
plot(p_i2_seg.mesh*p_i2_seg.period,p_i2_seg.profile(3:4,:),'LineWidth',2)
yline(0.5,'k--','LineWidth',2,'FontSize',14)
plot(0.07,0.73,'.','Color',clrs2(2,:),'MarkerSize',mft)
set(gca, set_bif{:})
xlim([0,p_i2_seg.period])
ylim([0,1])
xticks([p_i2_seg.period/2,p_i2_seg.period])
yticks([0,0.5,1])
xlabel('time (s)','FontSize',ft,'FontWeight','normal')
title('D','FontSize',ft,'FontWeight','normal')
legend('$u_\mathrm{A}$','$u_\mathrm{B}$','$s_\mathrm{A}$','$s_\mathrm{B}$','','Interpreter','latex','FontSize',ft,'FontWeight','bold')
