%% Computation of non-symmetric solutions (original case)
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
 load('branch_of_sympos_original_case_and_threshold_crossing_coj.mat')
 %%
figure(2);clf;hold on;grid on
plot(rp2_x(nunst2_sym==0),max2_y(nunst2_sym==0),'bo',rp2_x(nunst2_sym>=1),max2_y(nunst2_sym>=1),'kx')%,'LineWidth',1)
plot(rp2_x(nunst2_sym==0),min2_y2(nunst2_sym==0),'bo',rp2_x(nunst2_sym>=1),min2_y2(nunst2_sym>=1),'kx')%,'LineWidth',1)
%
%%
chang2_stb_sym=find(diff(nunst2_sym));
% branching off to the asymmetric solutions continuation in PR 
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,3)','condprojmat',[1,0,0,0,0,0]);
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[pffuncs,nonsymper,suc_v]=SetupPsol(funcs_audi,po2_symmetry_wbifs,chang2_stb_sym(2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.PR,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.01; 0,0.05]);
%%
nonsymper.parameter.max_step(4)=0.05;
nonsymper.parameter.max_step(end)=0.05;
figure(98)
clf;
hold on
nonsymper=br_contn(pffuncs,nonsymper,180);
%%
[nonsymper_wbifs,uns_um,uns_bifs,uns_bifind]=MonitorChange(pffuncs,nonsymper,'print_residual_info',0);
%[uns_um,uns_dom,uns_triv]=GetStability(nonsymper_wbifs,'exclude_trivial',true);
rp_unsym=arrayfun(@(x)x.parameter(in.PR),nonsymper_wbifs.point);
y_unsym=arrayfun(@(x)max(x.profile(1,:)),nonsymper_wbifs.point);
y2_unsym=arrayfun(@(x)min(x.profile(1,:)),nonsymper_wbifs.point);
%%
[~,it_intg]=min(abs(rp2_x-3));
[~,it_seg]=min(abs(rp2_x-8.15));
it_bis=find(diff(sign(rp_unsym-4.0)));
figure(212);clf;hold on;%grid on
v1=plot(rp2_x(nunst2_sym==0),max2_y(nunst2_sym==0),'b.',rp2_x(nunst2_sym>=1),max2_y(nunst2_sym>=1),'kx','MarkerS',10)%,'LineWidth',1)
%v2=plot(rp2_x(it_intg),max2_y(it_intg),'r.','MarkerSize',15)%plot(rp2_x(nunst2_sym==0),min2_y2(nunst2_sym==0),'bo',rp2_x(nunst2_sym>=1),min2_y2(nunst2_sym>=1),'kx')%,'LineWidth',1)
%
figure(212)
hold on
v3=plot(rp_unsym(uns_um==0),y_unsym(uns_um==0),'r.','MarkerS',10);
   % rp_unsym(uns_um>=1),y_unsym(uns_um>=1),'ko','LineWidth',1)
   v4=plot(rp2_x(it_intg),max2_y(it_intg),'m.',rp2_x(it_seg),max2_y(it_seg),...
       rp_unsym(it_bis),y_unsym(it_bis),'g.','MarkerSize',35);
   plot(rp2_x(it_seg),max2_y(it_seg),'m.','MarkerSize',35)
      %plot(rp_unsym(it_bis),y_unsym(it_bis),'m.','MarkerSize',25)
     
legend([v1(1),v1(2),v3(1)],...
 {'stable symmetric POs','unstable symmetric POs', 'stable asymmetric POs'})
xlabel('r_p')
ylabel('max u_A')
yticks([])
set(gca,'FontWeight','bold')
%% 
save('branching_towrds_asymmetric_sols_original_case_coj.mat')
%save('branching_off_to_unsymmetris_sols_original_case.mat')
