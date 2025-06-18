clear;
 addpath('\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool\',...
     '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_extra_psol\',...
     '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_utilities\',...
     '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_extra_rotsym\',...
    '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_extra_symbolic\',...
    '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_extra_nmfm\',...
   '\Users\aa1045\mydesktop\dde_biftool_snapshot\ddebiftool_coco\')
format compact
%%                                            Parameters
parnames={'a','b','tau','tau_i','D','PR','df','TD','lambda','theta','c', 'm', 'ph1', 'ph2'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
%%      vecotr parameters  is required when fixed df=0.1
par([in.a, in.b, in.tau, in.tau_i, in.D, in.PR, in.df, in.TD,in.lambda,in.theta, in.c, in.m, in.ph1, in.ph2])=...
    [2,    2.8,   0.0025,   0.25,   0.04, 15,     0.2,  0.022,   30,       0.6,    5.5, 6, 1, pi];
%%      vecotr parameters  is required when fixed df=0.73

%%                                             initial condition(s) for dde23
x0=[1 0 1 0 0 1];
%%                                                 Sigmoid functions
sig=@(x)sigmfunc(x,par(in.lambda));
SF=@(x)sig(x-par(in.theta)); %
%% dde23 solution with fixed df=0.1, where Hopf normal form is included %%
lag=[par(in.D), par(in.TD)];
tspan=linspace(0,5,1e6);
dde23_d_large=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,par,in),lag,x0,tspan, ddeset('RelTol',1e-7,'AbsTol',1e-7, 'Events',@(t,y,z)event(y)));
%%
% %% funcs_symb is required when only using ''sys_cond''
% funcs_symb=set_symfuncs(@sym_auditory_with_Hopf_normal_form,'sys_tau',@()[in.D, in.TD],...
%     'sys_cond',@sys_cond);
% %% funcs_symb is required when using ''sys_cond_symmetry''
% funcs_symb_sm=set_symfuncs(@sym_auditory_sm,'sys_tau',@()[in.D, in.TD],...
%     'sys_cond',@sys_cond_symmetry);
%% save the file 
save('sol23_delaylarge.mat')
%%
figure(8)
clf; hold on
plot(dde23_d_large.x,dde23_d_large.y(1:4,:),'LineWidth',2)
yline(par(in.theta),'--','LineWidth',1)
xlim([1.85,2])


















