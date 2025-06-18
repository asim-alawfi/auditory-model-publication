clear
  base=[pwd(),'\..\..\ddebiftool_snapshot_23October2022\'];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco']);
format compact
%%                                            Parameters
parnames={'a','b','tau','tau_i','D','PR','df','TD','lambda','theta','c', 'm', 'ph1', 'ph2'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
%%      vecotr parameters  is required when fixed df=0.1
par([in.a, in.b, in.tau, in.tau_i, in.D, in.PR, in.df, in.TD,in.lambda,in.theta, in.c, in.m, in.ph1, in.ph2])=...
    [2,    2.8,   0.025,   0.25,   0.015, 10,     0.1,  0.022,   30,       0.5,    5.5, 6, 1, pi];
%%      vecotr parameters  is required when fixed df=0.73
par_tinytau=par;
par_tinytau(in.df)=0.73;
par_tinytau(in.tau)=0.0025;
par_tinytau(in.tau_i)=0.25;
par_tinytau(in.PR)=25;
%%                                             initial condition(s) for dde23
x0=[1 0 1 0 0 1];
tspan=linspace(0,3,1e7);
%%                                                 Sigmoid functions
sig=@(x)sigmfunc(x,par(in.lambda));
SF=@(x)sig(x-par(in.theta)); %
%% dde23 solution with fixed df=0.1, where Hopf normal form is included %%
lag=[par(in.D), par(in.TD)];
grid on
%%
sol23_df_tinytau2=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,par_tinytau,in),lag,x0,tspan, ddeset('RelTol',1e-7,'Events',@(t,y,z)event(y)));
%% the computation when using force function
figure(3032)
clf; hold on
plot(sol23_df_tinytau2.x,sol23_df_tinytau2.y(1,:),sol23_df_tinytau2.x,sol23_df_tinytau2.y(2,:),...
    sol23_df_tinytau2.x,sol23_df_tinytau2.y(3,:),sol23_df_tinytau2.x,sol23_df_tinytau2.y(4,:),'LineWidth',1);
grid on
ylim([0,1.5])
%% funcs_symb is required when only using ''sys_cond''
% funcs_symb=set_symfuncs(@sym_auditory_with_Hopf_normal_form,'sys_tau',@()[in.D, in.TD],...
%     'sys_cond',@sys_cond);
% %% funcs_symb is required when using ''sys_cond_symmetry''
% funcs_symb_sm=set_symfuncs(@sym_auditory_sm,'sys_tau',@()[in.D, in.TD],...
%     'sys_cond',@sys_cond_symmetry);
%% save the file 
save('dde23_small_tau.mat')




















