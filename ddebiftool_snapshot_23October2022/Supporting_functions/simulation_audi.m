clear;
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
par1=par;
par1(in.df)=0.73;
par2=par;
par2(in.df)=0.3;
%%                                             initial condition(s) for dde23
x0=[1 0 1 0 0 1];
%%                                                 Sigmoid functions
sig=@(x)sigmfunc(x,par(in.lambda));
SF=@(x)sig(x-par(in.theta)); %
%% dde23 solution with fixed df=0.1, where Hopf normal form is included %%
lag=[par(in.D), par(in.TD)];
sol23=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,par,in),lag,x0,[0,4], ddeset('RelTol',1e-7,'Events',@(t,y,z)event(y)));
%
%% dde23 solution with fixed df=0.73, where Hopf normal form is included %%
period1=sol23.xe(end)-sol23.xe(end-1);

figure(301)
clf; hold on
plot(sol23.x,sol23.y(1:4,:),'LineWidth',2);
%xlim([4*period1+period1,4*period1+period1])
grid on
%%
sol23_df=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,par1,in),lag,x0,[0,4], ddeset('RelTol',1e-7,'Events',@(t,y,z)event(y)));
%% the computation when using force function
figure(303)
clf; hold on
plot(sol23_df.x,sol23_df.y(1:4,:),'LineWidth',2);
grid on
sol23_df2=dde23(@(t,y,SD)audi_hopf_rhs(y,SD,SF,sig,par2,in),lag,x0,[0,4], ddeset('RelTol',1e-7,'Events',@(t,y,z)event(y)));
%% the computation when using force function
figure(304)
clf; hold on
plot(sol23_df2.x,sol23_df2.y(1:4,:),'LineWidth',2);
grid on
save('audi_simulations.mat')