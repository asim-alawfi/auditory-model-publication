clear;
    base=[pwd(),'\..\..\ddebiftool_snapshot_23October2022\'];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco']);
%%
ntau=2;
parnames={'a','b','tau','tau_i','D','PR','df','TD','lambda','theta','c', 'm','ph1','ph2'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
%% Define system using symbolic algebra
% define arbitrary variable names
u_A=sym('u_A',[1,ntau+1]); % u_A
u_B=sym('u_B',[1,ntau+1]); % u_B
s_a=sym('s_a',[1,ntau+1]); % s_A
s_b=sym('s_b',[1,ntau+1]); % s_B
% Hopf variables
x=sym('x',[1,ntau+1]); % x
y=sym('y',[1,ntau+1]); % x
syms(parnames{:});
par=sym(parnames);
%% Sigmoid function
sig=@(x)sigmfunc(x,lambda);
S=@(x)sig(x-theta); %sigmoid function with threshold activity \theta
%%
d=c*(1-nthroot(df,m)); % rescale 'd' according to the formal in the paper
iA=c*sig(x(1))*sig(-x(3))+d*sig(-x(1))*sig(x(3));
iB=d*sig(x(1))*sig(-x(3))+c*sig(-x(1))*sig(x(3));
%%
rhs=[(-u_A(1)+S(a*u_B(1)-b*s_b(2)+iA))/tau;... 
     (-u_B(1)+S(a*u_A(1)-b*s_a(2)+iB))/tau;...
    S(u_A(1))*(1-s_a(1))/tau - s_a(1)/tau_i;...
    S(u_B(1))*(1-s_b(1))/tau - s_b(1)/tau_i;...
    ph1*x(1)-pi*PR*y(1)- x(1)*(x(1)^2 + y(1)^2);...
    pi*PR*x(1)+ph1*y(1)- y(1)*(x(1)^2 + y(1)^2)];
%% Defferentiate and generate code, exporting it to sym_auditory ( creat file with the same name, and continue)
[fstr,erives]=dde_sym2funcs(rhs,[u_A;u_B;s_a;s_b;x;y],par,'filename','symbolic_auditory_with_symmetry_version',...
    'maxorder',2,'directional_derivative',true);

