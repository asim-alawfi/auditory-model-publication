function [trfuncs,extra_freepar,initfuncs]=set_torusfuncs(funcs,point,method,varargin)
%% set up extended systems, numbering and values of additional parameters and artificial delays
%% default extra conditions
standardcond={'Torus_norm',true,'Torus_phase_condition',true};
%% process options
default={'sys_deri',1e-4,'sys_dtau',1e-4,...
    'hjac',1e-4,'biftype','torus',...
    'usercond',cell(0,1),...
    'initcond',{},...
    standardcond{:},...
    'nullparind',zeros(0,1)}; %#ok<CCAT>
options=dde_set_options(default,varargin,'pass_on');
if isempty(options.initcond)
    options.initcond=options.usercond;
end
%% set up numbering and values of additional parameters
ip.dim=size(point.profile,1);    % dimension of original problem
ip.xrg=1:ip.dim;
ip.re=ip.dim+(1:ip.dim);
ip.im=ip.re(end)+(1:ip.dim);
ip.null=[ip.re,ip.im];
ip.nuserpar=length(point.parameter); % number of original system parameters
ip.nullparind(:,1)=options.nullparind;
ip.omega=ip.nuserpar+1;             % location of add. parameter omega
ip.period=ip.omega+1;            % location of add. parameter (equal to period)
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
%% set up functions of extended system
user_lhs_num=funcs.lhs_matrix(ip.dim);
lhs_num=kron(eye(3),user_lhs_num);
fun_args={'x_vectorized',funcs.x_vectorized,...
    'lhs_matrix',lhs_num};
if ~funcs.tp_del 
    %% constant delays
    ip.orig_tau=funcs.sys_tau();
    ip.orig_ntau=length(ip.orig_tau);
    sys_rhs=@(x,p)sys_rhs_TorusBif(x,p,ip,funcs,options.sys_deri,user_lhs_num);
    sys_deri=@(x,p,nx,np,v)...
        sys_deri_TorusBif(x,p,nx,np,v,options.hjac,ip,funcs,struct('sys_rhs',sys_rhs));
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',@()ip.orig_tau,...
        'sys_deri',sys_deri,fun_args{:});
else
    %% state dependent delay
    ip.orig_ntau=funcs.sys_ntau();  % number of state-dependent delays
    % additional delays needed for extended system:
    sys_rhs=@(x,p)sys_rhs_SD_TorusBif(x,p,ip,funcs,options.sys_deri,options.sys_dtau);
    sys_ntau=@()ip.orig_ntau*2;
    sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p,ip,funcs);
    %sys_dtau=@(itau,x,p,nx,np)sys_dtau_SD_PObif(itau,x,p,nx,np,funcs.sys_dtau,...
    %    dim,xtau_ind);
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,'sys_ntau',sys_ntau,...
        fun_args{:},'hjac',options.hjac);
    trfuncs.delayed_derivs=[zeros(1,ip.orig_ntau),ones(1,ip.orig_ntau)];
end
get_comp=@(p,component)extract_from_tr(p,component,options.biftype,ip);
trfuncs.sys_cond_reference=true;
trfuncs.get_comp=get_comp;
trfuncs.kind=options.biftype;
trfuncs.userfuncs=funcs;
trfuncs.usermethod=method;
trfuncs.ip=ip;
initfuncs=trfuncs;
sys_cond_extra={@(p,pref)sys_cond_Torus_phase_condition(p,pref,ip.dim),...
    @(p,pref)sys_cond_Torus_norm(p,ip.dim)};
trfuncs.sys_cond=@(p,pref)dde_sys_cond_collect(trfuncs,p,pref,...
    {sys_cond_extra{[options.Torus_phase_condition,options.Torus_norm]},...
    @(p,pref)sys_cond_coll_fixperiod(p,ip.period),...
    options.usercond{:}}); %#ok<CCAT>
initfuncs.sys_cond=@(p,pref)dde_sys_cond_collect(initfuncs,p,pref,...
    options.initcond); 
extra_freepar=[ip.omega,ip.period];
end