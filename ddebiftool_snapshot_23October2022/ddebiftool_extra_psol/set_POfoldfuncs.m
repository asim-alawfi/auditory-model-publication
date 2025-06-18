function [pfuncs,extra_freepar,initfuncs]=set_POfoldfuncs(funcs,point,method,varargin)
%% set up extended systems, numbering and values of additional parameters and artificial delays
%% default extra conditions
standardcond={'POEV1_norm',true,'POEV1_phase_condition',true};
%% process options
default={'sys_deri',1e-4,'sys_dtau',1e-4,'hjac',1e-4,'df_deriv',false,...
    'nullparind',zeros(0,1),'usercond',cell(0,1),'initcond',cell(0,1),...
    standardcond{:}}; %#ok<CCAT>
options=dde_set_options(default,varargin,'pass_on');
if isempty(options.initcond)
    options.initcond=options.usercond;
end
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
if isfield(point,'profile')
    ip.dim=size(point.profile,1);    
elseif isfield(point,'x')
    ip.dim=length(point.x);
end
ip.extdim=ip.dim;
ip.nuserpar=length(point.parameter); % number of original system parameters
%% extend problem
ip.beta=ip.nuserpar+1;       % location of add. parameter beta
ip.period=ip.nuserpar+2;     % location of copy of period
nnull=length(options.nullparind);
ip.nullparind=NaN(nnull,2); % location of parameters for nullspace vector
ip.nullparind(:,1)=options.nullparind;
ip.nullparind(:,2)=ip.period+(1:nnull)';
get_comp=@(p,component)extract_from_POEV1(p,component,ip);
user_lhs_num=funcs.lhs_matrix(ip.dim);
lhs_num=kron(eye(2),user_lhs_num);
fun_args={'x_vectorized',funcs.x_vectorized,...
    'lhs_matrix',lhs_num,'hjac',options.hjac};
% indices of additional delays and relations needed for extended system:
delay_cond={};
ip.ext_tau=[];
if ~funcs.tp_del % constant delay
    ip.orig_tau=funcs.sys_tau(); % system delays
    % additional delays needed for extended system:
    ip.orig_ntau=length(ip.orig_tau);
    ip.ext_tau=ip.period+nnull+(1:length(ip.orig_tau));
    delay_cond={@(p,pref)sys_cond_POEV1_extradelay(p,ip.orig_tau,ip.orig_tau,ip.ext_tau)};
    %% set up functions of extended system
    sys_rhs=@(x,p)sys_rhs_POEV1(x,p,ip,funcs,options.sys_deri);
    sys_tau=@()[ip.orig_tau,ip.ext_tau];
    %sys_deri=@(x,p,nx,np,v)...
    %    sys_deri_POEV1(x,p,nx,np,v,options.hjac,...
    %    beta_ind,period_ind,dim,xtau_ind,funcs,sys_rhs,options.df_deriv);
    pfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,...%disabled forn now: 'sys_deri',sys_deri,...
        fun_args{:});
else % state-dependent delay
    % indices of additional delays needed for extended system:
    ip.orig_ntau=funcs.sys_ntau();  % number of state-dependent delays
    %% set up functions of extended system
    sys_rhs=@(x,p)sys_rhs_SD_POEV1(x,p,ip,funcs,options.sys_deri,options.sys_dtau);
    sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p,ip,funcs);
    sys_ntau=@()ip.orig_ntau*2;
    pfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,...
        'sys_ntau',sys_ntau,...
        fun_args{:},'hjac',options.hjac);
end
%% required amendments of structures for extended system
pfuncs.delayed_derivs=[zeros(1,ip.orig_ntau),ones(1,ip.orig_ntau)];
nullfreepar=[ip.beta,ip.nullparind(:,2)'];
extra_freepar=[nullfreepar,ip.period,ip.ext_tau];
%% extended delays are for derivatives
pfuncs.sys_cond_reference=true;
pfuncs.get_comp=get_comp;
pfuncs.kind='POfold';
pfuncs.userfuncs=funcs;
pfuncs.usermethod=method;
pfuncs.ip=ip;
%% add requested and standard extra conditions
sys_cond_extra={@(p,pref)sys_cond_POEV1_phase_condition(p,ip.dim),...
    @(p,pref)sys_cond_POEV1_norm(p,ip.dim,nullfreepar,'res',1,'period',false)};
psys_cond=@(p,pref)dde_sys_cond_collect(pfuncs,p,pref,...
    {@(p,pref)sys_cond_coll_fixperiod(p,ip.period),...
    delay_cond{:},...
    sys_cond_extra{[options.POEV1_phase_condition,options.POEV1_norm]},...
    options.usercond{:}}); %#ok<CCAT>
initfuncs=pfuncs;
isys_cond=@(p,pref)dde_sys_cond_collect(initfuncs,p,pref,...
    {sys_cond_extra{[options.POEV1_phase_condition,false]},...
    options.initcond{:}}); %#ok<CCAT>
pfuncs.sys_cond=psys_cond;
initfuncs.sys_cond=isys_cond;
end
