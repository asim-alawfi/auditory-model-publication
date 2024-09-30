function [poev1ini,sv,vpoint]=POEV1Init(funcs,point,method,varargin)
%% crude initial guess for fold of periodic orbits
%
%
ip=funcs.ip;
internalpar=ip.nullparind(:,1)';
nint=length(internalpar);
default={'v_scal',@(p,pref)sys_cond_POEV1_norm(p,ip.dim,[ip.beta,ip.nullparind(:,2)'],...
    'period',false,'res',0),'nulldim',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ptini=point;
ptini.profile(ip.dim+(1:ip.dim),:)=0;
if ~funcs.tp_del
    ptini.parameter(ip.ext_tau)=point.parameter(ip.orig_tau);
end
ivarpar=[ip.beta,ip.nullparind(:,2)'];
ptini.parameter(ivarpar)=0;
ptini.parameter(ip.period)=ptini.period;
[f,x0]=p_correc_setup(funcs,ptini,ivarpar,method.point,...
    'remesh_flag',0,'previous',ptini,'output','J',pass_on{:});
ix=dde_ind_from_point(ptini,ivarpar);
[Jext,dum,ieq]=f(x0); %#ok<ASGLU>
vec=@(v)reshape(v,1,[]);
iev=vec(ix.profile(ip.dim+(1:ip.dim),:));
J=Jext([iev,ieq.extra],[iev,ix.parameter]);
%dde_psol_jac_res(funcs,pfoldini,orig_free_par,method.point);
[nullvecs,adj,sv]=dde_svdspaces_lr(J,options.nulldim,'nullspaces',false); %#ok<ASGLU>
template=point;
template.parameter=zeros(1,nint);
vpoint=dde_point_from_x(nullvecs,template,1:nint);
poev1ini=ptini;
%normv2=p_dot(vpoint,vpoint,'free_par_ind',orig_freepar,'period',true);
%normv=sqrt(normv2);
%vpoint=p_axpy(1/normv,vpoint,[]);
poev1ini.profile(ip.dim+(1:ip.dim),:)=vpoint(1).profile;
poev1ini.parameter(ip.beta)=vpoint(1).period;
poev1ini.parameter(ip.nullparind(:,2))=vpoint(1).parameter;
normv2=options.v_scal(poev1ini,poev1ini);
poev1ini.profile(ip.dim+1:end,:)=poev1ini.profile(ip.dim+1:end,:)/sqrt(normv2);
poev1ini.parameter([ip.beta,ip.nullparind(:,2)'])=...
    poev1ini.parameter([ip.beta,ip.nullparind(:,2)'])/sqrt(normv2);
end
