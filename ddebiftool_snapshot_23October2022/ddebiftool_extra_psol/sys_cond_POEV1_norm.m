function [vtvres,vtvJ]=sys_cond_POEV1_norm(point,dim,free_par_ind,varargin)
%% obtain condition that nullvector/eigenvector has length sqrt(res)
default={'res',1,'period',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
vpoint=p_axpy(0,point,[]);
vpoint.profile(dim+1:end,:)=point.profile(dim+1:end,:);
vpoint.parameter=point.parameter;
[vtv,vtvJ]=p_dot(vpoint,vpoint,'free_par_ind',free_par_ind,pass_on{:});
vtvJ=p_axpy(2,vtvJ,[]);
vtvres=vtv-options.res;
end
