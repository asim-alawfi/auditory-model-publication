function [r,J]=dde_extreme_cond(p,c,ind_d, ind_t)
%% provide condition for extrema of cx(t0)-d=0, with t0 a free parameter
% for periodic orbits
t0=p.parameter(ind_t);
d=p.parameter(ind_d);
[x,Jx]=dde_coll_eva(p.profile,p.mesh,mod(t0,1),p.degree,'kron',true);
[dx,dJx]=dde_coll_eva(p.profile,p.mesh,mod(t0,1),p.degree,'kron',true,'diff',1);
d2x=dde_coll_eva(p.profile,p.mesh,mod(t0,1),p.degree,'kron',true,'diff',2);
r=[c*x-d;c*dx];
J=repmat(p_axpy(0,p,[]),2,1);
J(1).profile=reshape(c*Jx,size(p.profile));
J(1).parameter(ind_d)=-1;
J(1).parameter(ind_t)=c*dx;
J(2).profile=reshape(c*dJx,size(p.profile));
J(2).parameter(ind_t)=c*d2x;
end

