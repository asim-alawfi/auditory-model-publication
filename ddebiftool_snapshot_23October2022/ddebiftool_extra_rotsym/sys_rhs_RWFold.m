function y=sys_rhs_RWFold(x,p,ip,orig_rhs,orig_dirderi)
%% rhs of extended DDE for fold of relative equilibria
%
% x extended state (orig state and null vector)
% p user parameters
% ip.omrga, ip.rho rotation frequency and its derivative
% orig_rhs user r.h.s
% ip.dim original system dimension
dim=ip.dim;

x0=x(1:dim,:,:);
v=x(dim+1:end,:,:);
y0=orig_rhs(x0,p(1:ip.omega));
dp=0*p(1:ip.omega);
dp(ip.omega)=p(ip.rho);
%% add partial derivatives of all terms
y1=orig_dirderi(x0,p(1:ip.omega),v,dp);
y=[y0;y1];
end
