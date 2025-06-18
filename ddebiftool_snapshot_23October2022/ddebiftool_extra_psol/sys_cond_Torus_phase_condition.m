function [r,J]=sys_cond_Torus_phase_condition(point,pref,dim)
%% ensure that real and imaginary part of Floquet mode are orthogonal
% to reference
inner_matrix=kron([0,0,0;0,0,-1;0,1,0],speye(dim));
[r,Worth]=dde_coll_profile_dot(point,pref,'inner_matrix',inner_matrix);
J=p_axpy(0,point,[]);
J.profile(:)=Worth*pref.profile(:);
% Jtemplate=p_axpy(0,point,[]);
% irgx=1:dim;
% irgu=dim+(1:dim);
% irgv=dim*2+(1:dim);
% irguv=[irgu,irgv];
% uvpoint=Jtemplate;
% uvpoint.profile(irguv,:)=point.profile(irguv,:);
% vupoint=pref;
% vupoint.profile(irgx,:)=0;
% vupoint.profile([irgu,irgv],:)=vupoint.profile([irgv,irgu],:);
% [r,Worth]=dde_coll_profile_dot(uvpoint,vupoint);
% J=Jtemplate;
% J.profile(:)=Worth*vupoint.profile(:);
end