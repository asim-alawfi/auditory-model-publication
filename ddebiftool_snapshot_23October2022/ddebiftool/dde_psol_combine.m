function E=dde_psol_combine(point,S,shift)
ns=size(S,1);
nx=size(S,2);
nd=size(shift,2);
t0=point.mesh(:);
nt=size(point.mesh,2);
E=sparse(ns*nt,nx*nt);
for i=1:nd
    ts=t0+shift(i,1)/shift(i,2);
    wrap=double(ts<0|ts>1);
    wrap(ts<0)=floor(wrap(ts<0));
    wrap(ts>1)=floor(wrap(ts>1));
    ts(ts<0)=ts(ts<0)-wrap(ts<0);
    ts(ts>1)=ts(ts>1)-wrap(ts>1);
    E1=dde_coll_eva(point.profile,point.mesh,ts(:)',point.degree,...
        'output','matrix','kron',false);
    E1(wrap~=0,end)=E1(wrap~=0,end)+wrap(wrap~=0);
    E1(wrap~=0,  1)=E1(wrap~=0,  1)-wrap(wrap~=0);
    E=E+kron(speye(nt),S(:,:,i))*kron(E1,speye(nx));
end
end