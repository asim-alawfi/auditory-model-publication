function [r,J]=dde_pd_cond(p,xdim,t0)
Id=eye(xdim);
Z=zeros(xdim);
for i=length(t0):-1:1
    trafo=cat(1,[(1+cos(pi*t0(i)))*Id,(1-sin(pi*t0(i)))*Id],...
        [Z,Z]);
    [rc{i},Jc{i}]=dde_psol_lincond(p,xdim,'profile','trafo',trafo,'shift',[0,1],...
        'condprojint',t0(i)*[1,1],'condprojmat',[Id,Id],'stateproj',[Z,Id,Z;Z,Z,Id]);
end
r=cat(1,rc{:});
J=cat(1,Jc{:});
end
