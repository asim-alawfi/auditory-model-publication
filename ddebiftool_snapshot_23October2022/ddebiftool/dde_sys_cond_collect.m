function [r,J]=dde_sys_cond_collect(funcs,p,pref,conds)
%% collect a cell array of extra conditions
% $Id$
%%
res=cell(1,0);
Jac=cell(1,0);
for i=length(conds):-1:1
    [res{i},Jac{i}]=conds{i}(p,pref);
    res{i}=res{i}(:)';
    Jac{i}=Jac{i}(:)';
end
[ru,Ju]=sys_cond_user(p,pref,funcs);
r=[ru(:)',res{:}]';
J=[Ju(:)',Jac{:}]';
end
function [r,J]=sys_cond_user(point,pref,funcs)
%% call user defined extra conditions and embed derivative into extended point
% for psol extensions only at the moment
% $Id$
%%
[r,J]=call_sys_cond(funcs,point,pref);
%% check if this is extended system (works for coll/psol)
if ~isfield(funcs,'get_comp') ||~isfield(funcs,'userfuncs')
    return
else
    r0=r;
    J0=J;
end
userpoint=funcs.get_comp(point,'solution');
userfuncs=funcs.userfuncs;
userref=funcs.get_comp(pref,'solution');
[r,userJ]=call_sys_cond(userfuncs,userpoint,userref);
J=dde_apply({'dde_',point.kind,'_extendblanks'},userJ,point);
r=[r0;r];
J=[J0;J];
end
