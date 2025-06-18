function [r,J]=call_sys_cond(funcs,p,pref)
if nargin<3
    pref=p;
end
if isfield(funcs,'sys_cond_reference') && funcs.sys_cond_reference
    [r,J]=funcs.sys_cond(p,pref);
else
    [r,J]=funcs.sys_cond(p);
end
end
