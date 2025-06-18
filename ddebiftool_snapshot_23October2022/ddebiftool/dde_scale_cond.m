function [r,J]=dde_scale_cond(pt,pref,cond,varargin)
default={'res',1,'matrix',1,'ref',false};
options=dde_set_options(default,varargin,'pass_on');
[r0,J0]=cond(pt,pref);
if ~isempty(r0)
    if options.ref
        [rref,Jref]=cond(pref,pref);
        fac=1;
    else
        rref=r0;
        Jref=J0;
        fac=2;
    end
    prefac=rref'*options.matrix;
    r=prefac*r0-options.res;
    J=p_axpy(0,Jref(1),[]);
    for i=1:numel(J0)
        J=p_axpy(fac*prefac(i),Jref(i),J);
    end
else
    r=r0;
    J=J0;
end
end
