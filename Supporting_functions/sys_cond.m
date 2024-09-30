function [r,J]=sys_cond(p)
if strcmp(p.kind,'psol')
    r=p.profile(5,1);
    J=p_axpy(0,p,[]);
    J.profile(5,1)=1;
else
    r=[];
    J=p_axpy(0,p,[]);
    J=repmat(J,0,1);
end
end