function [r,J]=sys_cond_psoluser(point,ind,get_comp,userfuncs,pref)
%% call user defined extra conditions and embed derivative into extended point
% for psol extensions only at the moment
% $Id: sys_cond_psoluser.m 369 2019-08-27 00:07:02Z jansieber $
%%
userpoint=get_comp(point,'solution');
userref=get_comp(pref,'solution');
[r,userJ]=call_sys_cond(userfuncs,userpoint,userref);
nuserpar=length(userpoint.parameter);
Jtemplate=p_axpy(0,point,[]);
J=repmat(Jtemplate,length(userJ),1);
%% append artificial parameters and components as zero in derivative
for i=1:length(userJ)
    J(i).parameter(1:nuserpar)=userJ(i).parameter;
    J(i).profile(1:ind.dim,:)=userJ(i).profile;
    J(i).period=userJ(i).period;
end
end