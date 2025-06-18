function [r,J]=sys_cond_coll_fixperiod(point,period,ispar)
%% fix period to given parameter (visible to user fcns)
% implemented as user condition
%% 
if nargin<3
    ispar=true;
end
J=p_axpy(0,point,[]);
J.period=-1;
if ispar
    ind_period=period;
    period=point.parameter(ind_period);
    J.parameter(ind_period)=1;
end
r=period-point.period;
end
