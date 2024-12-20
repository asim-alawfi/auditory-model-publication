%% dde23
function dy=audi_hopf_rhs(y,SD,SF,sig,par,in)
%SD:= delay variable
%parnames={'a','b','tau','tau_i','D','PR','df','TD','lambda','theta','c','m,...
% 'ph1',ph2};
%cind=[parnames;num2cell(1:length(parnames))];
%in=struct(cind{:});
% par:=column vectors contains parameters of the auditory model
%SF:=Sigmoid or Heaviside function with threshold activity \theta
% sig:= Sigmoid function
%% using Hopf normal form to mimic the periodic forcing
% we use TD as a delay to use it in the periodic force formla 
Sl1=SD(:,1);
Sl2=SD(:,2);
sin_indx=5;
d=par(in.c)*(1-nthroot(par(in.df),par(in.m))); % rescale 'd' according to the paper
iA=par(in.c)*sig(y(sin_indx))*sig(-Sl2(sin_indx))+d*sig(-y(sin_indx))*sig(Sl2(sin_indx));
iB=d*sig(y(sin_indx))*sig(-Sl2(sin_indx))+par(in.c)*sig(-y(sin_indx))*sig(Sl2(sin_indx));
%% Right hand side of auditory model ( system (1) in the paper)
% u_A=y(1); u_B=y(2); s_A=y(3); s_B=y(4);
%y(5),y(6) are  Hopf normal form varibles 
dy=[(-y(1)+SF(par(in.a)*y(2)-par(in.b)*Sl1(4)+iA))/par(in.tau);...
    (-y(2)+SF(par(in.a)*y(1)-par(in.b)*Sl1(3)+iB))/par(in.tau);...
    SF(y(1))*(1-y(3))/par(in.tau) - y(3)/par(in.tau_i);...
    SF(y(2))*(1-y(4))/par(in.tau) - y(4)/par(in.tau_i);...
    par(in.ph1)*y(5)-pi*par(in.PR)*y(6)- y(5)*(y(5)^2 + y(6)^2);...
    pi*par(in.PR)*y(5)+ par(in.ph1)*y(6)- y(6)*(y(5)^2 + y(6)^2)];
end
