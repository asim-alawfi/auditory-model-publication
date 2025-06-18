function [e_sym,m_sym]=err_plot_half(pt,ax,ntfine)
if nargin<2
    ax=[];
end
if nargin<3
    ntfine=1001;
end
err=auto_eqd(pt);
tcoarse=pt.mesh(1:pt.degree:end);
nt=length(tcoarse);
err_int=interp1(tcoarse,err,'linear','pp');
m_int=interp1(tcoarse,(0:nt-1)/nt,'linear','pp');
e_max=err(end);
err_est=fnder(err_int);
tfine=linspace(pt.mesh(1),pt.mesh(end),ntfine);
t1=tfine(tfine<0.5);
t2=tfine(tfine>=0.5);
e1=fnval(err_est,t1)/e_max;
e2=fnval(err_est,t2)/e_max;
se1=fnval(err_int,t1);
se2=fnval(err_int,t2(1:end-1));
me1=fnval(m_int,t1);
me2=fnval(m_int,t2(1:end-1));
e_sym=max(abs(se2-se1-se2(1)))/e_max;
m_sym=max(abs(me2-me1-me2(1)));
if isempty(ax)
    return
end
ish=ishold(ax);
clr=lines();
p1=plot(ax,t1,e1,'o','color',clr(2,:),...
    'DisplayName','$\hat{e}(t)$ (scaled)');
hold(ax,'on');
p2=plot(ax,t2-0.5,e2,'k.','MarkerSize',8 ,...
    'DisplayName','$\hat{e}(t+0.5)$ (scaled)');
if ~ish
    hold(ax,'off');
end
mxe=max([e1,e2]);
mne=min([e1,e2]);
ylim([mne*0.9,mxe*1.2]);
xlabel(ax,'$t$','Interpreter','latex');
legend(ax,[p1,p2],'Interpreter','latex','EdgeColor',0.5*[1,1,1])
set(ax,'FontSize',16,'FontWeight','bold','FontName','Courier',...
    'LineWidth',2)
end