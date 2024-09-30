function [trini,sv,uvpoint]=TorusInit(initfuncs,point,method,varargin)
%% crude initial guess for start of torus bifurcation from Floquet mode
%
%
ip=initfuncs.ip;
default={'v_scal',@(p,pref)sys_cond_Torus_norm(p,ip.dim,'res',0),'closest',[],...
    'nremove',1,'initmethod','eig','nulldim',2,'align',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.closest) || strcmp(options.initmethod,'eig')
    [eigval,eigprofile]=mult_crit(initfuncs.userfuncs,point,method.stability,...
        options.nremove,options.closest);
end
if isempty(options.closest)
    omega=atan2(imag(eigval),real(eigval))/pi;
else
    omega=atan2(imag(options.closest),real(options.closest))/pi;
end
trini=point;
trini.profile(ip.null,:)=0;
trini.profile(ip.xrg,:)=point.profile;
trini.parameter(1:ip.nuserpar)=trini.parameter;
trini.parameter(ip.omega)=omega;
trini.parameter(ip.period)=trini.period;
switch options.initmethod
    case 'eig'
        t=repmat(point.mesh,size(point.profile,1),1);
        % convert Floquet multiplier mode to Floquet exponent mode (periodic
        % function)
        eigprofile=eigprofile.*exp(-log(eigval).*t);
        upoint=p_axpy(0,point,[]);
        upoint.profile=reshape(real(eigprofile),size(point.profile));
        vpoint=upoint;
        vpoint.profile=reshape(imag(eigprofile),size(point.profile));
        utu=dde_coll_profile_dot(upoint,upoint);
        vtv=dde_coll_profile_dot(vpoint,vpoint);
        utv=dde_coll_profile_dot(upoint,vpoint);
        r=1/sqrt(utu+vtv);
        gamma=atan2(2*utv,vtv-utu)/2;
        qr=r*(upoint.profile*cos(gamma)-vpoint.profile*sin(gamma));
        qi=r*(upoint.profile*sin(gamma)+vpoint.profile*cos(gamma));
        trini.profile(ip.null,:)=[qr;qi];
        uvpoint=initfuncs.get_comp(trini,'eigenvector');
        sv=NaN(2,1);
    case 'svd'
        trini.parameter(ip.omega)=omega;
        mth=setfield(method.point,'preprocess',''); %#ok<SFLD>
        [f,x0]=p_correc_setup(initfuncs,trini,[],mth,...
            'remesh_flag',0,'previous',trini,'output','J',pass_on{:});
        [Jext,dum,ieq]=f(x0); %#ok<ASGLU>
        vec=@(v)reshape(v,1,[]);
        ix=dde_ind_from_point(trini,[]);
        iev=vec(ix.profile(ip.dim+(1:2*ip.dim),:));
        J=Jext([iev,ieq.extra],iev);
        [nullvec,adj,sv]=dde_svdspaces_lr(J,options.nulldim,'nullspaces',false); %#ok<ASGLU>
        trini.profile(iev)=nullvec(:,1);
        normv2=options.v_scal(trini,trini);
        trini.profile(iev)=trini.profile(iev)/sqrt(normv2);
        uvpoint=initfuncs.get_comp(repmat(trini,1,options.nulldim),'eigenvector');
        for i=size(nullvec,2):-1:1
            uvpoint(i).profile=reshape(nullvec(:,i),size(uvpoint(i).profile));
        end
end
end
