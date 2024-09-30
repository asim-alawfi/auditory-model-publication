function [v,sv,point]=dde_stst_critical_space(funcs,point,varargin)
default={'initcond',[],'eigenvalue',[],'kind','fold','output','space',...
    'tolerance',1e-6};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[A,taus]=dde_stst_linearize_rhs(funcs,point);
taus=taus(:).';
m=size(A,3);
om=ones(m,1);
xx=point.x(:,om);
n=size(xx,1);
lhs_mat=funcs.lhs_matrix(n);
if strcmp(options.kind,'fold')
    options.eigenvalue=0;
end
if ~isempty(options.eigenvalue)
    z=options.eigenvalue;
else
    [freq,point]=dde_stst_critical_freq(funcs,point,pass_on{:});
    if freq>0
        z=1i*freq;
    else
        z=0;
    end
end
D=dde_stst_ch_matrix(A,taus,z,'lhs_matrix',lhs_mat);
pini=feval(['dde_',options.kind,'_create'],'point',point,'omega',imag(z));
pini.v=0*pini.x;
D_ext=zeros(0,length(pini.v));
if ~isempty(options.initcond)
    if ~iscell(options.initcond)
        options.initcond={options.initcond};
    end
    Jac=cell(1,0);
    for i=length(options.initcond):-1:1
        [dum,Jac{i}]=options.initcond{i}(pini,pini); %#ok<ASGLU>
        Jac{i}=Jac{i}(:)';
    end
    J=[Jac{:}]';
    for i=length(J):-1:1
        D_ext(i,:)=J(i).v.';
    end
end
if strcmp(options.kind,'hopf')
    D=[real(D), -imag(D);
       imag(D),  real(D)];
    D_ext=[real(D_ext),imag(D_ext)];
end
A=cat(1,D,D_ext);
[U,S,V]=svd(A); %#ok<ASGLU>
sv=diag(S);
v=V(:,end);
if strcmp(options.kind,'hopf')
    v=v(1:n)+1i*v(n+1:2*n);
end
pini.v=v;
switch options.output
    case {'sv','singular_values'}
        [v,sv]=deal(sv,v);
    case {'dim','dimension'}
        sv=sum(sv<=options.tolerance);
        [v,sv]=deal(sv,v);
    case {'pini'}
        v=pini;
    case {'point'}
        [v,sv,point]=deal(point,v,sv);
end
end
