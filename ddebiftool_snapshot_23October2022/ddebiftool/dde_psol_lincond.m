function [r,J]=dde_psol_lincond(point,varargin)
if length(varargin)==1
    S=varargin{1};
else
    S=dde_lincond_struct(varargin{:});
end
if ~strcmp(S.fieldname,'profile')
    S=dde_psol_lincond_correct(point,S);
end
%% linear condition for kind psol
% impose linear condition on solution
ns=size(S.condprojmat,1);
fn=S.fieldname;
if ~isfield(point,fn) || ns==0
    r=zeros(0,1);
    J=repmat(point,0,1);
    return
end
x=point.(fn);
[xdim,nt]=size(x);
[np,nxs]=size(S.stateproj);
if nxs<xdim
    S.stateproj=[S.stateproj,zeros(size(S.stateproj,1),xdim-nxs)];
end
op=ones(np,1);
frac=S.rotation(1)/S.rotation(2);
rot=exp(2*pi*1i*frac);
if 2*frac==round(2*frac)
    rot=sign(real(rot));
    matrot=diag(rot(op,1));
else
    op=op(1:end/2);
    rrot=diag(real(rot(op,1)));
    irot=diag(imag(rot(op,1)));
    matrot=[rrot,-irot; irot,rrot];
end
coeffs=cat(3,S.condprojmat*S.trafo*matrot*S.stateproj,...
    -S.condprojmat*S.stateproj);
Es=dde_psol_combine(point,coeffs,cat(1,S.shift,[0,1]));
ms=S.condprojint;
isint=find(ms(:,1)<ms(:,2));
ispt=find(ms(:,1)==ms(:,2));
Jprofsum=zeros(ns,size(ms,1),ns*nt);
sp=point;
sp.(fn)=zeros(ns,nt);
%% integral parts of profile measure
for i=1:length(isint)
    idv=kron(ones(1,nt),eye(ns));
    [dum,W]=dde_coll_profile_dot(sp,sp,'bd',ms(isint(i),:)); %#ok<ASGLU>
    Jint=idv*W;
    Jprofsum(:,isint(i),:)=reshape(Jint,ns,1,ns*nt);
end
%% discrete parts of profile measure
Ept=dde_coll_eva(sp.(fn),point.mesh,ms(ispt,1)',point.degree,...
    'output','matrix','kron',false,'sparse',false);
for i=1:length(ispt)
    Jpt=kron(Ept(i,:),eye(ns));
    Jprofsum(:,ispt(i),:)=reshape(Jpt,ns,1,ns*nt);
end
%% residual
Jprof=reshape(Jprofsum,ns*size(ms,1),ns*nt)*Es;
r=Jprof*reshape(point.(fn),[],1);
if ~isempty(S.res_parameters)
    par_sel=S.res_parameters~=0;
    par=zeros(size(r));
    par(par_sel)=point.parameter(S.res_parameters(par_sel));
    r=r-par;
end
%% Jacobian as point array
Jprof=reshape(full(Jprof),ns,size(ms,1),xdim,nt);
J=p_axpy(0,point,[]);
J=repmat(J,ns,size(ms,1));
for k=1:size(ms,1)
    for i=1:ns
        J(i,k).(fn)=reshape(Jprof(i,k,:,:),xdim,nt);
        if ~isempty(S.res_parameters)&& par_sel(i,k)
            J(i,k).parameter(S.res_parameters(i,k))=-1;
        end
    end
end
J=J(:);
end
