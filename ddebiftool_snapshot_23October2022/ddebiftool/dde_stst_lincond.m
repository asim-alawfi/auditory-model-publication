function [r,J]=dde_stst_lincond(point,varargin)
if length(varargin)==1
    S=varargin{1};
else
    S=dde_lincond_struct(varargin{:});
end
ns=size(S.condprojmat,1);
fn=S.fieldname;
if ~isfield(point,fn) || ns==0
    r=zeros(0,1);
    J=repmat(point,0,1);
    return
end
x=point.(fn);
[np,nxs]=size(S.stateproj);
if nxs<size(x,1)
    S.stateproj=[S.stateproj,zeros(size(S.stateproj,1),size(x,1)-nxs)];
end
frac=S.rotation(1)/S.rotation(2);
rot=exp(2*pi*1i*frac);
if 2*frac==round(2*frac)
    rot=sign(real(rot));
end
matrot=diag(rot(ones(np,1),1));
matsym=S.condprojmat*(S.trafo*matrot-eye(size(matrot)))*S.stateproj;
r=matsym*x;
if ~isempty(S.res_parameters)
    par_sel=S.res_parameters~=0;
    par=zeros(ns,1);
    par(par_sel)=point.parameter(S.res_parameters(par_sel));
    r=r-par;
end
J0=p_axpy(0,point,[]);
J=repmat(J0,ns,1);
for i=1:ns
    J(i).(fn)=reshape(matsym(i,:)',size(x,1),1);
    if ~isempty(S.res_parameters)&& par_sel(i)
        J(i).parameter(S.res_parameters(i))=-1;
    end
end
ind=dde_ind_from_point(point,1:length(point.parameter));
ix=ind.(fn);
if isstruct(ix)
    r=[real(r);imag(r)];
    J(ns+(1:ns))=J0;
    for i=ns:-1:1
        J(i+ns).(fn)=reshape(1i*J(i).(fn),size(x,1),1);
    end
end
end
