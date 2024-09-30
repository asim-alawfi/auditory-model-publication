function ph=coll_phase(po,icomp,sg)
if nargin<3
    sg=-1;
end
if nargin<2
    icomp=2;
end
v=po.profile(icomp,:);
ind=find(diff(sign(v))*sg>0);
isnip=floor((ind-1)/(po.degree+1))+1;
ph=root_interp(po.tbp(:,isnip)',v((isnip-1)*(po.degree+1)+(1:po.degree+1)));
end