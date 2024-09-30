function r=root_interp(tvals,xvals,restrict)
t=2*(tvals-tvals(1))/(tvals(end)-tvals(1))-1;
d=length(t)-1;
t=t(:)';
Aint=[ones(1,d+1);t;NaN(d-1,d+1)];
for i=3:d+1
    Aint(i,:)=2*t.*Aint(i-1,:)-Aint(i-2,:);
end
ccheb=xvals/Aint;
c=ccheb(:,1:end-1)./repmat(ccheb(:,end),1,d);
cm0=diag([ones(d-2,1);2],-1)+diag(ones(d-1,1),1);
nx=size(xvals,1);
r=cell(nx,1);
if nargin<3
    restrict=2;
end
for i=1:nx
    cm=cm0;
    cm(1,:)=cm(1,:)-c(end:-1:1);
    r{i}=eig(cm/2);
    if restrict>0
        r{i}=r{i}(imag(r{i})==0);
    end
    if restrict>1
        r{i}=r{i}(r{i}<=1&r{i}>=-1);
    end
    r{i}=(r{i}+1)*0.5*(tvals(end)-tvals(1))+tvals(1);
end
if nx==1
    r=r{1};
end
end