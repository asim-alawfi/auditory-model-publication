function newbranch=br_remove_points(branch,indxparameter,indx_threshold)

if branch.point(end).parameter(indxparameter) < branch.point(1).parameter(indxparameter)
    branch=br_rvers(branch);
end
newbranch=branch;
newbranch.point=branch.point(1);
%N = numel(branch.point);
%Startindx=1;
%i=Startindx;
for i=1:indx_threshold
    newbranch.point(i)=branch.point(i);
end
end