function branch=br_reorder(branch,indx)
if branch.point(end).parameter(indx) > branch.point(1).parameter(indx)
    branch=br_rvers(branch);
end
