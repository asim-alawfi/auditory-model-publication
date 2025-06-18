function out=dde_setupoutput(fcnstruc,branch,success,outputfuncs)
if nargin>3 && outputfuncs
    out={fcnstruc,branch,success};
else
    out={branch,success,fcnstruc};
end
end
