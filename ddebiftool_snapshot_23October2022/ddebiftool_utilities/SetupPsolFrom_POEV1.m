%% Branch off at Branch point (point is either psol near BP or BP)
%%
function varargout=SetupPsolFrom_POEV1(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |psolbranch|: branch of |'psol'| psol solutions from which one wants to
% branch off
% * |ind|: index in |'point'| field of |psolbranch| which is close to
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'POfold.correction'|: whether POfold point should be corrected
% * |'contpar'|: index of continuation parameters if different from
% |'psolbranch'|
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * (|pfuncs|: modified rhs functions. This output is given only
% if requested by optional argument |'outputfuncs', true|.)
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
%
%%
default={'radius',0.01,'contpar',[],'corpar',[],'usercond',[],'outputfuncs',false};
[pofoldargs,args]=dde_options_filter('SetupPOEV1',varargin);
[options,pass_on]=dde_set_options(default,args,'pass_on');
% create branch per of periodic solutions branching off from an
% approximate BP point ind on a branch of periodic orbits (or POfold's)
if isempty(options.contpar)
    options.contpar=psolbranch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
[pofuncs,pfbr,suc]=SetupPOEV1(funcs,psolbranch,ind,'outputfuncs',true,...
    'contpar',options.contpar,'corpar',options.corpar,pofoldargs{:},pass_on{:});
if suc==0
    varargout=dde_setupoutput(funcs,[],suc,options.outputfuncs);
    return
end
pofold=pfbr.point(1);
deg_psol=pofuncs.get_comp(pofold,'solution');
dev=pofuncs.get_comp(pofold,'nullvector');
devnorm2=p_dot(dev,dev);
dev_psol=p_axpy(options.radius/sqrt(devnorm2),dev,deg_psol);
per=replace_branch_pars(setfield(psolbranch,'point',[deg_psol,dev_psol]),...
    options.contpar,pass_on); %#ok<SFLD>
funcs=dde_add_cond('SetupPsolFrom_POEV1',funcs,options);
[psol,suc]=p_correc(funcs,dev_psol,options.contpar,dev,per.method.point,1, ...
    dev_psol);
if suc==0
    varargout=dde_setupoutput(funcs,[],suc,options.outputfuncs);
    return;
end
per.point(2)=psol;
varargout=dde_setupoutput(funcs,per,suc,options.outputfuncs);
end
%%

