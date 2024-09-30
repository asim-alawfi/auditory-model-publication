%% Initialize continuation of equilibrium bifurcation (fold or Hopf)
function varargout=SetupStstBifurcation(funcs,branch,ind,kind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which bifurcation was discovered
% * |ind|: index of approximate bifurcation point
% * |type|: tpye of bifurcation, string: at the moment either |'fold'| or |'hopf'|
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
%  (replacing free pars of branch)
% * |'correc'| (logical, default true): apply |p_correc| to first points on
% bifurcation branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially
% along bifurcation branch (bifbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
%
% All other named arguments are passed on to fields of |bifbranch|
%% Outputs
% 
% * |bifbranch|: branch of bifurcation points with first point (or two points)
% * |suc|: flag whether corection was successful
% * |funcs|: same as |funcs|, unless branch point
% * (if option |'outputfuncs'| is |true|, first output is |funcs|)
%
% Parameter limits for bifbranch etc are inherited from branch, unless overridden by
% optional input arguments.
%
% $Id: SetupStstBifurcation.m 309 2018-10-28 19:02:42Z jansieber $
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3,...
    'usercond',[],'outputfuncs',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% initialize branch of bifurcations (bifbranch)
bifbranch=branch;
bifbranch=replace_branch_pars(bifbranch,options.contpar,[{'branchtype',kind},pass_on]);
point=branch.point(ind);
if ~isfield(point,'stability') || isempty(point.stability)
    point.stability=p_stabil(funcs,point,branch.method.stability);
end
%% create initial guess for correction
pini0=feval(['p_to',kind],funcs,point,'method',bifbranch.method,pass_on{:});
if options.outputfuncs
    pini0.nvec=[];
    pini0.nmfm=[];
end
%% add new sys_cond if given by user
funcs=dde_add_cond('SetupStstBifurcation',funcs,options);
%% correct and add 2nd point if desired
[bifbranch,suc]=correct_ini(funcs,bifbranch,pini0,...
    options.dir,options.step,options.correc);
varargout=dde_setupoutput(funcs,bifbranch,suc,options.outputfuncs);
end
