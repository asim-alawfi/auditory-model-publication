%% SetupTorusBifurcation - Initialize continuation of torus or period doubling bifurcations
%%
function [trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which bifurcation was discovered
%
% optional inputs
%
% * |contpar| (integer default |[]|):  set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on
%    torus branch
% * |dir| (integer, default |[]|): which parameter to vary initially along
%    torus branch (trbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to trbranch.method.continuation,
% trbranch.method.point and trbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |trbranch|: bifurcation branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% <html>
% $Id: SetupTorusBifurcation.m 374 2019-09-14 14:02:58Z jansieber $
% </html>
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3,'nremove',1,...
    'biftype','torus','jacobian_nonsquare',false,'stop_1_2',false,...
    'postprocess','','align',[]};
[init_args,others]=dde_options_filter('TorusInit',varargin);
[options,pass_on]=dde_set_options(default,others,'pass_on');
branch.point=branch.point(ind);
%% If branch is psol bifurcation remove extra args
[funcs,branch,addremove]=PsolFromPsolbif(funcs,branch);
options.nremove=[options.nremove(:).',addremove(:).'];
%% initialize branch of torus bifurcations (trbranch) and pass on optional args
trbranch=replace_branch_pars(branch,options.contpar,pass_on);
trbranch.method.point.postprocess=options.postprocess;
if options.jacobian_nonsquare
    trbranch.method.point.preprocess='dde_jac2square_preprocess';
end
%% set up numbering and values of additional parameters
point=dde_psol_create(trbranch.point);
[trfuncs,extra_freepar,initfuncs]=set_torusfuncs(funcs,point,branch.method,...
    'biftype',options.biftype,pass_on{:});
%% required amendments of structures for extended system
free_par=[trbranch.parameter.free,extra_freepar];
[dum,is]=unique(free_par); %#ok<ASGLU>
trbranch.parameter.free=free_par(sort(is));
trbranch.method.point.extra_condition=1;
%% create initial guess for correction
trini0=TorusInit(initfuncs,point,trbranch.method,'nremove',options.nremove,init_args{:});
if ~isempty(options.align)
    trini0=dde_TorusBifurcation_postprocess(trini0,struct('previous',options.align));
end

[trbranch,suc]=correct_ini(trfuncs,trbranch,trini0,...
    options.dir,options.step,options.correc);
if strcmp(options.biftype,'torus') && options.stop_1_2
    trbranch.parameter.max_bound=cat(1,trbranch.parameter.max_bound,...
        [trfuncs.ip.omega,1]);
    trbranch.parameter.min_bound=cat(1,trbranch.parameter.min_bound,...
        [trfuncs.ip.omega,-1]);
end
end
