function [branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries,varargin)
%% extend DDE-BIFTOOL branch
% function [c_branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries)
% INPUT:
%   funcs problem functions
%	branch initial branch (contains method and initial points)
%	max_tries maximum number of tries
%   optional (named):
%   'plotaxis' (default []): if plotting is on
%   (branch.method.continuation.plot>=1) then the user may specify an axis
%   into which to plot, default [] chooses gca
% OUTPUT:
%	branch extended branch
%	succ number of succesfull corrections
%	fail number of failed corrections
%	rjct number of failures which ended in rejecting a point

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: br_contn.m 369 2019-08-27 00:07:02Z jansieber $
%
%% introduce optional argument to set plotting axis
default={'plotaxis',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch=replace_branch_pars(branch,branch.parameter.free,pass_on);
tries=0;
fail=0;
rjct=0;
successful=1;
bound=0;
bound_tau=0;
stop=0;
stop_tau=0;
kontinue=1;
tp_del=funcs.tp_del;

l=length(branch.point);
trim=@(point)dde_trim_point(point,branch.point(1));
if l<=1
    error('br_contn:start','BR_CONTN: could not start branch, length=%d.',l);
end
l_usetangent=branch.method.continuation.use_tangent && isfield(branch.point(end),'nvec');

method=branch.method.point;
free_par=branch.parameter.free;
max_step=branch.parameter.max_step;
min_bound=branch.parameter.min_bound;
max_bound=branch.parameter.max_bound;

growth_factor=branch.method.continuation.steplength_growth_factor;

% prepare plotting:

if branch.method.continuation.plot>0
    if branch.method.continuation.plot>=1 && isempty(options.plotaxis)
        options.plotaxis=gca;
    end
end

while kontinue && (tries<=max_tries || bound || bound_tau)
    tries=tries+1;
    l=length(branch.point);
    if l==1
        error('br_contn:fail',['BR_CONTN: could not continue branch,',...
            '%d points, %d fails, %d rejected.'],tries-fail,fail,rjct);
    end
    branch.point(l)=add_tangent(funcs,branch.method,branch.point(l),branch.point(l-1),free_par);
    branch.point(l-1)=add_tangent(funcs,branch.method,branch.point(l-1),branch.point(l),free_par);
    prev_point=branch.point(l-1);
    last_point=branch.point(l);
    %% check boundaries
    if bound % if we are on the boundary, we should stop, unless we crossed more
        stop=1;
    end
    if bound_tau % if delay crossed zero, we should stop
        stop_tau=1;
        bound_tau=0;
    end
    if (stop && bound_tau==0) || stop_tau
        tries=tries-1;
        break
    end
    [bound,bound_fraction,bound_parameter]=check_parameter_bound(...
        last_point,prev_point,{min_bound,max_bound},[-1,1]);
    if tp_del~=0
        % check sign of delays
        [delay_nr,t_z]=p_tsgn(funcs,last_point);
        if ~isempty(delay_nr)
            bound_tau=1;
        end
    end
    
    if bound && bound_tau
        bound_tau=0; % we first treat case bound~=0
    end
    
    if (tries>max_tries) && bound==0 && bound_tau==0
        break
    end
    if bound % we already crossed a boundary
        bound_secant=p_axpy(0,last_point,[]);
        bound_secant.parameter(bound)=1;
    end
    
    % predict and determine steplength
    
    if l==2 || branch.method.continuation.prediction==1
        % linear prediction
        secant=p_axpy(-1,last_point,prev_point);
        dist=p_norm(secant);
        if successful % use extrapolation
            steplength=growth_factor*dist;
        else
            if l_usetangent
                steplength=steplength/2;
            else % use interpolation
                steplength=-dist/2;
            end
            % break if steplength falls below minimum (if set)
            if isfield(branch.method.continuation,'steplength_minimum') && ...
                    abs(steplength)<branch.method.continuation.steplength_minimum
                break
            end
        end
        if bound
            steplength=-bound_fraction*dist;
        end
        new_point=p_axpy(-steplength/dist,secant,last_point);
        % check for maximal steplengths
        fraction=1;
        for j=1:size(max_step,1)
            if max_step(j,1)==0
                dp=p_norm(p_axpy(-1,last_point,new_point));
            elseif max_step(j,1)==-1
                if isfield(new_point,'period')
                    lp=setfield(last_point,'period',new_point.period); %#ok<SFLD>
                else
                    lp=last_point;
                end
                dp=p_norm(p_axpy(-1,lp,new_point));
            else
                dp=abs(new_point.parameter(max_step(j,1))-last_point.parameter(max_step(j,1)));
            end
            if dp>max_step(j,2)
                f=max_step(j,2)/dp;
                if f<fraction
                    fraction=f;
                end
            end
        end
        if fraction<1
            steplength=steplength*fraction;
            new_point=p_axpy(-steplength/dist,secant,last_point);
        end
        if bound_tau %negative delay
            % new_point=(tau_p*last_point-tau_n*prev_point)/(tau_p-tau_n)
            tau_n=p_tau(funcs,last_point,delay_nr,t_z);
            tau_p=p_tau(funcs,prev_point,delay_nr,t_z);
            del_tau=tau_p-tau_n;
            pp=p_axpy(-tau_n/del_tau,prev_point,[]);
            new_point=p_axpy(tau_p/del_tau,last_point,pp);
        end
    else
        err=[branch.method.continuation.prediction];
        error('br_contn:pred',['BR_CONTN: only linear prediction',...
            'is currently implemented, prediction=%d.'],err);
    end
    %% choose reference point for correction r.h.s
    %if l==2
        ref_point=last_point;
    %else
    %    ref_point=prev_point;
    %end
    pred_point=new_point;
    %% plot
    br_online_plot(options.plotaxis,branch,new_point,last_point,...
        length(branch.point)+1,'g','pred');
    %% check if user has instructed to stop
    dostop=check_stops('predictor',branch.method.continuation,[ref_point,pred_point]);
    if dostop
        break
    end
    %% correct
    if bound
        disp('BR_CONTN warning: boundary hit.');
        new_point.parameter(bound)=bound_parameter;
        pred_point.parameter(bound)=bound_parameter;
        [new_point,success]=...
            p_correc(funcs,pred_point,free_par,bound_secant,method,tries+1,ref_point);
    elseif bound_tau
        s=strcat('BR_CONTN warning: delay number_',num2str(delay_nr),' becomes negative.');
        disp(s);
        [fcnz,mthz,frz]=dde_point_delay_zero_prep(funcs,method,free_par,delay_nr,t_z);
        [new_point,success]=p_correc(fcnz,pred_point,frz,[],mthz,0,ref_point);
    else
        if ~branch.method.continuation.steplength_condition
            secant=[];
        else % normalize secant
            secant=p_secant(secant,p_norm(pred_point));
        end
        stpcond=secant;
        if l_usetangent
            stpcond=last_point.nvec.tangent.vector;
        end
        [new_point,success]=...
            p_correc(funcs,pred_point,free_par,stpcond,method,tries+1,ref_point);
        
    end
    new_success=success;
    %JS: check for nan, inf
    pdiff=p_axpy(-1,last_point,new_point);
    %% compute angle between points
    if ~isfinite(p_norm(pdiff))
        new_success=false;
    end
    %% if minimal angle requested: call too large angle unsuccessful
    % take into account that after previously not successful step the new
    % point is between last and previous point
    if l_usetangent && new_success && ~bound && ~bound_tau && ...
            isfield(branch.method.continuation,'minimal_angle') 
        cosang=p_dot(pdiff,last_point.nvec.tangent.vector)/p_norm(pdiff);
        if ~(cosang>=branch.method.continuation.minimal_angle)
            new_success=false;
        end
    end
    if new_success
        % do some normalisations (interferes with general normalization
        % conditions for problems with symmetry, eg, pitchfork bifs)
        %new_point=p_normlz(new_point);
        % plot
        br_online_plot(options.plotaxis,branch,new_point,last_point,...
            length(branch.point)+1,...
            'b','                        cor');
    else
        fail=fail+1;
    end
    
    % keep new or throw away and maybe throw away last branch point too
    
    if new_success
        if bound || bound_tau
            branch.point(l)=trim(new_point);
        elseif successful || l_usetangent
            branch.point(l+1)=trim(new_point);
        else
            branch.point(l+1)=branch.point(l);
            branch.point(l)=trim(new_point);
        end
        %% check if user has instructed to stop
        dostop=check_stops('corrector',branch.method.continuation,...
            branch.point(max(1,end-1):end));
        if dostop
            break
        end
    elseif ~successful || bound || bound_tau
        bound=0;
        bound_tau=0;
        if branch.method.continuation.halt_before_reject==0 
            if~l_usetangent
                branch.point=branch.point(1:l-1);
            end
        else
            kontinue=0;
        end
        rjct=rjct+1;
    end
    
    successful=new_success;
    
end
succ=tries-fail;
end
%% check if stop function is installed
% by default stop functions are checked after correction. If structure is
% provided, it may have fields 'predictor' and 'corrector'. Their values
% can be cell arrays of functions accepting a point.
function dostop=check_stops(state,mth,p)
dostop=false;
if ~isfield(mth,'stops')
    return
end
stopfcn=mth.stops;
switch state
    case 'predictor'
        if isstruct(stopfcn) && isfield(stopfcn,'predictor')
            tests=stopfcn.predictor;
        else
            return
        end
    case 'corrector'
        if isstruct(stopfcn) && isfield(stopfcn,'corrector')
            tests=stopfcn.corrector;
        else
            tests=stopfcn;
        end
end
for i=1:length(tests)
    dostop=dostop || tests{i}(p);
    if dostop
        fprintf('BR_CONTN warning: boundary %d hit.\n',i);
        break
    end
end
end
%% check if parameter bound has been crossed
function [bd,bd_frac,bd_param]=check_parameter_bound(last_point,prev_point,bdlist,sg)
bd=0;
bd_frac=NaN;
bd_param=NaN;
for i=1:length(bdlist)
    bds=bdlist{i};
    for j=1:size(bds,1)
        if sg(i)*(last_point.parameter(bds(j,1))-bds(j,2))>0 % over maximum
            bd=bds(j,1);
            param=last_point.parameter(bds(j,1));
            bd_frac=(param-bds(j,2) ) / ...
                (param-prev_point.parameter(bds(j,1)));
            bd_param=bds(j,2);
            break
        end
    end
end
end
%% add tangent
function point=add_tangent(funcs,mth,point,prev,free_par)
if ~mth.continuation.use_tangent || ~isfield(point,'nvec')
    return
end
cnd=@(p)isfield(p.nvec,'tangent') && ...
    length(p.nvec.tangent.free_par)==length(free_par) && ...
    all(p.nvec.tangent.free_par==free_par);
if ~cnd(point)
    stpcond=p_tangent(funcs,mth.point,point,free_par);
    stpcond=p_axpy(1/p_norm(stpcond),stpcond,[]);
    point.nvec.tangent=struct('free_par',free_par,'vector',stpcond);
end
if cnd(prev)
    cosang=p_dot(point.nvec.tangent.vector,prev.nvec.tangent.vector);
else
    secant=p_axpy(-1,prev,point);
    cosang=p_dot(point.nvec.tangent.vector,secant)/p_norm(secant);
end
if cosang<0
    point.nvec.tangent.vector=p_axpy(-1,point.nvec.tangent.vector,[]);
end
end
%% check use of tangent
function usetan=use_tangent(mth,point)
usetan=mth.use_tangent && isfield(point,'nvec');
end