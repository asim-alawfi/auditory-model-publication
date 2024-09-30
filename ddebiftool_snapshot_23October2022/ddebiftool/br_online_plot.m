function br_online_plot(ax,branch,new_point,last_point,step,clr,txt)
%% Online plotting during continuation
% this function is called during |br_contn|.
%
% *Inputs*
%
% * |ax|: axis object where plot is updated
% * |branch|: branch structure maintained during continuation. Field
% |'method.continuation'| contains fields relevant for plotting such as
% |plot_measure|.
% * |new_point|: predicted or corrected point (end point of line thatis
% drawn)
% * |last_point|: previous/reference point (starting point of line)
% * |step|: integer, counting the step (used in title for user information)
% * |clr|: color used for line (hard coded in |br_contn| at the moment),
% different colors indicate predictor step and correced step.
% * |txt|: replacement text used for printout if
% method.continuation.plot_progress>0 but <1 (hard coded as |'pred'| and |'
% cor'| in |br_contn| for now).
%%
mth=branch.method.continuation;
if mth.plot<=0
    return
end
m_isfun=false;
if isempty(mth.plot_measure)
    [x_m,y_m]=df_measr(0,branch);
elseif iscell(mth.plot_measure)
    x_m=mth.plot_measure{1};
    y_m=mth.plot_measure{2};
    m_isfun=true;
else
    x_m=mth.plot_measure.x;
    y_m=mth.plot_measure.y;
end
if ~isempty(x_m) && isstruct(x_m)
    x1=p_measur(last_point,x_m);
    x2=p_measur(new_point,x_m);
elseif m_isfun
    x1=feval(x_m,last_point);
    x2=feval(x_m,new_point);
else
    x1=step;
    x2=step+1;
end
if ~isempty(y_m) && isstruct(y_m)
    y1=p_measur(last_point,y_m);
    y2=p_measur(new_point,y_m);
elseif m_isfun
    y1=feval(y_m,last_point);
    y2=feval(y_m,new_point);
else
    y1=step;
    y2=step+1;
end
if mth.plot>=1
    try
        ish=ishold(ax);
        hold(ax,'on');
        if ~dde_isoctave()
            args={'HandleVisibility','off'};
        else
            args={};
        end
        plot(ax,[x1 x2],[y1 y2],clr,args{:});
        plot(ax,x2,y2,[clr,'.'],args{:});
        steplength=p_norm(p_axpy(-1,last_point,new_point));
        title(ax,sprintf('step %d, steplength=%5.3e',step,abs(steplength)));
        if ~ish
            hold(ax,'off');
        end
    catch ME
        warning('br_contn:online_plot',...
            'br_contn: online plotting requested but plot command failed.\n%s',ME.message);
    end
    if mth.plot_progress
        drawnow;
    end
elseif mth.plot_progress
    fprintf('%s: (%g, %g)\n',txt,x2,y2);
end
end
