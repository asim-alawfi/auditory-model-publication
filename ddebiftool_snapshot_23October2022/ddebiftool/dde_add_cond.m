function funcs=dde_add_cond(name,funcs,options,condname,output)
if nargin<4
    condname='usercond';
end
if nargin<5
    output=true;
end
%% add new sys_cond if given by user
if ~isempty(options.(condname))
    if output && ~options.outputfuncs
        error([name,':sys_cond'],[name,': ',...
            'new extra conditions change functions structure. Set\n',...
            '''outputfuncs'' to  true.']);
    end
    funcs.sys_cond=@(p,pref)dde_sys_cond_collect(funcs,p,pref,options.(condname));
    funcs.sys_cond_reference=true;
end
end
