function [freq,stst,ind]=dde_stst_critical_freq(funcs,stst,varargin)
default={'excludefreqs',[],'includehopf',false,...
    'method',getfield(df_mthod('stst'),'stability'),'evfield','l1','closest',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
excludefreqs=abs(options.excludefreqs(:).');
if isfield(stst,'omega') && ~isempty(stst.omega) && ~options.includehopf
    excludefreqs=[abs(stst.omega),excludefreqs];
end
if isempty(stst.stability)
    if isempty(funcs)
        error('dde_stst_critical_freq:arguments',...
            'dde_stst_critical _freq: eigenvectors not present in stst and r.h.s. not provided');
    end
    stst.stability=p_stabil(funcs,stst,options.method,pass_on{:});
end
evals=stst.stability.(options.evfield);
selimp=1:length(evals);
%% remove frequencies to be excluded
if ~isempty(excludefreqs)
    excludefreqs=[excludefreqs,-excludefreqs(excludefreqs>0)];
    ind=dde_match_complex(excludefreqs(:),imag(evals(selimp)));
    selimp(ind)=[];
end
selimp=selimp(find(imag(evals(selimp))>=0));
if isempty(selimp)
    error('dde_hopf_freq:eigenvalues',['dde_hopf_freq: ',...
        'no eigenvalues found.']);
end
%% select eigenvalue closest to imaginary axis
evals=evals(selimp);
if isempty(options.closest) || ~isnumeric(options.closest)
    [i1,i2]=min(abs(real(evals))); %#ok<ASGLU>
    freq=abs(imag(evals(i2)));
else
    [i1,i2]=min(abs(evals-options.closest)); %#ok<ASGLU>
    freq=abs(imag(evals(i2)));
end
ind=selimp(i2);
end
