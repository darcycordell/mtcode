function plot_objfun(objfun,outpath)
%==========================================================================

nit=size(objfun,1);
cobjfun={'chisq';'rms';'roughness_vs_tau';'roughness';'closeness';...
         'objective_function';'back to main'};
while 1
    choice=menu('Plot',cobjfun);
    if choice==7 return; end
    figure
    plot([0:nit-1],objfun(:,choice));
    title(char(cobjfun(choice)));
    eval(['print -dpsc ',outpath,'/',char(cobjfun(choice)),'.ps']);
end
end