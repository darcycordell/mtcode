function h = logerrorbar(X,Y,dY,linetype1,linetype2) 
%Function to plot error bars on log-log plot and avoid negative numbers
%
% Usage: h = logerrorbar(X,Y,dY,linetype1,linetype2)
%
%   X = x data
%   Y = y data
%   dY = y error
%   linetype1 = string for plotting linetype of data. Example: linetype = '-.r'
%   linetype2 = string for plotting the linetype of the errorbar lines (usually should be '-' with some color).
%%

h = loglog(X,Y,linetype1); hold on

lower = Y-dY;
lower(lower<=0) = 10^-40; 
plot([X X]',[lower Y+dY]',linetype2); hold on %Plot vertical error bar


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


end %END LOGERRORBAR