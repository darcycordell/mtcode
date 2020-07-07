function plot_misfit_rms_by_site(dobs,dpred,is,s)
% Function which plots rms misfit for two MT datasets for all sites as well
% as a highlighted site specified by "is".
% Site rms is plotted as a function of site number
%
% Usage: plot_misfit_rms_by_site(dobs,dpred) OR plot_rms_by_site(dobs,dpred,is)
%
% "dobs" is observed MT data structure
% "dpred" is predicted MT data structure (no longer required since s is an
% input)
% "s" is a stats structure containing data misfit
% "is" is an OPTIONAL input to highlight a specific site index
%
%
%Note: This function is usually only called by other functions (see plot_data.m)
% It is usually placed in the lower right subplot of the currently open
% figure (2x2 or 2x3 subplots of data plots)
%



%Get detailed statistics and rms by site
% [all_stats,rms_site,~,~,~,~,~] = detailed_statistics(dobs,dpred);

if length(findobj(gcf,'type','axes')) > 4 %If diagonals exist plot in the 6th subplot position
    subplot(2,3,6)
elseif length(findobj(gcf,'type','axes')) == 3 || length(findobj(gcf,'type','axes')) == 4 %If diagonals do not exist, plot in the 4th subplot location
    subplot(2,2,4)
else %Otherwise just plot normally
    close all
    subplot(1,1,1)
end

if ~exist('s','var')
    s = detailed_statistics(dobs,dpred);
end

% title('RMS Error Plot'); hold on
plot(s.rms_site,'*k'); hold on
if exist('is','var')
    plot(is,s.rms_site(is),'*r','MarkerSize',12);
end
plot([0,dobs.ns+1],[s.rms, s.rms],'-r')
xlabel('Site Number')
ylabel('Root Mean Square Misfit')
axis([0 dobs.ns+1 0 max(s.rms_site)+1])
text(((max(dobs.ns)/2)+2),0.5,num2str(s.rms,'%.2f'),'Fontsize',10);
text(2,0.5,'Overall R.M.S. misfit:','Fontsize',10);
box on

if exist('is','var')
    annotation('textbox', [0 0.9 1 0.08], ...
        'String', strrep(['Compare Responses at Site ',num2str(is,'%03.0f'),': ',dobs.site{is},', R.M.S. = ',num2str(s.rms_site(is),'%.2f')],'_','\_'), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
end

                 

end