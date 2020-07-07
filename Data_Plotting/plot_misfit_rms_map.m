function plot_misfit_rms_map(dobs,dpred,s)
% Function which plots rms misfit between two MT datasets in map view.
% Overall rms for each site is plotted as a colored circle at that site
% location. Color is the rms value.
%
% Usage: 
% plot_misfit_rms_map(dobs,dpred)
% plot_misfit_rms_map(dobs,dpred,s)
%
% "dobs" is observed MT data. Errors are taken from this data structure
% "dpred" is predicted MT data. Errors are ignored (or assumed to be NaN)
% "s" is an OPTIONAL input that is the stats structure output from
% detailed_statistics or detailed_staistics_2D. If s is not an input, the
% misfit stats are computed from impedances by default.
%
%%
u = user_defaults;
[L] = load_geoboundary_file_list;

%Get detailed statistics and rms by site
if ~exist('s','var')
    disp('No stats structure input. Calculating rms misfit from impedances')
    s = detailed_statistics(dobs,dpred);
end

for i = 1:dobs.ns %Plot each site rms as a colored circle at the site location in map view

    th = 0:pi/20:2*pi;
    r = abs(dobs.lim(4)-dobs.lim(3))/u.rmsscale; %This is the radius of the circle. It is 1/50 of the survey area
    xp = r * cos(th) + dobs.loc(i,2);
    yp = r * sin(th)*0.85+ dobs.loc(i,1);

    m_fill(xp,yp,s.rms_site(i)); hold on;
    shading flat;

end

plot_geoboundaries(L);
caxis(u.rmslim);
% %---------------------------------------------------
% data = rms_site;
% indexValue = 2.37;     % value for which to set a particular color
% topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
% indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
% bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
% % Calculate where proportionally indexValue lies between minimum and
% % maximum values
% largest = max(u.rmslim);
% smallest = min(u.rmslim);
% index = dobs.ns*abs(indexValue-smallest)/(largest-smallest);
% % Create color map ranging from bottom color to index color
% % Multipling number of points by 100 adds more resolution
% customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
%             linspace(bottomcolor(2),indexColor(2),100*index)',...
%             linspace(bottomcolor(3),indexColor(3),100*index)'];
% % Create color map ranging from index color to top color
% % Multipling number of points by 100 adds more resolution
% customCMap2 = [linspace(indexColor(1),topColor(1),100*(dobs.ns-index))',...
%             linspace(indexColor(2),topColor(2),100*(dobs.ns-index))',...
%             linspace(indexColor(3),topColor(3),100*(dobs.ns-index))'];
% customCMap = [customCMap1;customCMap2];  % Combine colormaps
% colormap(customCMap)
% %--------------------------------------------------
colormap(jet)
hcb = colorbar;
hcb.Label.String = 'Root Mean Square Misfit';
m_grid('box','fancy','tickdir','in','xlabeldir','end'); hold on
title(['Root Mean Square By Station. Overall RMS = ',num2str(s.rms)]);


