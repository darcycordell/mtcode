function plot_rho_pha_map_det_by_site(d,ip)
%
% Function which plots apparent resistivity and phase for a single period in map view
% for the determinant apparent resistivity and phase, station by station.
% The code is based on plot_rho_pha_map_determinant and plot_misfit_rms_map
%
% Usage: plot_rho_pha_map_det_by_site(d,ip)
%
% "d" is an MT data structure
% "ip" is the index of the period to plot.

close all
u = user_defaults;
[d] = set_map_projection(d);
[L] = load_geoboundary_file_list;

if isnan(d.Zerr) %If all errors are NaN, then data is predicted data (from an inversion)
    datatype = 'Modelled';
else
    datatype = 'Observed';
end

% Calculate determinant impedance, apparent resistivity and phase
[detZ, detZerr] = calc_determinant(d.Z,d.Zerr);
[detrho, detpha, detrhoerr, detphaerr] = calc_rho_pha(detZ,detZerr,d.T);

set_figure_size(1);

%%
% Plot apparent resistivity calculated from determinant of Z

subplot(1,2,1) % rho

plot_geoboundaries(L);

for i = 1:d.ns % Plot each site rho as a colored circle at the site location in map view

    th = 0:pi/20:2*pi;
    r = abs(d.lim(4)-d.lim(3))/u.rmsscale; % This is the radius of the circle. It is 1/50 of the survey area
    xp = r * cos(th) + d.loc(i,2);
    yp = r * sin(th)*0.85 + d.loc(i,1);
    
    m_fill(xp,yp,log10(detrho(ip,:,i))); hold on
    shading flat;
end

m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
colormap(gca,u.cmap);
caxis(u.colim);
add_rho_colorbar(u);
set(gca,'Layer','top')
title('Determinant App. Res.')

%%
% Plot phase calculated from determinant of Z

subplot(1,2,2) % pha

plot_geoboundaries(L);

for i = 1:d.ns % Plot each site pha as a colored circle at the site location in map view

    th = 0:pi/20:2*pi;
    r = abs(d.lim(4)-d.lim(3))/u.rmsscale; % This is the radius of the circle. It is 1/50 of the survey area
    xp = r * cos(th) + d.loc(i,2);
    yp = r * sin(th)*0.85 + d.loc(i,1);

    m_fill(xp,yp,detpha(ip,:,i)); hold on
    shading flat;
end

m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
colormap(gca,flipud(u.cmap))
caxis(u.phalim)
hcb = colorbar;
hcb.Label.String = 'Phase (degrees)';
hcb.TickLabels = num2cell([0 15 30 45 60 75 90]);
hcb.Ticks = [0 15 30 45 60 75 90];
set(gca,'Layer','top')
title('Determinant Phase')

%%
% Add a figure title

annotation('textbox', [0 0.9 1 0.08], ...
    'String', [datatype,' App. Res. and Phase @ T = ',num2str(d.T(ip)),' s'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
end
