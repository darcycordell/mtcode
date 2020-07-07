function plot_rho_pha_map_determinant(d,ip)
%
% Function which plots interpolated apparent resistivity and phase for a single period in map view
% for the determinant apparent resistivity and phase. The code is very
% similar in structure to plot_rho_pha_map but with fewer subplots
%
% Usage: plot_rho_pha_map_determinant(d,ip)
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

%Set up interpolation grid
xgrid = d.lim(1):u.dx:d.lim(2);     
ygrid = d.lim(3):u.dy:d.lim(4); 
[X,Y] = meshgrid(xgrid,ygrid);
  
set_figure_size(1);

rhogrid = zeros(length(ygrid),length(xgrid),4);
phagrid = zeros(length(ygrid),length(xgrid),4);

[~,J] = find(squeeze(isnan(d.rho(ip,:,:)))); % find stations missing data
missing = unique(J);
[~,K] = find(squeeze(~isnan(d.rho(ip,:,:)))); % find stations missing data
notmissing = unique(K);

[detZ, detZerr] = calc_determinant(d.Z,d.Zerr);
[detrho, detpha, detrhoerr, detphaerr] = calc_rho_pha(detZ,detZerr,d.T);

% plot apparent resistivity calculated from determinant of Z

F = scatteredInterpolant(d.loc(notmissing,2),d.loc(notmissing,1), squeeze(log10(detrho(ip,:,notmissing))),'linear','none');
rhogrid = F(X,Y);

subplot(1,2,1) % rho av
m_pcolor(X,Y,rhogrid); shading flat; hold on
m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
m_plot(d.loc(notmissing,2),d.loc(notmissing,1),'k.','markersize',6);
plot_geoboundaries(L);

m_plot(d.loc(missing,2),d.loc(missing,1),'ko','markersize',2);

colormap(gca,u.cmap);

caxis(u.colim);

add_rho_colorbar(u);

set(gca,'Layer','top')

title(['Determinant App. Res.'])

%%
% plot phase calculated from determinant of Z

F = scatteredInterpolant(d.loc(notmissing,2),d.loc(notmissing,1), squeeze(detpha(ip,:,notmissing)),'linear','none');
phagrid = F(X,Y);

subplot(1,2,2) % pha av
m_pcolor(X,Y,phagrid); shading flat; hold on
m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
m_plot(d.loc(notmissing,2),d.loc(notmissing,1),'k.','markersize',6);
plot_geoboundaries(L);

m_plot(d.loc(missing,2),d.loc(missing,1),'ko','markersize',2);

colormap(gca,flipud(u.cmap))
caxis(u.phalim)
hcb = colorbar;
hcb.Label.String = 'Phase (degrees)';
hcb.TickLabels = num2cell([0 15 30 45 60 75 90]);
hcb.Ticks = [0 15 30 45 60 75 90];
set(gca,'Layer','top')

title(['Determinant Phase'])


%%
annotation('textbox', [0 0.9 1 0.08], ...
    'String', ['Interpolated ',datatype,' App. Res. and Phase @ f = ',num2str(d.f(ip)),' Hz'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
end
