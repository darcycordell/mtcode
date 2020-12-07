function plot_model_topography(m,d)
%
% Function which plots a model's topography in pcolor view
%
% Usage: plot_model_topography(m,d)
%
% Inputs:
%   m is a standard model structure
%   d is an OPTIONAL data structure (used to plot station locations only)
%
% Note: If you want to plot data, then the model and data must be linked
%       Run [m,d] = link_model_data(m,d) to link them.
%


u = user_defaults;

[L] = load_geoboundary_file_list;

%Set plotting limits
if length(u.xylims)==4 %If xy limits were specified in user_defaults
    xind = nearestpoint(u.xylims(1),m.cx/1000,'previous'):nearestpoint(u.xylims(2),m.cx/1000,'next');
    yind = nearestpoint(u.xylims(3),m.cy/1000,'previous'):nearestpoint(u.xylims(4),m.cy/1000,'next');
    axlims = [u.xylims(3) u.xylims(4) u.xylims(1) u.xylims(2)];
    
else %If the limits specified are not a vector of numbers with length 4 then just take non-padding cells
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
    axlims = [sort([m.y(min(yind)) m.y(max(yind))]) sort([m.x(min(xind)) m.x(max(xind))])]/1000;
end

if any(isnan(xind))
    xind = 1:m.nx;
end

if any(isnan(yind))
    yind = 1:m.ny;
end

%PLOT MODEL
y = [m.y(yind); m.y(yind(end)+1)]/1000;
x = [m.x(xind); m.x(xind(end)+1)]/1000;
C = m.Z(xind,yind);
C = horzcat(C,C(:,end));
C = vertcat(C,C(end,:));
pcolor(y,x,C); hold on;

if strcmp(u.gridlines,'off')
    shading flat
end

%If data exists, plot site locations
if exist('d','var')
    
    plot(d.y/1000,d.x/1000,'k.','markersize',12); hold on; axis equal
    plot_geoboundaries(L,d.origin,0)

end

%If the chosen slice is all above topo (i.e. all NaN), then it is not
%possible to plot contours
contour_vec = -500:250:max(-m.Z(:));
if u.plot_contours %Plot contours if plot_contour flag is true
    contour(m.X(xind,yind)/1000,m.Y(xind,yind)/1000,-m.Z(xind,yind),contour_vec,'-k','ShowText',u.contour_text);
end

axis equal
axis([sort([m.cy(min(yind)) m.cy(max(yind))]) sort([m.cx(min(xind)) m.cx(max(xind))])]/1000)
set(gca,'Layer','top')
xlabel('Distance East-West (km)')
ylabel('Distance North-South (km)')
title(['Model Topography Surface']);

cb = colorbar;
ylabel(cb,'Elevation (m)') % this works in MATLAB 2013a, 2016b
[cmap] = demcmap(u.elev_colim);

colormap(cmap)

    