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
    xind = nearestpoint(u.xylims(1),m.cx/1000):nearestpoint(u.xylims(2),m.cx/1000);
    yind = nearestpoint(u.xylims(3),m.cy/1000):nearestpoint(u.xylims(4),m.cy/1000);
    
else %If the limits specified are not a vector of numbers with length 4 then just take non-padding cells
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
end

%PLOT MODEL
pcolor(m.X(xind,yind)/1000,m.Y(xind,yind)/1000,-m.Z(xind,yind)); shading flat; hold on

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

    