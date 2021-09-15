function plot_slice(m,id,d)
%
% Function to plot a slice through a model in model coordinates
%
% Usage: plot_slice(m,id,d) OR plot_slice(m,id)
%
% "m" is the model structure
% "id" is the index of the slice to be plotted
% "d" is an OPTIONAL data structure to plot site locations on the slice
if id > m.nz
    warning('Layer index greater than number of layers! Plotting bottom layer.')
    id = m.nz;
end

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
%pcolor(m.X(xind,yind)/1000,m.Y(xind,yind)/1000,log10(m.A(xind,yind,id))); hold on; shading flat;
y = [m.y(yind); m.y(yind(end)+1)]/1000;
x = [m.x(xind); m.x(xind(end)+1)]/1000;
C = log10(m.A(xind,yind,id));
C = horzcat(C,C(:,end));
C = vertcat(C,C(end,:));
pcolor(y,x,C); hold on;%view(0,90)

if strcmp(u.gridlines,'off')
    shading flat
end

%If data exists, plot site locations
if exist('d','var')
    dind = d.z<=m.z(id);
    plot(d.y(dind)/1000,d.x(dind)/1000,'k.','markersize',12); hold on; axis equal
    plot_geoboundaries(L,d.origin,0)

end

%If the chosen slice is all above topo (i.e. all NaN), then it is not
%possible to plot contours
if u.plot_contours %Plot contours if plot_contour flag is true
    if all(isnan(C(:)))~=1 && length(unique(C(~isnan(C))))>1
        contour(m.X(xind,yind)/1000,m.Y(xind,yind)/1000,C(1:end-1,1:end-1),u.contours,'-k','ShowText',u.contour_text);
    end
end


colormap(u.cmap); caxis(u.colim);
add_rho_colorbar(u);

axis equal
axis(axlims)
set(gca,'Layer','top')
xlabel('Distance East-West (km)')
ylabel('Distance North-South (km)')
title(['Layer ',num2str(id),' | Depth = ',num2str(m.z(id)/1000),' to ',num2str(m.z(id+1)/1000),' km b.s.l.']);

end