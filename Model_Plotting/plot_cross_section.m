function [x,z,rho] = plot_cross_section(x,z,rho)
%Function to plot a cross section of resistivity in log scale
%
% Usage: plot_cross_section(x,z,rho)
%
% x is the x vector to plot (km)
% z is the depth vector to plot (km)
% rho is the nx by nz resistivity matrix (Ohm m)

u = user_defaults;

xlabel('Distance Along Profile (km)')
ylabel('Depth Below Sealevel (km)')
colormap(u.cmap); caxis(u.colim); 
add_rho_colorbar(u);

set(gca,'DataAspectRatio',[u.ve 1 1]);

%Create meshgrid of x and z values and plot as logarithmic resistivity
[X, Z] = meshgrid(x,z);
pcolor(X,Z,squeeze(log10(rho))'); shading flat; hold on; axis ij
axis([min(x) max(x) u.zmin u.zmax])
set(gca,'Layer','top')
set(gca,'Box','on');
%set(gca,'YScale','log')

if u.plot_contours
    contour(X,Z,squeeze(log10(rho))',u.contours,'-k','ShowText',u.contour_text);
end