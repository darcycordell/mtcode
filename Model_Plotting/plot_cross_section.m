function [x,z,rho] = plot_cross_section(x,z,rho)
%Function to plot a cross section of resistivity in log scale
%
% Usage: plot_cross_section(x,z,rho)
%
% x is the x vector to plot (km)
% z is the depth vector to plot (km)
% rho is the nx by nz resistivity matrix (Ohm m)

u = user_defaults;

%Create meshgrid of x and z values and plot as logarithmic resistivity
pcolor(x,z,rho); hold on; axis ij

xlabel('Distance Along Profile (km)')
ylabel('Depth Below Sealevel (km)')
colormap(u.cmap); caxis(u.colim); 
add_rho_colorbar(u);

set(gca,'DataAspectRatio',[u.ve 1 1]);
set(gca,'Layer','top')
set(gca,'Box','on');
%set(gca,'YScale','log')

if u.plot_contours
    contour(x,z,rho,u.contours,'-k','ShowText',u.contour_text);
end