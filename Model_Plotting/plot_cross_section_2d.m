function [x,y,z,rho] = plot_cross_section_2d(m,d)

u = user_defaults;
%close all

if exist('d','var') & d.loc ~= 0
    figure(1);
    plot_site_map(d);
    print_figure(['vertical_profiles'],['Profile_Model_',d.niter,'_MAP']);
end

yind = m.npad(2)+1:m.ny-m.npad(2);
zind = 1:nearestpoint(u.zmax*1000,m.cz);

figure(3); hold on
xlabel('Distance Along Profile (km)')
ylabel('Depth Below Sealevel (km)')
colormap(u.cmap); caxis(u.colim); 
add_rho_colorbar(u);

set(gca,'DataAspectRatio',[u.ve 1 1]);
    
[Y, Z] = meshgrid(m.cy(yind)/1000,m.cz(zind)/1000);
pcolor(Y,Z,squeeze(log10(m.A(1,yind,zind)))'); shading flat; hold on; axis ij
if exist('d','var')
    plot(d.y/1000,d.z/1000,'vk','MarkerFaceColor','k')
end
axis([min(m.cy(yind))/1000 max(m.cy(yind))/1000 u.zmin u.zmax])
title(['Cross Section Through 2D Resistivity Profile']);
print_figure(['vertical_profile'],['Profile_Model_',m.niter]);



y = m.cy(yind)/1000;    
x = ones(length(y),1).*m.cx(1)/1000;
z = m.cz(zind)/1000;
rho = squeeze(m.A(1,yind,zind));

end

