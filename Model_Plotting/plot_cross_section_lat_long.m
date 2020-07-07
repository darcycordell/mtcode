function plot_cross_section_lat_long
% Function which takes a latitude, longitude, depth and resistivity values
% and plots them as a cross-section
%
% Usage: plot_cross_section_lat_long
%
% Inputs: None
%
% For example, S3D can output a matfile which contains lat,long,z,rho
% values for a cross-section. This matfile can then be loaded into this
% function and plotted. Useful for debugging to make sure that S3D is
% outputting the file correctly.
%
% This has not been thoroughly debugged.

    u = user_defaults;

    curdir = pwd;

    [filename, filepath]=uigetfile({'*.dat'},'Choose file which contains a cross section in longitude and latitude'); if filename == 0; return; end

    cd(filepath)
    section = load(filename);
    cd(curdir)

    z = unique(section(:,3));
    nz = length(z);
    np = length(section);
    ndist = np/nz; %<<<<<<<-----This is incorrect

    plot(section(:,1),section(:,2),'.k')

    res_plot = reshape(section(:,4),ndist,nz);

    res_plot(res_plot==10^9)=NaN;

    figure(1);
    subplot(1,2,1)
    pcolor(1:ndist,z,log10(res_plot')); axis ij; shading flat; hold on
    colormap(u.cmap); caxis(u.colim);
    add_rho_colorbar(u);
        
    xlabel(['Distance along Profile (km)']);
    ylabel('Depth (km)');
    set(gca,'Layer','top')

    subplot(1,2,2)
    plot(section(:,1),section(:,2),'.k'); axis equal

end
