function plot_slice_map(m,id,d)
%
% Function which plots a model slice in map view (lat-long) using m_map
%
% Usage: plot_slice_map(m,id,d) OR plot_slice_map(m,id)
%
% "m" is the model structure
% "id" is the index of the slice to plot
% "d" is the data structure and must be included in this function

if id > m.nz
    warning('Layer index greater than number of layers! Plotting bottom layer.')
    id = m.nz;
end

u = user_defaults;
[L] = load_geoboundary_file_list;
[d] = set_map_projection(d);

if ~isfield(m,'LON')
    error('Your model has not been converted to lat-long coordinates yet! See "link_model_data.m" for more info')
end


m_pcolor(m.LON,m.LAT,log10(m.A(:,:,id))); shading flat; hold on
m_grid('box','fancy','xlabeldir','end','tickdir','in');

m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
plot_geoboundaries(L);

if u.plot_contours %Plot contours as well if plot_contours flag is true
    if all(all(isnan(m.A(:,:,id))))~=1
        m_contour(m.LON,m.LAT,log10(m.A(:,:,id)),u.contours,'-k','ShowText',u.contour_text);
    end
end

colormap(u.cmap); caxis(u.colim);
add_rho_colorbar(u);
title(['Layer ',num2str(id),' | Depth = ',num2str(m.z(id)/1000),' to ',num2str(m.z(id+1)/1000),' km b.s.l.']);

end