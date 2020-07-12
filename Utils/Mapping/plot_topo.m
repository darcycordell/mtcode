function plot_topo(d,flag)
% Function which plots SRTM topography downloaded from download_srtm.m
%
% Usage: plot_topo(d,flag)
%
% Inputs: d is a standard data structure
%         flag = 1 plots topography using geoshow
%         flag = 2 plots topography using UTM coordinates in 3-D (using surf)
%         flag = 3 plots topography using m_map
%
% The topography file is specified in user_defaults. If none is specified, then
% the topography is downloaded automatically. This can be slow, so it is
% preferable to specify the file in user_defaults.
%

u = user_defaults;
if ~u.plot_topo
    return
end

if strcmp(u.topo_file,'none')
    if ~isfield(d,'lim')
        d = set_map_projection(d);
    end

    maplims = [d.lim(3)-1 d.lim(4)+1 d.lim(1)-1 d.lim(2)+1];

    srtm = download_srtm(maplims);
else
    load(u.topo_file);
end




[srtm.LON,srtm.LAT] = meshgrid(srtm.lon,srtm.lat);
%[x,y] = geo2utm(long,lat,cent_long,orig_lat)

if flag == 1
    geoshow(srtm.LAT,srtm.LON,-srtm.z/1000,'DisplayType','texturemap')
elseif flag == 2

    %Convert y model coordinates to longitude referenced from the center of
    %mesh origin in lat-long
    [srtm.Y,~]=geo2utm(srtm.LON, d.origin(1)*ones(length(srtm.lat),1),d.origin(2), d.origin(1));

    %Convert x model coordinates to latitude
    [~,srtm.X]=geo2utm(d.origin(2)*ones(1,length(srtm.lon)),srtm.LAT, d.origin(2), d.origin(1));
    srtm.Y = srtm.Y-500000;

    surf(srtm.Y/1000,srtm.X/1000,-srtm.z/1000); shading flat;

elseif flag == 3
    if isempty(get_projection)
        disp('Map projection not yet initialized. Defaults being set')
        m_proj('mercator','long',[-180 180],'lat',[-90 90]);
    end
    h = m_pcolor(srtm.LON,srtm.LAT,-srtm.z/1000);
    set(h,'EdgeColor','none');
else

end

cmap = flipud(demcmap(u.elev_colim/1000));

%Set it up so that the ocean is plotted white
if ~isempty(find(isnan(srtm.z),1))%
    cmap(1,:) = 1;
    colormap(cmap)
    caxis(sort(-[-80 u.elev_colim(2)])/1000)
else
    colormap(cmap)
    caxis(sort(-u.elev_colim)/1000)
end

cb = colorbar;
ylabel(cb,'Elevation (km)') % this works in MATLAB 2013a, 2016b
set(cb,'YDir','reverse');

alpha(u.topo_transparency)


