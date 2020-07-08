function plot_stations_map_srtm(d)
%
% Function which plots map of MT station locations along with SRTM
% topography. The code looks on the MATLAB path for the relevant SRTM files
% (*.hgt) and, if it cannot find them, it downloads them to the current
% directory. It is recommended that, if files are downloaded, they are
% copied to a folder on the users path so that they do not have to be
% downloaded again.
%
% Usage: plot_stations_map_srtm(d)
%
% "d" is MT data structure
%

%%
u = user_defaults;

% Output locations to file
str='station_list.dat';
fid1=fopen(str,'w+');
for is=1:d.ns
fprintf(fid1,'%s',d.site{is});
fprintf(fid1,'%9.3f %9.3f %9.3f \n',[d.loc(is,1),d.loc(is,2), d.loc(is,3)]);
end   
fclose(fid1);

% Plot station map with SRTM elevations
if strcmp(u.maplims,'default')
    %maplims is defined as [minlat maxlat minlon maxlon]
    u.maplims = [d.lim(3)-1 d.lim(4)+1 d.lim(1)-1 d.lim(2)+1];
end

dlat = round(4*abs(u.maplims(2)-u.maplims(1))/10)/4; % interval for latitude ticks and labels
dlon = round(4*abs(u.maplims(4)-u.maplims(3))/10)/4; % interval for longitude ticks and labels
prec = -2; % precision for axis labels, to the power of 10


srtm = download_srtm(u.maplims); %Download SRTM tiles if necessary. This function outputs an elevation.mat file to the current directory

srtm.z(srtm.z==0) = NaN;

[qx,qy] = meshgrid(srtm.lon,srtm.lat);

set_figure_size(1);
ax = axesm('MapProjection',u.projection,'grid','on', ...
    'frame', 'on','maplonlimit',[d.lim(1) d.lim(2)],'maplatlimit',[d.lim(3) d.lim(4)],...
    'mlabellocation',dlon,'plabellocation',dlat,'mlinelocation',dlon,'plinelocation',dlat,'mlabelround',prec,...
    'plabelround',prec);
mlabel;plabel;

% plot topography using geoshow
geoshow(qy,qx,srtm.z,'displaytype','texturemap'); hold on
tightmap

cb = colorbar;
ylabel(cb,'Elevation (m)') % this works in MATLAB 2013a, 2016b
[cmap] = demcmap(u.elev_colim);

%Set it up so that the ocean is plotted white
if ~isempty(find(isnan(srtm.z),1))%
    cmap(1,:) = 1;
    colormap(cmap)
    caxis([-80 u.elev_colim(2)])
else
    colormap(cmap)
end


[L] = load_geoboundary_file_list;
d = set_map_projection(d);
% plot stations
p = geopoint(d.loc(:,1),d.loc(:,2));
geoshow(p,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8); hold on
plot_geoboundaries_geoshow(L);
for i = 1:d.ns
    textm(d.loc(i,1),d.loc(i,2),strrep(d.site{i},'_','\_')); % correct underscores
end