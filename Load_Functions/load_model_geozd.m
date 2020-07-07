function [m,d] = load_model_geozd(name,elev)
% Function which loads some arbitrary model with columns of longitude, 
%latitude, elevation and model parameter. The format is as follows:
%
% Longitude Latitude Elevation (mbsl) Data (e.g. density, velocity)
%  -70.0       36.0     2500            4.3
%
% Usage: [m,d] = load_model_geozd(name,elev)
%
% Inputs: "name" is a filename string of the model to load
%         "elev" is an OPTIONAL argument to pin the model to a particular 
%               elevation in meters above sealevel.
%               The default is 0.
%
% Outputs: "m" is the standard model structure.
%       "d" is the standard data structure but in this case only contains
%       the central latitude and longitude reference and dummy station
%       locations for geo-referencing the model.
%
% Note: the elevation surface ("m.Z") is not calculated. Also, this code
% was primarily written to load the Laguna del Maule Seismic Model and has not
% been widely tested for other models.


% [name,modpath]=uigetfile({'*'},'Choose model to load'); if name == 0; return; end
% 
% cd(modpath)
mod = table2array(readtable(name));

if ~exist('elev','var')
    elev = 0;
end


m.lon = flipud(unique(mod(:,1)));
m.lat = flipud(unique(mod(:,2)));
m.cz = 1000*unique(mod(:,3))-elev;

m.A = reshape(mod(:,4),length(m.lat),length(m.lon),length(m.cz));

[m.LON, m.LAT] = meshgrid(m.lon,m.lat);

d.origin(1) = (max(m.lat)+min(m.lat))/2;
d.origin(2) = (max(m.lon)+min(m.lon))/2;

[m.cy,~] = geo2utm(m.lon,m.lat(round(length(m.lat)/2)),d.origin(2),d.origin(1));
m.cy = m.cy-500000;
[~,m.cx] = geo2utm(m.lon(round(length(m.lon)/2)),m.lat,d.origin(2),d.origin(1));

m.nx = length(m.cx);
m.ny = length(m.cy);
m.nz = length(m.cz);

m.npad = [0 0];

[m.X,m.Y] = meshgrid(m.cy,m.cx);

m.x = (m.cx(1:end-1)+m.cx(2:end))/2;
m.x = [m.cx(1)-m.x(1)+m.cx(1); m.x];
m.x = [m.x; m.cx(end)-m.x(end)+m.cx(end)];
m.dx = diff(m.x);

m.y = (m.cy(1:end-1)+m.cy(2:end))/2;
m.y = [m.cy(1)-m.y(1)+m.cy(1); m.y];
m.y = [m.y; m.cy(end)-m.y(end)+m.cy(end)];
m.dy = diff(m.y);

m.z = (m.cz(1:end-1)+m.cz(2:end))/2;
m.z = [m.cz(1)-m.z(1)+m.cz(1); m.z];
m.z = [m.z; m.cz(end)-m.z(end)+m.cz(end)];
m.dz = diff(m.z);

m.origin = [min(m.x) min(m.y) min(m.z)];
m.name = name;





