function [d] = set_map_projection(d)
%Function which sets the map projection information for m_map given the
%latitude and longitude coordinates of the MT sites
%
% Usage: [d] = set_map_projection(d)
%
% Inputs: d is a standard data structure
%
% Outputs: d is a standard data structure with map projection information
%        appended
%
% Currently this function is called by link_model_data.m

u = user_defaults;

%Mapping constants which are added to the "user" data structure.

if length(u.maplims)==4 %If xy limits were specified in user_defaults
    % Geographic map limits
    d.buffer=0;  
    d.border=d.buffer*111*0.5; % size of borders on maps/sections (in kilometers)
    minlon = u.maplims(3)-d.buffer;
    maxlon = u.maplims(4)+d.buffer;
    minlat = u.maplims(1)-d.buffer;
    maxlat = u.maplims(2)+d.buffer;
    
else %If no xy limits were specified then calculate the limits
    
    if length(d.loc(:,1))~=1
        d.buffer=sqrt((max(d.loc(:,2))-min(d.loc(:,2)))^2 + (max(d.loc(:,2))-min(d.loc(:,2)) )^2 )*.1; % buffer is 10% of length of profile/grid
    else
        d.buffer = 0.1;
    end
    
    d.border=d.buffer*111*0.5; % size of borders on maps/sections (in kilometers)
    minlon = min(d.loc(:,2))-d.buffer;
    maxlon = max(d.loc(:,2))+d.buffer;
    minlat = min(d.loc(:,1))-d.buffer;
    maxlat = max(d.loc(:,1))+d.buffer;
        
end

d.lim=[minlon maxlon minlat maxlat];
%d.lim = [-121, -114, 48, 54];
% set map projection and limits for map plotting (using m_map)
m_proj(u.projection,'long',[d.lim(1) d.lim(2)],'lat',[d.lim(3) d.lim(4)]); % define the projection, and borders of the maps
