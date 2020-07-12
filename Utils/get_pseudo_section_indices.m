function [sidx,midpoint,azimuth] = get_pseudo_section_indices(d)
% 
% Usage: [sidx,midpoint,rot_ang] = get_pseudo_section_indices(d)
%
% Inputs: 
%       d:  data structure
% Outputs:
%       sidx: a list of station indices to pull sites from (e.g. to get the
%           names of the sites that are along the profile type "d.site{sidx}" )
%       midpoint: a [lon lat] pair to define the midpoint of the profile
%       azimuth: the profile azimuth (positive degrees is clockwise from north)
%
% updated May 2020
%
% Function has two main menus:
%
% 1) Choose a profile to plot
%    Option 1: Load an azimuth from user_defaults. The default of 0 degrees is an E-W profile.
%    Option 2: Draw a profile on the map
%    Option 3: Load profile endpoints from a text file. Format is two lines:
%               Lat1 Lon1
%               Lat2 Lon2
%
% 2) Choose stations to plot on the profile
%    Option 1: Load a text file from user_defaults. The default is to plot all stations
%    Option 2: Select stations on the map
%    Option 3: Load station indices from a text file. Format is one station
%    index per line:
%               1
%               2
%               3
%
% These stations along the profile are then used to plot a pseudo-section
% of some data parameter (e.g. rho, phase, phase tensor, etc.)
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)



profile_menu = menu('Choose a profile to plot','Load from user_defaults','Plot Along Line (Clicking)','Load Profile Endpoints from Text File');

u = user_defaults;
%%
if profile_menu == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load azimuth from user_defaults  
    
    midpoint = [d.origin(2) d.origin(1)]; % lon, lat 
    azimuth = u.profile_azimuth;
    lat = [d.origin(1) d.origin(1)]; % original before rotation by azimuth
    long = [d.lim(1) d.lim(2)];
    
    c = cosd(azimuth);      s = sind(azimuth);     R = [ c, -s ; s, c]'; % Tranpose to get a positive-clockwise rotation matrix

    loc = [long-midpoint(1);lat-midpoint(2)];
    loc_rot = R*loc + [midpoint' midpoint']; % [midpoint(1);midpoint(2)]; % addition needs same array size, version 2016a and older?
    long = loc_rot(1,:);
    lat = loc_rot(2,:);
    
    disp(['Profile azimuth of ',num2str(u.profile_azimuth),' degrees loaded from user_defaults'])
    disp('Profile midpoint is the center of your stations defined by d.origin')
    
    if isfield(d,'loc') % check for geographic coordinates before plotting map
        
        d = set_map_projection(d);
        [L] = load_geoboundary_file_list;

        figure(4)
        plot_geoboundaries_geoshow(L)
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(lat,long,'DisplayType','line','Color','r')
        axis(d.lim)
        xlabel('Longitude'); ylabel('Latitude')
        box on  
        
    end
    
elseif profile_menu == 2 || profile_menu == 3    
    
    if isfield(d,'loc') % check for geographic coordinates before plotting map
    
        d = set_map_projection(d);
        [L] = load_geoboundary_file_list;

        figure(4)
        plot_geoboundaries_geoshow(L)
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        axis(d.lim)
        xlabel('Longitude'); ylabel('Latitude')
        box on       

        if profile_menu == 2 % click on map
            title('Left-click the two endpoints of a profile')
            i = 1;
            long = [0,0]; lat = [0,0];
            while i<=2
            [long(i),lat(i)] = ginput(1);

            geoshow(lat(i),long(i),'DisplayType','point','Marker','*','MarkerEdgeColor','r','MarkerSize',15)
            i = i+1;
            end
            geoshow(lat,long,'DisplayType','line','Color','r')
            title(['Profile from (',num2str(long(1)),', ',num2str(lat(1)),') to (',num2str(long(2)),', ',num2str(lat(2)),')'])

            midpoint = [min(long) + (max(long)-min(long))/2, min(lat) + (max(lat)-min(lat))/2];
            azimuth = -atan2d(lat(2)-lat(1),long(2)-long(1)); % range [-180, 180] with positive degrees clockwise from North

        else % load endpoints from text file

            curdir = pwd;
            [profile_points_file, filepath]=uigetfile({'*.txt'},'Choose text file which contains profile endpoints'); 
            if profile_points_file == 0 
                return; 
            end
            cd(filepath)

            profile_loaded=load(profile_points_file);
            fclose all;
            cd(curdir)
            long =profile_loaded(:,2)'; % transpose just to keep dimensions consistent with clicking on map
            lat = profile_loaded(:,1)'; 

            geoshow(lat,long,'DisplayType','point','Marker','*','MarkerEdgeColor','r','MarkerSize',15)
            geoshow(lat,long,'DisplayType','line','Color','r','linewidth',1)
            title(['Profile from (',num2str(long(1)),', ',num2str(lat(1)),') to (',num2str(long(2)),', ',num2str(lat(2)),')'])                    

            midpoint = [min(long) + (max(long)-min(long))/2, min(lat) + (max(lat)-min(lat))/2];
            azimuth = -atan2d(lat(2)-lat(1),long(2)-long(1)); % range [-180, 180] with positive degrees clockwise from North
        end
        
    else % no locations exist in d structure
        disp('No geographic coordinates detected in data structure. Using default profile.')
        midpoint = [d.origin(2) d.origin(1)]; % lon, lat 
        azimuth = 0;
    end
    
end % end profile_menu
%%
stations_menu = menu('Choose stations to plot on the profile','Load from user_defaults','Select Stations on Map','Load Station Indices from Text File');

if stations_menu == 1 % Default all stations
    
    if strcmp(u.profile_stations,'default') % if set to default, use all stations
        sidx = 1:d.ns;
        disp('u.profile_stations set to default. Using all MT stations')
    else      
        if exist(u.profile_stations,'file') % check for file
            stations_loaded = load(u.profile_stations);
            fclose all;
            if any(mod(stations_loaded(:),1)) || any(stations_loaded==0) % check for non-integers or zero values
                error('The station indices file contained non-integer or zero values. Station indices must be non-zero integers!')
            end
            sidx = stations_loaded;
            disp(['Station indices loaded from: ',u.profile_stations])
        else 
            error(['Unable to find file ',u.profile_stations])
        end            
    end
    
elseif stations_menu == 2 % Select on map
    
    if ~isfield(d,'loc')
        error('Cannot select stations without coordinates in d structure')
    end
    
    if profile_menu == 1 % need to draw the map
        d = set_map_projection(d);
        [L] = load_geoboundary_file_list;

        figure(4)
        plot_geoboundaries_geoshow(L)
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(lat,long,'DisplayType','line','Color','r','linewidth',1)
        axis(d.lim)
        xlabel('Longitude'); ylabel('Latitude')
        box on
    end
    
    selected_ind = false(d.ns,1); % all stations are initially unselected
    while 1
        title('Click and Drag to Select Stations. (Select a station again to de-select it.)')
        rect = getrect(gca);
        in_ind = inpolygon(d.loc(:,2),d.loc(:,1),[rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1)],[rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2) rect(2)]);
        selected_ind = selected_ind + in_ind;
        selected_ind(selected_ind>1) = 0;
        selected_ind = logical(selected_ind);
        cla  
        geoshow(lat,long,'DisplayType','line','Color','r','linewidth',1)
        geoshow(d.loc(selected_ind,1),d.loc(selected_ind,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r')
        geoshow(d.loc(~selected_ind,1),d.loc(~selected_ind,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k')
        
        title('Left-click to continue selection. Right-click to Exit.')
        [~,~,button] = ginput(1);       
        if button == 3
            title(['Profile from (',num2str(long(1)),', ',num2str(lat(1)),') to (',num2str(long(2)),', ',num2str(lat(2)),')'])
            break
        end  
        
    end
    
    sidx = [1:d.ns]'.*selected_ind;
    sidx(sidx==0) = [];
    
    save_station_indices_menu = menu('Save station indices to text file?','Yes','No');
    
    if save_station_indices_menu == 1 % save to file
        curdir = pwd;
        file_default_name = ['station_indices_for_profile_', num2str(midpoint(1)), '_', num2str(midpoint(2)), '_azimuth_', num2str(azimuth), '.txt'];
        [output_station_indices_file, output_station_indices_path] = uiputfile('*.txt','Save station indices file as',file_default_name);       
        cd(output_station_indices_path)
        fid = fopen(output_station_indices_file,'w');
        fprintf(fid,'%d\n',sidx);
        fclose(fid);
        cd(curdir)
    end
    
elseif stations_menu == 3 % load station indices from text file
    
    curdir = pwd;
    [station_indices_file, filepath]=uigetfile({'*.txt'},'Choose text file which contains station indices'); 
    if station_indices_file == 0 
        return; 
    end
    cd(filepath)

    stations_loaded = load(station_indices_file);
    fclose all;
    cd(curdir)
    
    if any(mod(stations_loaded(:),1)) || any(stations_loaded==0) % check for non-integers or zero values
        error('The station indices file contained non-integer or zero values. Station indices must be non-zero integers!')
    end
    
    sidx = stations_loaded;
       
end % end stations_menu


end