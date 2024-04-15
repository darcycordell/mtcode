function M3_add_topo(hObject, ~, ~)
    
    H=guidata(hObject);
    
    if ~isfield(H,'m')
        warndlg('Make or load a model before adding topography!')
        return
    elseif isfield(H,'AAt')
        warndlg('Topography has already been added to this mesh!')
        return
    end
%%
    aircellthk = str2double(get(H.airspacing,'string'));
    def_res = str2double(get(H.res,'string'));

    maplim = floor([min(H.d.loc(:,2)) max(H.d.loc(:,2)) min(H.d.loc(:,1)) max(H.d.loc(:,1))]);
    
    topo_source_menu = menu('Get Topography From:','SRTM (.hgt) files (Default; No bathymetry)','Grid files (.grd) in current directory (Bathymetry)','Geotiff File','xyzt file','Elevation matfile in Current Directory','Station Elevations (interpolates between)');
    
    if topo_source_menu == 6 % use stations elevations
        %This can be done as a "last resort" if you have no other
        %topography data available. It can also be useful when doing
        %synthetics.
        if ~isfield(H,'dat') % check if data has been loaded
            error('Error: No data loaded! Load data first.');
        end       
        H.topo_lon = H.d.loc(:,2);
        H.topo_lat = H.d.loc(:,1);
        H.topo_z = -H.d.loc(:,3); % negative sign needed to convert to m a.s.l.
    elseif topo_source_menu == 5 % existing elevation.mat file  
        %MTplot outputs an elevation.mat file. This file will only contain
        %elevations within the array boundaries (i.e will not include
        %padding cells). This should also be used only as a last resort if
        %no other topography is available. 
        disp('Loading topography from elevation.mat file in current directory.')
        [file,~] = uigetfile('*.mat','Select elevation mat file');
        mat = load(file); % contains structure with latitude vector, longitude vector, elevation matrix
        elev_field = cell2mat(fieldnames(mat)); % get the topo structure regardless of structure name
        elev = mat.(elev_field);
%         elev = getfield(mat,cell2mat(fieldnames(mat))); % another way to do it
        H.topo_lon = elev.lon;
        H.topo_lat = elev.lat;
        H.topo_z = elev.z; % this is in m a.s.l.
    elseif topo_source_menu == 3 % Geotiff file
        %Global Geotiff elevation data can be downloaded from https://maps.ngdc.noaa.gov/viewers/bathymetry/
        % new link: https://www.ncei.noaa.gov/maps/bathymetry/
        %
        curdir = pwd;
        [file,path] = uigetfile({'*.tif','*.tif'},'Select GeoTiff file');
        cd(path)
        
        [A,R] = readgeoraster(file);
        
        H.topo_lat = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1));
        H.topo_lon = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2));
        
        H.topo_z = double(A(end:-1:1,:)); %this is the topo matrix
        
        cd(curdir);

    elseif topo_source_menu == 4
        %Load a lat/long/topo file from e.g. GeoMapApp which allows such
        %exports as a standard text file
        %
        curdir = pwd;
        [file,path] = uigetfile({'*'},'Select Topography xyz file');
        cd(path)

        A = load(file);

        H.topo_lat = unique(A(:,2),'stable');
        H.topo_lon = unique(A(:,1),'stable');
        H.topo_z = reshape(A(:,3),length(H.topo_lon),length(H.topo_lat))';

        %H.topo_lat = H.topo_lat+[diff(H.topo_lat); diff(H.topo_lat(end-1:end))];
        %H.topo_lon = H.topo_lon+[diff(H.topo_lon); diff(H.topo_lon(end-1:end))];

        cd(curdir);
        
    elseif topo_source_menu == 2 % existing .grd files 
        %If you have a GMT format grid file, this can be loaded directly. 
        %This could be an SRTM grid file or any other grd type file
        %If you want to include bathymetry in your model, normal SRTM data
        %is not sufficient. 
        %
        % For global datasets, one option is NOAA: 
        % https://maps.ngdc.noaa.gov/viewers/wcs-client/
        % New link: https://www.ncei.noaa.gov/maps/grid-extract/
        %   Select region you want
        %   Select Output Format as GMT NetCDF
        %   Download and put into your current directory
        
        curdir = pwd;
        [file,path] = uigetfile('*.grd','Select grid file');
        cd(path)
        [H.topo_lon,H.topo_lat,H.topo_z] = grdread2(file);
        
        cd(curdir);
        H.topo_z = double(H.topo_z); %this is a matrix
        
    else % need to use get_srtm to find topo data for area
        disp('Downloading topography data for mesh area.')
        temp = M3_get_srtm(maplim);
        H.topo_lon = temp.lon;
        H.topo_lat = temp.lat;
        H.topo_z = temp.z; % this is in m a.s.l. this is a matrix
    end
    
    if length(H.topo_lon)*length(H.topo_lat)>10^7
       
       size_menu = menu('Elevation matrix is very large and may be very slow and/or crash MATLAB. Downsample?','Yes','No');
       
       if size_menu
           
           %Find the optimal number of indices to skip in the lat and lon
           %vectors: sqrt(n*m/10^8)
            indskip = ceil(sqrt(length(H.topo_lon)*length(H.topo_lat)/10^7));
            
            H.topo_lon = H.topo_lon(1:indskip:end); H.topo_lat = H.topo_lat(1:indskip:end);
            H.topo_z = H.topo_z(1:indskip:end,1:indskip:end);
       end
           
        
    end
            
    % convert topography to mesh coords in km       
    long = H.d.loc(:,2);
    lat  = H.d.loc(:,1); 
    H.cent_lat  = (max(lat)+min(lat))/2;
    H.cent_long = (max(long)+min(long))/2;
    
    if size(H.topo_lon) == size(H.topo_lat) % these are vectors from get_srtm
        [H.topo_y,H.topo_x] = geo2utm(H.topo_lon,H.topo_lat,H.cent_long,H.cent_lat); % this is converted relative to center of stations
    else % if these are not the same dimensions, need to reorganize 
        [H.topo_y,~] = geo2utm(H.topo_lon,ones(size(H.topo_lon)).*H.cent_lat,H.cent_long,H.cent_lat);
        [~,H.topo_x] = geo2utm(ones(size(H.topo_lat)).*H.cent_long,H.topo_lat,H.cent_long,H.cent_lat);
    end
       
    H.topo_y = (H.topo_y./1000)-500; % these are topo coords not including entire mesh area
    H.topo_x = H.topo_x./1000;
          
    if topo_source_menu ~= 6 % use grids for SRMT, .grd, .tiff or .mat file where topo data are in a regular grid. if data is from station locations, keep as vectors
        [H.topo_y,H.topo_x] = meshgrid(H.topo_y,H.topo_x); % try putting into matrices early on
    end
    
    while 1 % check to see if mesh bounds are outside bounds of defined topo... this could be made more efficient  
        
        if H.mesh_rot ~= 0 % mesh rotated to align with strike
            c = cosd(H.mesh_rot);    s = sind(H.mesh_rot); % rotate mesh to angle to see what topo is needed
            R=[c s;-s c];
%             lim_rot = R*[[min(H.XX) min(H.XX) max(H.XX) max(H.XX)]; [min(H.YY) max(H.YY) max(H.YY) min(H.YY)]];
%             lim_chk = [ min(lim_rot(1,:)) < min(min(H.topo_x))  max(lim_rot(1,:)) > max(max(H.topo_x))  min(lim_rot(2,:)) < min(min(H.topo_y))  max(lim_rot(2,:)) > max(max(H.topo_y))];
            lim_rot = R*[[min(H.YY) min(H.YY) max(H.YY) max(H.YY)]; [min(H.XX) max(H.XX) max(H.XX) min(H.XX)]];
            lim_chk = [ min(lim_rot(1,:)) < min(min(H.topo_y))  max(lim_rot(1,:)) > max(max(H.topo_y))  min(lim_rot(2,:)) < min(min(H.topo_x))  max(lim_rot(2,:)) > max(max(H.topo_x))];
        else
%             lim_chk = [ min(H.XX) < min(min(H.topo_x))  max(H.XX) > max(max(H.topo_x))  min(H.YY) < min(min(H.topo_y))  max(H.YY) > max(max(H.topo_y))];
            lim_chk = [ min(H.YY) < min(min(H.topo_y))  max(H.YY) > max(max(H.topo_y))  min(H.XX) < min(min(H.topo_x))  max(H.XX) > max(max(H.topo_x))];
        end
        
        if lim_chk(1) == 1
            maplim(1) = maplim(1) - 1; % only add one degree of SRTM data at a time
            disp('The west edge of mesh is outside the topo limits')
        end
        if lim_chk(2) == 1
            maplim(2) = maplim(2) + 1;
            disp('The east edge of mesh is outside the topo limits')
        end
        if lim_chk(3) == 1
            maplim(3) = maplim(3) - 1;
            disp('The south edge of mesh is outside the topo limits')
        end
        if lim_chk(4) == 1
            maplim(4) = maplim(4) + 1;
            disp('The north edge of mesh is outside the topo limits')
        end
        if topo_source_menu == 1 % if just downloading topo data   
            if sum(lim_chk) ~= 0 % need to download additional topo grid since part of mesh is outside of topo limits          
                [new_elev] = M3_get_srtm(maplim);
                [Y,X] = meshgrid(new_elev.lon,new_elev.lat);
                lat = reshape(X, length(new_elev.lat)*length(new_elev.lon), 1);
                lon = reshape(Y, length(new_elev.lat)*length(new_elev.lon), 1);
                % get coords of the new topo 
                [tmp_y,tmp_x] = geo2utm(lon,lat,H.cent_long,H.cent_lat); 
                tmp_x = reshape(tmp_x,length(new_elev.lat),length(new_elev.lon));
                tmp_y = reshape(tmp_y,length(new_elev.lat),length(new_elev.lon));

                H.topo_lon = new_elev.lon;
                H.topo_lat = new_elev.lat;
                H.topo_y = (tmp_y./1000)-500; % these are now matrices
                H.topo_x = tmp_x./1000; 
                H.topo_z = new_elev.z;
            else
                break
            end
        else % otherwise just give a warning if topo does not cover the entire mesh
            if sum(lim_chk) ~= 0
                disp('Warning - Topo will not be interpolated onto the entire mesh!')
            end
            break % no automated way to get more data in grd format
        end % end downloading additional topo for mesh
    
    end % end while
            
    if H.mesh_rot ~= 0 % mesh rotated to align with strike. need to rotate topo as well
        c = cosd(-H.mesh_rot);    s = sind(-H.mesh_rot); % stations get rotated opposite direction of the data
        R=[c s;-s c];

        topo_rot = R*[H.topo_y(:) H.topo_x(:)]'; % 2020-05 this should always work now since topo_y and topo_x are matrices of the same size
        
        H.topo_y = reshape(topo_rot(1,:),size(H.topo_y)); % this is a rotated grid
        H.topo_x = reshape(topo_rot(2,:),size(H.topo_x));        
    end
    
    
    if topo_source_menu ~= 6 % if using station elevations, topo not in matrix form yet
        
        plot_menu = menu('Plot Grid Interpolation?','Yes','No');
        
        if plot_menu == 1
        
            figure % figure to show final topo region and mesh boundaries
            pcolor(H.topo_y,H.topo_x,H.topo_z)
            hold on; plot(H.d.y./1000,H.d.x./1000,'ro','markersize',4)
            shading flat; axis equal

    %         if H.mesh_rot ~= 0 % mesh rotated to align with strike
    %             c = cosd(-H.mesh_rot);    s = sind(-H.mesh_rot); % stations get rotated opposite direction of the data
    %             R=[c s;-s c];
    %             lim_rot = R*[[min(H.XX) min(H.XX) max(H.XX) max(H.XX)]; [min(H.YY) max(H.YY) max(H.YY) min(H.YY)]];
    %             plot([lim_rot(1,:) lim_rot(1,1)],[lim_rot(2,:) lim_rot(2,1)],'k-')
    %         else       
                rectangle('position',[min(H.YY) min(H.XX) max(H.YY)-min(H.YY) max(H.XX)-min(H.XX)])
    %         end

            set(gca,'dataaspectratio',[1 1 1])
            title(['UTM After merging - rotated to ',num2str(H.mesh_rot),' degrees'])

            figure(100)
            [LON, LAT] = meshgrid(H.topo_lon,H.topo_lat);
            pcolor(LON,LAT,H.topo_z); hold on
            plot(H.d.loc(:,2),H.d.loc(:,1),'or'); axis equal
            shading flat
            title('Topo in True Geographic Coordinates')
        
        end
    end
        
%     figure % figure to compare to original lat/long coords
%     pcolor(new_elev.lon,new_elev.lat,H.topo_z)
%     hold on; plot(H.coords(:,1),H.coords(:,2),'ro')
%     shading flat
% %         set(gca,'dataaspectratio',[1 1 1])
%     title('Lat/Long After merging')
             %%             
    % create interpolant function using grid from topo and the actual elevations
    [my,mx] = meshgrid(H.YY,H.XX); % create grid from mesh
    disp('Interpolating topo data onto mesh...')
    % if using SRMT, .grd, or .mat, topo will be in grids. if using station locations, topo is in vectors
    elev_grid = griddata(H.topo_y,H.topo_x,H.topo_z,my,mx);   
    disp('Done')
                        
    if isnan(sum(sum(elev_grid))) % if some elevations are not defined in SRTM data
        elev_grid = inpaint_nans(elev_grid); % could be better if done on original dataset, but faster this way
        disp('Missing elevations have been interpolated')
    end
    
    elev_grid(elev_grid > max(max(H.topo_z))) = max(max(H.topo_z));
    elev_grid(elev_grid < min(min(H.topo_z))) = min(min(H.topo_z)); % avoid extrapolation toward mesh edges- 
    % this might otherwise result in elevations much larger or smaller than the original maximum elevation  

    figure % plot topo with NaN removed and replaced with interpolated values
    pcolor(my,mx,elev_grid); colorbar
    hold on; plot(H.d.y./1000,H.d.x./1000,'ro','markersize',4)
    title('Mesh Elevations')
    xlabel('Distance (km)')
    ylabel('Distance (km)')
    
    xpad = str2double(get(H.moutcX,'string'));
    ypad = str2double(get(H.moutcY,'string'));
        
    % find difference between max and min elevation, to see number of air cells needed
    %When adding grd files with bathymetry, we want to ignore all the
    %elevations below sea level in this calculation
    elev_grid_orig = elev_grid;
    elev_grid(elev_grid<=0)=0;
    H.numair = max([1 floor((max(max(elev_grid)) - min(min(elev_grid(xpad:H.nx-xpad,ypad:H.ny-ypad)))) / aircellthk)]);
    aircells = ones(size(H.AA,1),size(H.AA,2),H.numair).*def_res;
%%
    H.AAt = cat(3,aircells,H.AA);
        
    z_air = -fliplr(cumsum(ones(1,H.numair).*aircellthk))./1000; % the depths of the air cells

    Z_new = cat(2,z_air,H.Z)-round(min(min(elev_grid(xpad:H.nx-xpad,ypad:H.ny-ypad))))/1000; % new vector of depths b.s.l. including top and bottom of mesh
    Z_new = reshape(Z_new,1,1,length(Z_new));
    
    elev = repmat(Z_new,size(H.AA,1),size(H.AA,2)); % create matrix of cell elevations, same size as AA
    
    elev_grid = -repmat(elev_grid_orig,1,1,size(elev,3))./1000; % these are the land elevations
    
    %Note: when adding bathymetry the coastline may not look very nice
    %because of the interpolation of elevations onto a coarse grid. It
    %might be blocky, with some "ocean" cells stranded inland.
    
    %Find the model indices above (i.e. less than) the elevation grid
    ind_above_elev = find(elev<elev_grid);
    %Find the model indices above sea level
%     ind_asl = find(elev<=0);
    ind_asl = find(elev<0); % May 2020 - changed the inequality since elev is the TOP of each layer, in meters b.s.l.
    %Find the modle indices below sea level
%     ind_bsl = find(elev>0);
    ind_bsl = find(elev>=0);
    
    %Air cells are above elevation grid and above sea level
    ind_air = intersect(ind_above_elev,ind_asl);
    
    %Ocean cells are above elevation grid but below sea level
    ind_ocean = intersect(ind_above_elev,ind_bsl);
    
    %Replace air cells with 10^17 and ocean cells with 0.3 Ohm m
    H.AAt(ind_air) = 10^17;
    H.AAt(ind_ocean) = str2double(get(H.ocean_res,'string'));
    
    %H.AAt(elev<elev_grid) = 10^17; %DC Edit: ModEM needs 10^17 for air cells
    % right now the x and y are reversed, need to decide what format to keep it in


    H.Z = squeeze(Z_new)'; % replace old Z vector, transpose to keep as row vector
    H.nz = length(H.Z);
    H.top = min(H.Z)*1000;

    %Determine elevation topography surface
    H.Zsurf = zeros(H.nx,H.ny);
    for i = 1:H.nx
        for j = 1:H.ny
            
            ind = find(elev(i,j,:)>elev_grid(i,j,:),1,'first');
            
            if isempty(ind)
                ind = 1;
            end
            
            if ind == H.nz
                ind = H.nz;
            end
            
            H.Zsurf(i,j) = ind;
            
        end
    end
        
    % 3DGrid sets top of model to 0, so all station Z in data file are postiive. Here we will keep in depth b.s.l.
%     H.z = -H.d.loc(:,3); % in meters % also does not seem necessary
%     H.z_new = NaN(size(H.z)); % not used?
    
    H.elev_grid = elev_grid; % save interpolated mesh elevations
    M3_update_model_plot(H)
    
%     [~,I] = min(H.AAt,[],3);
%     figure % figure for debugging
%     meshz(mx',my',H.Z(I))
%     set(gca,'zdir','reverse')
%     set(gca,'dataaspectratio',[1 1 1])
%     hold on
%     plot3(H.x,H.y,H.z_new,'kv','markerfacecolor','k')
       
    guidata(hObject, H);
end % end add_topo_Callback