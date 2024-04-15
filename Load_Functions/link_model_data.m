function [m,d] = link_model_data(m,d)
% Function to link up a dataset stucture and a model structure by
% converting the model coordinates to latitude and longitude and the data
% site coordinates to x-y UTM if they don't yet exist. This function also
% does corrections to site locations if they are not consistent between
% data and model (likely due to old 85000 error in M3DET). This function
% also deals with meshes which have been rotated.
%
% Usage: [m,d] = link_model_data(m,d)
%
% "m" is model structure used to reference the x,y UTM locations of sites
% if they don't yet exist
% "d" is the data structure used to reference the model to lat-long.
%
% Variables m.lon, m.lat, d.x, d.y, d.z, d.idx, and d.idy are all appended
% to their respective structures if they don't yet exist.
%
% Note: If your model cells are very small in the x-y direction (e.g. <50
% m), there may be some problems in this code because of the precision of
% the lat/long coordinates. ModEM saves output data files with 3 digits of
% precision in lat/long which may not be sufficient for very small grid
% cells.
%
% To do:
% -Currently there is no similar function to map 2D NLCG profiles to lat-long
% space because the 2D NLCG model files have no lat-long reference point.
%
if d.origin == 0
    disp('Warning: Your reference latitude and longitude is (0,0). Are you sure this is correct?')
end

if isfield(d,'x')
    %-------------------------------------------------------------------------
    %CORRECT BAD SITE LOCATIONS IN X AND Y

    %Sometimes in previous code versions, there was a dubious
    %85000 added to the site locations (in M3DET_vxx).
    %Enci noticed this and added a fix. I have modified the fix slightly so
    %that it does not change the origin but only changes the site x and y locations 
    [long,lat] = utm2geo(d.y+500000,d.x,d.origin(2),d.origin(1));

    differenceinlong=(d.loc(:,2)-long); dlong=mean(differenceinlong);
    differenceinlat=(d.loc(:,1)-lat);  dlat=mean(differenceinlat);

    %Find the length in x and y directions that mesh was moved.
    %If this value is very close to 85000, then it is likely that error
    %that needs to be corrected
    [S] = (1/sqrt(2))*distance(d.origin(1),d.origin(2),d.origin(1)+dlat,d.origin(2)+dlong,6371000);


    if S>median(m.dx) && S>median(m.dy) && unique(d.zrot)==0 %If the distance is greater than a median grid cell, then correct it

        %Recalculate the new long and lat positions from the new site locations
        %and compare to the true lat and long positions
        [long,lat] = utm2geo(d.y-S+500000,d.x-S,d.origin(2),d.origin(1));

        differenceinlong=(d.loc(:,2)-long); dlong2=mean(differenceinlong);
        differenceinlat=(d.loc(:,1)-lat);  dlat2=mean(differenceinlat);

        [Scheck] = (1/sqrt(2))*distance(d.origin(1),d.origin(2),d.origin(1)+dlat2,d.origin(2)+dlong2,6371000);

        if Scheck<median(m.dx) && Scheck<median(m.dy) %If the new distance is less than the median grid cell

            if S>84000 && S<86000
                disp('+85000 Error is being corrected.')
            else
                disp('Station model coordinates are incorrect and have been shifted (But this is not due to the 85000 error)')
            end

            disp(['The difference of coords between the UTM2GEO result and the original was ',...
            num2str(dlong),' in longitude and ', num2str(dlat), ' in latitude.']);

            disp(['This was corrected and is now ',...
            num2str(dlong2),' in longitude and ', num2str(dlat2), ' in latitude.']);

            disp(['The mesh x and y coordinates were each moved by ',num2str(S/1000),' km before lat-long conversion']) 

            %Convert y model coordinates to longitude referenced from the center of
            %mesh origin in lat-long and subtract the 85,000 (e.g. "S") from
            %the model coordinates
            [m.lon,~]=utm2geo(m.cy+500000-S, d.origin(1)*ones(length(m.cy),1),d.origin(2), d.origin(1));

            %Convert x model coordinates to latitude with 85,000 subtracted
            [~,m.lat]=utm2geo(ones(length(m.cx),1)*500000,m.cx-S, d.origin(2), d.origin(1));

            %Build meshgrid of longitudes and latitudes.
            [m.LON,m.LAT]=meshgrid(m.lon,m.lat);            

        else %If the new distance is still greater than a grid cell then something else bad is going on
            disp('Your mesh coordinates and site coordinates do not match. Something may be wrong with your inversion files.')
            disp('If your model cells are very small (e.g. <50 m), you may need to correct the origin to more decimal places.')
        end
        
    else
        if unique(d.zrot)~=0
            disp(['The mesh is rotated ', num2str(unique(d.zrot)),' degrees']) 

            %(1) try rotating cell centers in m.X and m.Y first, then convert to LAT LON
            % also rotate station xy to lon/lat
            theta = -unique(d.zrot); % negative TO ROTATE CCW
            R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];

            YX = [m.Xc(:) m.Yc(:)]; % these are cell centers SWAPPED from load_model_modem X,Y. X IS E-W AND Y IS N-S
            rotYX=YX*R';
            Xqr = reshape(rotYX(:,1), size(m.Xc,1), []);
            Yqr = reshape(rotYX(:,2), size(m.Yc,1), []);

        %         m.X=Xqr;m.Y=Yqr; % new rotated cell centers .. switch these to keep compatibility? or just don't assign new m.X and m.Y     
        %         [m.lon,m.lat] = utm2geo(m.Y(:)+500000,m.X(:),d.origin(2),d.origin(1));
        %         m.LON = reshape(m.lon,size(m.Y));
        %         m.LAT = reshape(m.lat,size(m.X));

            [m.lon,m.lat] = utm2geo(Xqr(:)+500000,Yqr(:),d.origin(2),d.origin(1));
            m.LON = reshape(m.lon,size(Yqr));
            m.LAT = reshape(m.lat,size(Xqr));                            

            % commented below - is this necessary??? 
            
            % rotate and store lon/lat of rotated station coords - since they
            % probably won't be the same as the original lon/lat due to
            % rounding, cell-centering, etc.
             yx = [d.y d.x]; % these are cell centers
             rotyx=yx*R';
 
             roty = rotyx(:,1);
             rotx = rotyx(:,2);
 
             [rotlon,rotlat]=utm2geo(roty+500000,rotx,d.origin(2),d.origin(1));

             shiftlon = mean(rotlon-d.loc(:,2));
             shiftlat = mean(rotlat-d.loc(:,1));

             d.loc(:,1) = rotlat-shiftlat; % these won't match original d.loc
             d.loc(:,2) = rotlon-shiftlon;

             m.LON = m.LON-shiftlon;
             m.LAT = m.LAT-shiftlat;

        else

           %Convert y model coordinates to longitude referenced from the center of
            %mesh origin in lat-long
            [m.lon,~]=utm2geo(m.cy+500000, d.origin(1)*ones(length(m.cy),1),d.origin(2), d.origin(1));

            %Convert x model coordinates to latitude
            [~,m.lat]=utm2geo(ones(length(m.cx),1)*500000,m.cx, d.origin(2), d.origin(1));

            %Build meshgrid of longitudes and latitudes.
            [m.LON,m.LAT]=meshgrid(m.lon,m.lat);  
        end

    end
    
%--------------------------------------------------------------------------

else %If x,y,z station locations in model coordinates do not yet exist, then add them
    %Note: If this is done, then stations will no longer be centered in
    %cells
    disp('Warning: Your data structure did not yet have x-y-z site locations in model coordinates.')
    disp('These locations have been added, but your x-y site locations are no longer cell-centered')
    disp('and your z locations are no longer projected onto the topography')
    [d.y,d.x] = geo2utm(d.loc(:,2),d.loc(:,1),d.origin(2),d.origin(1));
    d.y = d.y-500000;
    d.z = d.loc(:,3);

end

%Find NS (idx) and EW (idy) indices where sites are located in the model
%mesh
d.idx = nearestpoint(d.x,m.cx);
d.idy = nearestpoint(d.y,m.cy);
d.idz = nearestpoint(d.z,m.cz,'next');




end % END function

