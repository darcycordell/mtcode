function [new_elev] = get_srtm(maplim)
    bound=floor(maplim); % round since all that matters is the whole number for downloading the data    
    if bound(3) >= 60 || bound(4) >=60 % check that boundaries are within limits of SRTM data:
        warning('SRTM data does not go above 60°N or below 60°S, therefore no elevation data is plotted at higher (or lower) latitudes')
        if bound(3) >=60
            bound(3)=59;
        end
        if bound(4) >=60
            bound(4)=59;
        end
    end
    
    new_elev=readhgt([bound(3):bound(4)],[bound(1):bound(2)],'merge','srtm3');
    % decimate the data:
    n = ceil(sqrt(numel(new_elev.z))/2401); % try 2401 to keep larger grids
    if n > 1
        new_elev.lon = new_elev.lon(1:n:end);
        new_elev.lat = new_elev.lat(1:n:end);
        new_elev.z = new_elev.z(1:n:end,1:n:end);
    end
    z=double(new_elev.z);
    clear elev.z;
    ind=find(z==-32768);
    z(ind)=NaN;
    new_elev.z=z;
        
end