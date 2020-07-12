function srtm = download_srtm(maplims)
%-------------------------------------------------------------------
% Function to download srtm3 data.
%
% Usage: srtm = download_srtm(maplims)
%
%
% "maplims" is the area to download: [min(lat),max(lat),min(lon),max(lon)]
% "srtm" is a structure containing the latitudes, longitudes and elevations. NaN values in the SRTM data are
% interpolated.
%
% The function also saves an "elevation.mat" file which contains "srtm"
% variable.
%---------------------------------------------------------------------


bound=floor(maplims); % round since all that matters is the whole number for downloading the data
% check that boundaries are within limits of SRTM data:
if bound(1) >= 60 || bound(2) >=60
    warning('SRTM data does not go above 60°N or below 60°S, therefore no elevation data is plotted at higher (or lower) latitudes')
    if bound(1) >=60
        bound(1)=59;
    end
    if bound(2) >=60
        bound(2)=59;
    end
end

% Download the SRTM data
% Note that this will not work if you do not have an internet
% connection unless the SRTM HGT files are already on your MATLAB path
% or in the current directory
srtm=readhgt(bound,'merge','srtm3');

% Decimate the data:
n = ceil(sqrt(numel(srtm.z))/1201);
if n > 1
    srtm.lon = srtm.lon(1:n:end);
    srtm.lat = srtm.lat(1:n:end);
    srtm.z = srtm.z(1:n:end,1:n:end);
end
z=double(srtm.z);
clear srtm.z;
z(z==-32768)=NaN;
z=inpaint_nans(z);
srtm.z=z;
save('elevation.mat','srtm')

end


