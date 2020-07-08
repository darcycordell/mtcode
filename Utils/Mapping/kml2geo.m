function kml2geo(kmlfile,geofile)
% Script to read kmlfile with placemarks, and output longitude and latitude
% into MTcode style geofile
%
% Usage: kml2geo(kmlfile,geofile)
%
% Inputs: kmlfile = name of kml file, geofile = name of output file
% can also hardcode file names at top of function when using 0 inputs
%
% only tested with kml file output from Google Earth
if nargin == 0
    kmlfile = 'west shore.kml';
    geofile = 'Kinbasket_Lake_west_coords.txt'; % name of output file
end
%%

fid = fopen(kmlfile);

nlon = 0;
nlat = 0;
filetype = 'Point';
while 1
    
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    tline = strtrim(tline);
    
    if strfind(tline,'<Style id="inline">')
        filetype = 'Line';
    end
    
    if strcmp(filetype,'Line')
        
        if strfind(tline,'<coordinates>')
            tline = fgetl(fid);
            dat = strsplit(tline,' ');
            
            for i = 1:length(dat)-1
                dat2 = strsplit(dat{i},',');
                lon(i,1) = str2double(dat2{1});
                lat(i,1) = str2double(dat2{2});
            end
            
            nlat = length(lat);
            
        end 
        
    else
    
        if strfind(tline,'<coordinates>')
            %dat = extractBetween(tline,'<coordinates>','</coordinates>');
            [~,remain] = strtok(tline,'>'); dat = strtok(remain,'<');
            dat = strsplit(dat(2:end),',');

            nlon = nlon + 1;
            nlat = nlat + 1;

            lon(nlon,1) = str2double(dat(1));
            lat(nlat,1) = str2double(dat(2));

        end
    end
    
end

fclose(fid);

fid = fopen(geofile,'w');
for i = 1:nlat
    fprintf(fid,'%12.5f %12.5f\n',lon(i),lat(i));
end
fclose(fid);