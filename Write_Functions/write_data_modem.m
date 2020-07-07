function write_data_modem(datafile,d)
% Function to write a data file for use with the modEM 3D inversion
%
% Usage: write_data_modem(datafile,d)
%
% Inputs: datafile is a string file name
%       d is a data structure
%
% Impedance units are assumed to be [V/m]/[T]
% sign convention is assumed to be -1, i.e.: exp(-i\omega t)
% orientation is assumed to be 0.0 (unrotated)
% origin is assumed to be 0,0,0
%
% siteName is a cell array of characters containing site names
% coord is a [ns x 3] matrix containing the site locations in x, y, and z
%       (in meters)
% data is a complex valued matrix of size [nr x nf x ns], where nr is
% either 4,8, or 12
% (diagonals, full tensor, full tensor and tipper)
%
% Mostly copied from writeZ_3D included with the modEM software
%
% BL 2019
% Now calculates stations and frequencies missing impedance and/or tipper
% data so that the correct numbers are written in headers.
%
% DC Edits, June 2, 2015:
%       Bug Fix: When you hit "Tipper Only", the data file was writing out
%       impedance only. I added some if statements to ensure that all the
%       different data types are written out correctly

% DC Edits, August 25, 2015
%       Added an all() logical check when writing out impedances and tipper
%       so that if there is no data, these data points are not written out.

% Greg Nieuwenhuis, July 2013
%-------------------------------------------------------------
disp('Writing ModEM Data File')

mu = 4*pi*10^-7;
d.Z = d.Z/mu/1000;
d.Zerr = d.Zerr/mu/1000;

d.Zerr(isnan(d.Zerr)) = 2*10^15;
d.tiperr(isnan(d.tiperr)) = 2*10^15;

if all([d.zrot(:); d.trot(:)]==d.zrot(1))
    d.rotation = d.zrot(1); 
else
    d.rotation = 0;
    disp(['WARNING: The rotations in your impedances and/or tipper data are different at different periods or sites. ', ...
        'Rotation set to zero.'])
end

% check stations and frequencies missing data for headers
noZ_s = numel(find(squeeze(sum(sum(isnan(d.Z),1),2))==d.nf*4)); % find number of stations with no impedance data
noZ_f = numel(find(squeeze(sum(sum(isnan(d.Z),3),2))==d.ns*4)); % find number of frequencies with no impedance data

notip_s = numel(find(squeeze(sum(sum(isnan(d.tip),1),2))==d.nf*2)); % find number of stations with no tipper data
notip_f = numel(find(squeeze(sum(sum(isnan(d.tip),3),2))==d.ns*2)); % find number of frequencies with no tipper data

%  description <= 80 char in length
header = 'written by Matlab function write_data_modem';

% Assume +iwt convention for all data
signstr = 'exp(+i\omega t)';

% Open the file
fid = fopen(datafile,'w');

r = char(d.responses)';
r = r(:)';

%Check if impedances exist
if ~isempty(strfind(r,'ZXY'))
    
    if ~isempty(strfind(r,'ZXX')) %If diagonal components exist then we have full impedance
        dataType = 'Full_Impedance';
    else %else we have only off-diagonal impedances
        dataType = 'Off_Diagonal_Impedance';
    end

    %  assume SI units by default as in ModEM; alternative is [mV/km]/[nT] 
    units = '[mV/km]/[nT]';
    
    %  file header: the info line
    fprintf(fid,'# %s\n',header);

    %  data type header: comment, then six lines
    comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
    fprintf(fid,'# %s\n',comment);
    fprintf(fid,'> %s\n',dataType);
    fprintf(fid,'> %s\n',signstr);
    fprintf(fid,'> %s\n',units);
    fprintf(fid,'> %.2f\n',d.rotation);
    fprintf(fid,'> %.4f %.4f\n',d.origin(1:2));
    fprintf(fid,'> %d %d\n',d.nf-noZ_f,d.ns-noZ_s); % subtract stations/frequencies with no impedance data

    %  now write all impedances line by line

    for k = 1:d.ns
        for i = 1:d.nf
            for j = 1:4
                if ~isnan(real(d.Z(i,j,k))) && real(d.Z(i,j,k))~= 0 
                    fprintf(fid,'%12.6E ',d.T(i)); % transmitter
                    fprintf(fid,'%s %9.4f %9.4f ',d.site{k},d.loc(k,1),d.loc(k,2)); % receiver
                    fprintf(fid,'%12.3f %12.3f %12.3f ',[d.x(k) d.y(k) d.z(k)]); % receiver x,y,z
                    fprintf(fid,'%s %15.6E %15.6E %15.6E\n',d.responses{j},real(d.Z(i,j,k)),imag(d.Z(i,j,k)),real(d.Zerr(i,j,k))); % data
                end
            end
        end
    end

end


%Check if tipper data exists
if ~isempty(strfind(r,'TX'))

    dataType = 'Full_Vertical_Components';
    units = '[]';

    %  file header: the info line
    fprintf(fid,'# %s\n',header);

    %  data type header: comment, then six lines
    comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
    fprintf(fid,'# %s\n',comment);
    fprintf(fid,'> %s\n',dataType);
    fprintf(fid,'> %s\n',signstr);
    fprintf(fid,'> %s\n',units);
    fprintf(fid,'> %.2f\n',d.rotation);
    fprintf(fid,'> %.4f %.4f\n',d.origin(1:2));        
    fprintf(fid,'> %d %d\n',d.nf-notip_f,d.ns-notip_s); % subtract stations/frequencies with no tipper data
    
    %  now write all tippers line by line
    for k = 1:d.ns
        for i = 1:d.nf
            for j = 1:2
                if ~isnan(real(d.tip(i,j,k))) && real(d.tip(i,j,k))~= 0 
                    fprintf(fid,'%12.6E ',d.T(i)); % transmitter
                    fprintf(fid,'%s %9.4f %9.4f ',d.site{k},d.loc(k,1),d.loc(k,2)); % receiver
                    fprintf(fid,'%12.3f %12.3f %12.3f ',[d.x(k) d.y(k) d.z(k)]); % receiver x,y,z
                    fprintf(fid,'%s %15.6E %15.6E %15.6E\n',d.responses{end-2+j},real(d.tip(i,j,k)),imag(d.tip(i,j,k)),real(d.tiperr(i,j,k))); % data
                end
            end
        end
    end

end

fclose(fid);
