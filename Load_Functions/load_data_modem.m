function [d] = load_data_modem(name)
%Function which loads a ModEM data or inversion response file and puts it
%into a consistent format.
%
% Usage: [d] = load_data_modem(name)
%
% "d" is output data structure
% "name" is a ModEM filename string which must be on the current path.
%
% If data does not exist at a particular period or component or site, it is
% set to NaN. If the input errors on all impedances is greater than 
% 10^10, then the errors are set to NaN as well (this is the case for
% inversion responses which do not have error bars).
%
%The data structure includes:
% T = period in seconds (increasing periods)
% f = frequency in Hertz (decreasing frequency)
% zrot = impedance data rotation in degrees (clockwise from north) (nf x ns)
% trot = tipper data rotation in degrees (clockwise from north) (nf x ns)
% site = MT site names
% loc = MT site locations in latitude, longitude and elevation
%           Note: elevations in ModEM data file must be in meters below sea
%           level
% ns, nf, and nr are the number of stations, frequencies and responses
% components = the string of components which are included.
% Z = impedance in SI units
% Zerr = impedance error in SI units
% rho = apparent resistivity in Ohm m
% rhoerr = apparent resistivity error in Ohm m
% pha = phase in degrees
% phaerr = phase error in degrees
% tip = tipper
% tiperr = tipper error
%
% This function also outputs the data locations in the x-y model space and
% the latitude and longitude origin of the stations which is included in
% the ModEM data file.
%
% d.x, d.y, d.z = MT site locations in x and y model coordinates (x is N-S)
% d.origin = The center of the model mesh in [lat, long]
%
% For more info on ModEM see:
%   Kelbert, A., Meqbel, N., Egbert, G., & Tandon, K. (2014). ModEM: a modular 
%   system for inversion of electromagnetic geophysical data. Computers & Geosciences, 66, 40–53. 
%
% To be added:----------------------------------------------------------
%
%
% Potential Bugs: ------------------------------------------------------
%
% -When should tipper be conjugated?
%
%
%------------------------------------------------------------------------

disp('Loading ModEM Data File')

mu = 4*pi*10^-7;

%Initial variables
[d] = make_nan_data;
d.niter = ''; %Number of iterations left blank and can be added later
d.name = name;
d.T = []; d.f = [];
d.zrot = []; d.trot = [];
d.site = {}; d.loc = [];
d.ns = 0; d.nf = 0; d.nr = 0;
d.responses = {};

%Read data file
fid = fopen(name,'r');
nlines = 0; hc = 0;
while 1
   
    line = fgetl(fid);
    
    if isempty(line)
        line = ' ';
    end
    
    if line==-1
        break
    end
    
    if ~isempty(str2num(line(1)))
        nlines = nlines+1;
        data = textscan(line,'%f %s %f %f %f %f %f %s %f %f %f');
        
        T(nlines) = [data{1}];
        site(nlines) = [data{2}];
        loc(nlines,1) = [data{3}];
        loc(nlines,2) = [data{4}];
        loc(nlines,3) = [data{7}];
        
        x(nlines) = [data{5}];
        y(nlines) = [data{6}];
        z(nlines) = [data{7}];
        
        responses(nlines) = [data{8}];
        D(nlines) = [data{9}]+1i*[data{10}];
        Derr(nlines) = [data{11}]+1i*[data{11}];
        
       
    else
        % read the header info
        tmp = textscan(fid,'> %s',6,'delimiter','\n');
        tmp = char(tmp{1});
          
        if ~isempty(tmp)
            hc = hc+1;
            header{hc} = tmp;
            % read the block header
            dataType{hc} = strtrim(header{hc}(1,:));
            signstr{hc}  = header{hc}(2,:);
            if strfind(signstr{hc},'-')
                isign(hc) = -1;
            else
                isign(hc) = +1;
            end
            units{hc} = strtrim(header{hc}(3,:));
            rotation(hc) = sscanf(header{hc}(4,:),'%f',1);
            origin(hc,:) = sscanf(header{hc}(5,:),'%f',2);
            nfns{hc}  = sscanf(header{hc}(6,:),'%d',2);
        end
        
    end
   
    
end


fclose(fid); %Done reading text file

%Site names are in data column 2
[d.site,inds] = unique(site);
d.site = d.site';
d.ns = length(d.site);

%Site locations in lat-long are in columns 4 (longitude) and 3 (latitude)
% Site elevations are in column 7
d.loc = loc(inds,:);

%Site locations in x-y model space are in columns 6 (NS) and 5 (EW)
d.x = x(inds)'; 
d.y = y(inds)'; 
d.z = z(inds)';

%Periods are in column 1
d.T = unique(T');
d.f = 1./d.T;
d.nf = length(d.f);

%Find all the responses that are included in this dataset
d.responses = unique(responses,'stable');
d.nr = length(d.responses);

components = {'ZXX','ZXY','ZYX','ZYY','TX','TY'};

%Initialize all the data variables in the "d" structure
d.Z = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns); d.Zerr = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns);
d.rho = nan(d.nf,4,d.ns); d.rhoerr = nan(d.nf,4,d.ns);
d.pha = nan(d.nf,4,d.ns); d.phaerr = nan(d.nf,4,d.ns);
d.tip = nan(d.nf,2,d.ns)+1i*nan(d.nf,2,d.ns); d.tiperr = nan(d.nf,2,d.ns)+1i*nan(d.nf,2,d.ns);

for i = 1:nlines
    
    [~, ifreq] = max(T(i)==d.T); %frequency index
    [~, isite] = max(strcmp(site(i),d.site)); %site index
    [~, icomp] = max(strcmp(responses(i),components)); %component index (xx,yy, etc.)
    
    if icomp == 5 || icomp == 6
        %If the component is 5 or 6, this signifies a tipper component
        d.tip(ifreq,icomp-4,isite) = D(i); %Tx component
        d.tiperr(ifreq,icomp-4,isite) = Derr(i); %Ty component
        
    else
        %Else the component is an impedance component
        d.Z(ifreq,icomp,isite) = D(i);
        d.Zerr(ifreq,icomp,isite) = Derr(i);
    end
    
end

%Convert to SI units if necessary
%Bug fix to index units(i) instead of units(hc).
for i = 1:hc
    if strcmp(units(i),'[mV/km]/[nT]') == 1
        d.Z = d.Z*mu*1000;
        d.Zerr = d.Zerr*mu*1000;
    elseif strcmp(units(i),'[V/m]/[T]') == 1 || strcmp(units(i),'[]') == 1
        %Leave in SI units if SI units are specified or if no units are
        %specified
    else
        error('There is something wrong with the units in your datafile')
    end
end

if ~all(isign==isign(1)) || ~all(rotation==rotation(1)) || ~all(origin(:,1)==origin(1,1)) || ~all(origin(:,2)==(origin(1,2)))
    error('Your sign convention, rotation or mesh origin are inconsistent between impedance/tipper data types. load_data_modem is not compatible')
else
    if rotation(1)~=0
        %ModEM data format files have only one rotation applied to all the data
        %In general this number is zero. Caution is required if rotation ~=0
        disp('Caution: The rotation in your ModEM data file is not zero. This is non-standard and may cause problems when using these codes.')
    end
    
    d.zrot = rotation(1)*ones(d.nf,d.ns);
    d.trot = d.zrot;
    d.origin = origin(1,:)';
    d.origin(3) = 0;
    
end


if isign(1) == -1
    %If the sign convention is e^{-iwt} then you must conjugate the
    %impedances to correct.
    d.Z = conj(d.Z);
    
    %Potential Bug: should tipper be conjugated as well?
    %d.tip = conj(d.tip);
    
elseif isign(1) == 1
    %If the sign convention is e^{+iwt} then do nothing
    
else
    error('Your sign convention must be either e^{+iwt} or e^{-iwt}')
end

%If all the errors are >10^10, then set them as NaN. This is for the case
%when you are loading an inversion response or synthetic data with no
%errors
if d.Zerr(isnan(d.Zerr)==0) > 10^10
    d.Zerr = NaN*d.Zerr;
end

if d.tiperr(isnan(d.tiperr)==0) > 10^10
    d.tiperr = NaN*d.tiperr;
end

%Calculate apparent resistivity and phase from impedances
[d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);

%Set up the map limits
d = set_map_projection(d);


end

