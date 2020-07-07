function [d] = load_data_wsinv(name,matfile,varargin)
%Function which loads a WSINV3D data or inversion response file and puts it
%into a consistent format
%
% Usage: [d] = load_data_wsinv(name,matfile,varargin)
%
% "d" is output data structure
% "name" is a WSINV filename string which must be on the current path.
% "matfile" is the matfile output from M3 which must be on the current
% path. this is necessary to get the stations in geographic coordinates
% because the wsinv data file only contains model coordinates. this code
% supports new and legacy format matfiles.. but not fully! debugged
%
% "varargin" can be either an inversion type or the name of the inversion startup
% file
% if using "invtype", it must be a scalar:
%   1 = full impedance
%   2 = Off-diagonal impedance only
%   3 = Tipper only
%   4 = Off-diagonal impedance and tipper
%   5 = Full impedance and tipper
%
% if using startup file, it must be a string 
%
% varargin is intended for flexibility so that this code can be called
% outside of S3D
%
% If data does not exist at a particular period or component or site, it is
% set to NaN. If the input errors on all impedances is greater than 
% 10^10, then the errors are set to NaN as well (this is the case for
% inversion responses which do not have error bars).
%
%The data structure includes:
% T = period in seconds (increasing periods)
% f = frequency in Hertz (decreasing frequency)
% zrot = impedance rotation in degrees (clockwise from north)
% trot = tipper rotation in degrees (clockwise from north)
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
% the WSINV data file.
%
% d.x, d.y, d.z = MT site locations in x and y model coordinates (x is N-S)
% d.origin = The center of the model mesh in [lat, long]
%
% Note: the WSINV data file and does not include any information
% about data rotations. As such, the data rotation information is taken
% from the matfile
%
%
% Potential Bugs: ------------------------------------------------------
%
% -When should tipper be conjugated?
% -It is assumed that the WSINV data file is in e^-iwt format. Is this always true?
% -It is assuemd that the WSINV data file is in SI units. Is this always true?
%
%
%------------------------------------------------------------------------

if nargin < 3
    error('Not enough input arguments')
elseif nargin > 3
    error('Too many input arguments')
end

if isnumeric(varargin{1})
    invtype = varargin{1};
    Z_err_floor = 0;
    tip_err_floor = 0;
else
    startupfile = varargin{1};
    [invtype, ~, Z_err_floor, tip_err_floor, ~, ~] = load_startup_wsinv(startupfile);
end

disp('Loading WSINV Data File')

mu = 4*pi*10^-7;

mat = load(matfile);
if isfield(mat,'d') % if matfile is new format
    d = mat.d;
    
else % if matfile is legacy format
    [d] = make_nan_data;
    d.loc = [mat.coord(:,2) mat.coord(:,1) mat.coord(:,3)];
    d.site = mat.site';
    d.origin(1)  = (max(d.loc(:,1))+min(d.loc(:,1)))/2;
    d.origin(2) = (max(d.loc(:,2))+min(d.loc(:,2)))/2;
    d.origin(3) = 0;   
    
end


d.niter = ''; %Number of iterations left blank and can be added later
d.name = name;

d.T = []; d.f = [];
d.ns = 0; d.nf = 0; d.nr = 0;
d.responses = {};

fid = fopen(name,'r');

tmp=fscanf(fid,'%f',3); % reading first line in the file that gives number of sites, period and responses
d.ns=tmp(1); d.nf=tmp(2); d.nr=tmp(3); %separating the above informations to individual information

fscanf(fid,'%s',2); 
d.x=fscanf(fid,'%f',d.ns);  %reading the ns site locations
fscanf(fid,'%s',2); 
d.y=fscanf(fid,'%f',d.ns);  %reading the ew site locations

d.z = zeros(d.ns,1); %assume there is no topography in WSINV so site locations are at zero elevation.

all_data = nan(d.nf,d.nr,d.ns);
all_err_map = nan(d.nf,d.nr,d.ns);
all_error = nan(d.nf,d.nr,d.ns);
d.T = zeros(d.nf,1);


for i=1:d.nf % for the remaining period/frequency - 2 to nf(max. number of period)
    fscanf(fid,'%s',1); 
    d.T(i)=fscanf(fid,'%f',1);
    all_data(i,:,:)=reshape(fscanf(fid,'%f',d.nr*d.ns),d.nr,d.ns);
end % end reading the impedances for all periods

while 1

    tmp = fscanf(fid,'%s',1);
    
    if isempty(tmp) || strcmp(tmp,'#Iteration')
        break
    else
        for i=1:d.nf % for the remaining period/frequency - 2 to nf(max. number of period)
            fscanf(fid,'%f',1);
            all_error(i,:,:)=reshape(fscanf(fid,'%f',d.nr*d.ns),d.nr,d.ns);
            fscanf(fid,'%s',1); 
        end

        for i=1:d.nf % for the remaining period/frequency - 2 to nf(max. number of period)
            fscanf(fid,'%s',1); 
            all_err_map(i,:,:)=reshape(fscanf(fid,'%f',d.nr*d.ns),d.nr,d.ns);
            fscanf(fid,'%s',1);
        end  

    end % end reading the impedances for all periods
    
end

all_data(all_err_map==999)=NaN;
all_error(all_err_map==999)=NaN;
all_error(isnan(all_err_map)) = NaN;

d.f = 1./d.T;

d.Z = nan(d.nf,4,d.ns);
d.Zerr = nan(d.nf,4,d.ns);
d.tip = nan(d.nf,2,d.ns);
d.tiperr = nan(d.nf,2,d.ns);

if invtype == 1
    
    if d.nr == 8
        d.Z = all_data(:,[1 3 5 7],:)+1i*all_data(:,[2 4 6 8],:);
        d.Zerr = all_error(:,[1 3 5 7],:)+1i*all_error(:,[2 4 6 8],:);
        d.responses = {'ZXX','ZXY','ZYX','ZYY'};
    else
        error('Data File Error: Inversion type specified is not compatible with data file supplied. Too many or too few responses');
    end
    
elseif invtype == 2
        
    if d.nr == 4
        d.Z(:,[2 3],:) = all_data(:,[1 3],:)+1i*all_data(:,[2 4],:);
        d.Zerr(:,[2 3],:) = all_error(:,[1 3],:)+1i*all_error(:,[2 4],:);
        d.responses = {'ZXY','ZYX'};
    else
        error('Data File Error: Inversion type specified is not compatible with data file supplied. Too many or too few responses');
    end
    
elseif invtype == 3
    
    if d.nr == 4
        d.tip = all_data(:,[1 3],:)+1i*all_data(:,[2 4],:);
        d.tiperr = all_error(:,[1 3],:)+1i*all_error(:,[2 4],:);
        d.responses = {'TX','TY'};
    else
        error('Data File Error: Inversion type specified is not compatible with data file supplied. Too many or too few responses');
    end
    
elseif invtype == 4
    
    if d.nr == 8
        d.Z(:,[2 3],:) = all_data(:,[1 3],:)+1i*all_data(:,[2 4],:);
        d.Zerr(:,[2 3],:) = all_error(:,[1 3],:)+1i*all_error(:,[2 4],:);
        d.tip = all_data(:,[5 7],:)+1i*all_data(:,[6 8],:);
        d.tiperr = all_error(:,[5 7],:)+1i*all_error(:,[6 8],:);
        d.responses = {'ZXY','ZYX','TX','TY'};
    else
        error('Data File Error: Inversion type specified is not compatible with data file supplied. Too many or too few responses');
    end
    
elseif invtype == 5
    
    if d.nr == 12
        d.Z = all_data(:,[1 3 5 7],:)+1i*all_data(:,[2 4 6 8],:);
        d.Zerr = all_error(:,[1 3 5 7],:)+1i*all_error(:,[2 4 6 8],:);
        d.tip = all_data(:,[9 11],:)+1i*all_data(:,[10 12],:);
        d.tiperr = all_error(:,[9 11],:)+1i*all_error(:,[10 12],:);
        d.responses = {'ZXX','ZXY','ZYX','ZYY','TX','TY'};
    else
        error('Data File Error: Inversion type specified is not compatible with data file supplied. Too many or too few responses');
    end
    
else
    error('Unknown invtype')
    
end

% Apply error floors if startup file was an input
% at this point, Zerr and tiperr are from the wsinv data file. real and
% imaginary parts should have the same error

Zerrf = squeeze(0.01.*Z_err_floor.*sqrt(abs(d.Z(:,2,:).*d.Z(:,3,:)))); % computed with value from startup file

for ic = 1:4 % loop over Z components
    tmp = squeeze(real(d.Zerr(:,ic,:))) ;
    [ind] = find(tmp<Zerrf);
    tmp(ind) = Zerrf(ind);
    tmp = tmp+1i.*tmp;
    d.Zerr(:,ic,:) = reshape(tmp,d.nf,1,d.ns);
end

tiperrf = squeeze(0.01.*tip_err_floor .* sqrt(abs((d.tip(:,1,:).^2) + (d.tip(:,2,:).^2)))); % not 100% sure if this is how wsinv does it

for ic = 1:2
    tmp = squeeze(real(d.tiperr(:,ic,:))) ;
    [ind] = find(tmp<tiperrf); % linear indices work here
%     [row,col] = ind2sub(size(tmp),ind);
    tmp(ind) = tiperrf(ind);
%     tmp(row,col) = tiperrf(row,col);
    tmp = tmp+1i.*tmp;
    d.tiperr(:,ic,:) = reshape(tmp,d.nf,1,d.ns);
end
   
d.nr = d.nr/2;
d.Z = conj(d.Z); %Note: We assume that Z is always in e^{+iwt} for MT data structure. WSINV uses e^{-iwt} so we need to cnjugate

%Calculate apparent resistivity and phase from impedances
[d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);  
if isfield(mat,'d')
    rot = unique(d.zrot);
    if numel(rot)>1
        warning('More than one rotation found in d structure')
    end
    d.zrot = [];
    d.trot = [];
    for is = 1:d.ns
        d.zrot(:,is) = rot.*ones(d.nf,1);
        d.trot(:,is) = rot.*ones(d.nf,1);
    end
else
    
    for is = 1:d.ns
        d.zrot(:,is) = mat.station(is).zrot*ones(d.nf,1);
        d.trot(:,is) = mat.station(is).trot*ones(d.nf,1);
    end
end

d = set_map_projection(d);

end
