function [d,units_code] = load_data_edi(filename)
%
% Function which takes a given EDI filename and loads it into the standard
% data structure format (SI unit impedance with e^{+iwt} convention).
%
% Usage: [d,units_code] = load_edi(filename)
%
% INPUTS
%
% filename is a string that must be on the current path
%
% OUTPUTS
%
% d is data structure (See below)
% unit_code = 0 (too small rho flag <1 ohm m), 1 (goldilocks rho), 2 (too high rho flag >10^6 ohm m)
%
%
% Note: This function *assumes* the EDI file is in field units (mV/km /nT) with
% E = Z*B. This function automatically converts the impedances to SI units (Volts / Amps)
% with E = Z*H.
%
% If EDI file does not contain impedances, they are calculated from the apparent
% resistivity and phase. The final outputted apparent resistivity and phase
% is ALWAYS calculated from the SI unit impedances.
%
% If you happen to have a rare EDI which is already in SI units, then you will have
% to use convert_edi_units.m to write out a new EDI in field units.
%
% OUTPUTS
%
% d is the data structure
%
% The data structure includes:
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
% Zerr = impedance standard error in SI units (Note this is converted from EDI variances with sigma = sqrt(variance))
% rho = apparent resistivity in Ohm m
% rhoerr = apparent resistivity error in Ohm m
% pha = phase in degrees
% phaerr = phase error in degrees
% tip = tipper
% tiperr = tipper error
%
%
%This is largely copy pasted from read_edi_imp.m (written by Dennis Rippe and
%Greg Nieuwenhuis)
%
% POTENTIAL BUGS:
%
% If/when to conjugate tipper components?
%
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)
%

% Reading in edi file
fid=fopen(filename,'rt'); linecount = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if strncmp(strtrim(tline),'EMPTY=',6)
        [~,empty_str]=strtok(tline,'=');
        ind=strfind(empty_str,'=');
        empty_str(ind)=[];
        empty=str2double(empty_str);
    end
    if strncmp(strtrim(tline),'LAT=',4)
        tline=strtrim(tline);
        lat=tline(5:length(tline));
    end
    if strncmp(strtrim(tline),'LONG=',5)
        tline=strtrim(tline);
        long=tline(6:length(tline));
    end
    if strncmp(strtrim(tline),'ELEV=',5)
        [~,tline]=strtok(tline,'=');
        ind=strfind(tline,'=');
        tline(ind)=[];
        elev=str2double(tline);
    end
    if strncmp(strtrim(tline),'REFLAT=',7)
        tline=strtrim(tline);
        reflat=tline(8:length(tline));
    end
    if strncmp(strtrim(tline),'REFLONG=',8)
        tline=strtrim(tline);
        reflong=tline(9:length(tline));
    end
    if strncmp(strtrim(tline),'REFELEV=',8)
        [~,tline]=strtok(tline,'=');
        ind=strfind(tline,'=');
        tline(ind)=[];
        refelev=str2double(tline);
    end
    if strncmp(strtrim(tline),'>FREQ',5)
        nfreq=find_nfreq(tline);
        f(:,1)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZROT',5)
        nfreq=find_nfreq(tline);
        zrot=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXXR',5)
        nfreq=find_nfreq(tline);
        if ~exist('zrot','var')
            zrot = str2double(tline(11));
        end
        zxxr=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXXI',5)
        nfreq=find_nfreq(tline);
        zxxi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXX.VAR',8)
        nfreq=find_nfreq(tline);
        Zvar(1,1,1:nfreq)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXYR',5)
        nfreq=find_nfreq(tline);
        ind1=strfind(tline,'ROT='); %Now make sure the impedance rotation is defined properly
        if ~isempty(ind1)
            ind1=ind1+4;
            ind2=strfind(tline,'//')-2;
            if strcmp(tline(ind1:ind2),'NORTH')
                zrot=zeros(nfreq); %if the rotation is north, then zrot=0
            elseif strcmp(strtrim(tline(ind1:ind2)),'ZROT')
                %don't do anything because ZROT should already be defined
                if ~exist('zrot','var')
                    error('zrot dosn''t appear to be defined, but is used!')
                end
            elseif strcmp(tline(ind1:ind2),'NONE')
                zrot=zeros(1,nfreq);
            else
                %if ~exist('zrot') == 1
                    error('The impedance rotation is undefined')
                %end
            end
        else
            warning('Not sure if the impedance rotation is defined properly!')
        end
        zxyr=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXYI',5)
        nfreq=find_nfreq(tline);
        zxyi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZXY.VAR',8)
        nfreq=find_nfreq(tline);
        Zvar(1,2,1:nfreq)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYXR',5)
        nfreq=find_nfreq(tline);
        zyxr=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYXI',5)
        nfreq=find_nfreq(tline);
        zyxi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYX.VAR',8)
        nfreq=find_nfreq(tline);
        Zvar(2,1,1:nfreq)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYYR',5)
        nfreq=find_nfreq(tline);
        zyyr=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYYI',5)
        nfreq=find_nfreq(tline);
        zyyi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>ZYY.VAR',8)
        nfreq=find_nfreq(tline);
        Zvar(2,2,1:nfreq)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TROT.EXP',9) || strncmp(strtrim(tline),'>TROT',5)
        nfreq=find_nfreq(tline);
        trot=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TXR.EXP',8) || strncmp(strtrim(tline),'>TXR',4) && ~strncmp(strtrim(tline),'>TXRVAR.EXP',11)
        nfreq=find_nfreq(tline);
        ind1=strfind(tline,'ROT='); %Now make sure the impedance rotation is defined properly
        if ~isempty(ind1)
            ind1=ind1+4;
            ind2=strfind(tline,'//')-2;
            if strcmp(tline(ind1:ind2),'NORTH')
                if exist('trot','var')
                    if all(trot ~= 0)
                        error('trot is defined as something other than zero, but ROT=NORTH.....nonsense!')
                    end
                else
                    trot=zeros(nfreq); %if the rotation is north, then zrot=0
                end
            elseif strcmp(strtrim(tline(ind1:ind2)),'TROT')
                %don't do anything because TROT should already be defined
                if ~exist('trot','var')
                    error('trot dosn''t appear to be defined, but is used!')
                end
            elseif strcmp(tline(ind1:ind2),'NONE')
                trot = zeros(1,nfreq);
            else
                error('Tipper rotation is undefined (ROT= does not match expected values)')
            end
        else %if it doesn't say what the rotation is at all
            if ~exist('trot','var') %check if trot exists, if not, there is a problem
                warning('Tipper rotation is undefined, assumed to be the same as impedance')
                trot=zrot;
            end
        end
        txr=fscanf(fid,'%f',nfreq);% read in the data
    end
    if strncmp(strtrim(tline),'>TXI.EXP',8) || strncmp(strtrim(tline),'>TXI',4)
        nfreq=find_nfreq(tline);
        txi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TXVAR.EXP',10) || strncmp(strtrim(tline),'>TX.VAR',7) || strncmp(strtrim(tline),'>TXRVAR.EXP',11)
        nfreq=find_nfreq(tline);
        Tvar(1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TYR.EXP',8) || strncmp(strtrim(tline),'>TYR',4) && ~strncmp(strtrim(tline),'>TYRVAR.EXP',11)
        nfreq=find_nfreq(tline);
        tyr=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TYI.EXP',8) || strncmp(strtrim(tline),'>TYI',4)
        nfreq=find_nfreq(tline);
        tyi=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>TYVAR.EXP',10) || strncmp(strtrim(tline),'>TY.VAR',7) || strncmp(strtrim(tline),'>TYRVAR.EXP',11)
        nfreq=find_nfreq(tline);
        Tvar(2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOROT',7)
        nfreq=find_nfreq(tline);
        rhorot=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOXX',6) && ~strncmp(strtrim(tline),'>RHOXX.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoa(1,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOXX.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoaerr(1,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOXY',6)&& ~strncmp(strtrim(tline),'>RHOXY.ERR',10)
        nfreq=find_nfreq(tline);
        ind1=strfind(tline,'ROT='); %Now make sure the impedance rotation is defined properly
        if ~isempty(ind1)
            ind1=ind1+4;
            ind2=strfind(tline,'//')-2;
            if strcmp(tline(ind1:ind2),'NORTH')
                if exist('rhorot','var')
                    if all(rhorot ~= 0)
                        error('rhorot is defined as something other than zero, but ROT=NORTH.....nonsense!')
                    end
                else
                    rhorot=zeros(nfreq); %if the rotation is north, then rhorot=0
                end
            elseif strcmp(tline(ind1:ind2),'RHOROT')
                %don't do anything because ZROT should already be defined
                if ~exist('rhorot','var')
                    error('rhorot dosn''t appear to be defined, but is used!')
                end
            else
                error('The resistivity/phase rotation is undefined')
            end
        else%if it doesn't say what the rotation is at all
            if ~exist('rhorot','var') %check if rhorot exists, if not, there is a problem
                error('Resistivity/Phase rotation is undefined')
            end
        end
        Rhoa(1,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOXY.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoaerr(1,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOYX',6) && ~strncmp(strtrim(tline),'>RHOYX.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoa(2,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOYX.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoaerr(2,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOYY',6) && ~strncmp(strtrim(tline),'>RHOYY.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoa(2,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>RHOYY.ERR',10)
        nfreq=find_nfreq(tline);
        Rhoaerr(2,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSXX',6) && ~strncmp(strtrim(tline),'>PHSXX.ERR',10)
        nfreq=find_nfreq(tline);
        Phs(1,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSXX.ERR',10)
        nfreq=find_nfreq(tline);
        Phserr(1,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSXY',6) && ~strncmp(strtrim(tline),'>PHSXY.ERR',10)
        nfreq=find_nfreq(tline);
        Phs(1,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSXY.ERR',10)
        nfreq=find_nfreq(tline);
        Phserr(1,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSYX',6) && ~strncmp(strtrim(tline),'>PHSYX.ERR',10)
        nfreq=find_nfreq(tline);
        Phs(2,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSYX.ERR',10)
        nfreq=find_nfreq(tline);
        Phserr(2,1,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSYY',6) && ~strncmp(strtrim(tline),'>PHSYY.ERR',10)
        nfreq=find_nfreq(tline);
        Phs(2,2,:)=fscanf(fid,'%f',nfreq);
    end
    if strncmp(strtrim(tline),'>PHSYY.ERR',10)
        nfreq=find_nfreq(tline);
        Phserr(2,2,:)=fscanf(fid,'%f',nfreq);
    end
    
    linecount = linecount+1;
end
fclose(fid);
% Creating data array
% Coordinates
if exist('lat','var')
    if lat(1)=='-'; latsign=-1; else latsign=1; end;
    lat=regexprep(lat,':',' ');
    lat=sscanf(lat,'%f');
end
if exist('long','var')
    if long(1)=='-'; longsign=-1; else longsign=1; end
    long=regexprep(long,':',' ');
    long=sscanf(long,'%f');
end
if ~exist('elev','var'); elev=0; end
if exist('reflat','var')
    if reflat(1)=='-'; reflatsign=-1; else reflatsign=1; end
    reflat=regexprep(reflat,':',' ');
    reflat=sscanf(reflat,'%f');
end
if exist('reflong','var')
    if reflong(1)=='-'; reflongsign=-1; else reflongsign=1; end
    reflong=regexprep(reflong,':',' ');
    reflong=sscanf(reflong,'%f');
end
% Convert coordinates to decimal degrees
if exist('lat','var') && numel(lat)==3; lat=latsign*(abs(lat(1))+lat(2)/60+lat(3)/3600); end
if exist('lat','var') && numel(lat)==2; lat=latsign*(abs(lat(1))+lat(2)/60); end
if exist('long','var') && numel(long)==3; long=longsign*(abs(long(1))+long(2)/60+long(3)/3600); end
if exist('long','var') && numel(long)==2; long=longsign*(abs(long(1))+long(2)/60); end
if exist('reflat','var') && numel(reflat)==3; reflat=reflatsign*(abs(reflat(1))+reflat(2)/60+reflat(3)/3600); end
if exist('reflat','var') && numel(reflat)==2; reflat=reflatsign*(abs(reflat(1))+reflat(2)/60); end
if exist('reflong','var') && numel(reflong)==3; reflong=reflongsign*(abs(reflong(1))+reflong(2)/60+reflong(3)/3600); end
if exist('reflong','var') && numel(reflong)==2; reflong=reflongsign*(abs(reflong(1))+reflong(2)/60); end
% Coordinates
coords=zeros(3,1);
if exist('reflat','var'); coords(1)=reflat; end
if exist('reflong','var'); coords(2)=reflong; end
if exist('refelev','var'); coords(3)=refelev; end
% lat, long, elev replace reflat, reflong, refelev when present
if exist('lat','var'); coords(1)=lat; end
if exist('long','var'); coords(2)=long; end
if exist('elev','var'); coords(3)=elev; end

if exist('refelev','var') && elev == 0
    coords(3) = refelev;
end


% Rotation angles
if ~exist('zrot','var'); zrot=NaN(nfreq,1); end
if ~exist('rhorot','var'); rhorot=NaN(nfreq,1); end
if ~exist('trot','var'); trot=NaN(nfreq,1); end


% Impedances
Z=NaN(2,2,length(f))+1i.*NaN(2,2,length(f));
if exist('zxxr','var'); Z(1,1,1:length(zxxr))=zxxr+(1i.*zxxi); end
if exist('zxyr','var'); Z(1,2,1:length(zxyr))=zxyr+(1i.*zxyi); end
if exist('zyxr','var'); Z(2,1,1:length(zyxr))=zyxr+(1i.*zyxi); end
if exist('zyyr','var'); Z(2,2,1:length(zyyr))=zyyr+(1i.*zyyi); end
if ~exist('Zvar','var'); Zvar=zeros(2,2,length(f))+1i.*zeros(2,2,length(f)); end
% Tipper
T=NaN(2,length(f));
if exist('txr','var'); T(1,:)=txr(1:length(f))+(1i.*txi(1:length(f))); end % do not record any tipper values which have no frequency associated with them ... (ABT dataset!)
if exist('tyr','var'); T(2,:)=tyr(1:length(f))+(1i.*tyi(1:length(f))); end
if ~exist('Tvar','var'); Tvar=NaN(2,length(f)); end
if ~exist('T','var'); T=NaN(2,length(f))+1i.*NaN(2,length(f)); end
% Apparent resistivities
if ~exist('Rhoa','var'); Rhoa=NaN(2,2,length(f)); end
if ~exist('Rhoaerr','var'); Rhoaerr=NaN(2,2,length(f)); end
% Phases
if ~exist('Phs','var'); Phs=NaN(2,2,length(f)); end
if ~exist('Phserr','var'); Phserr=NaN(2,2,length(f)); end
% remove empty data (make it NaN)
if exist('empty','var')
    if exist('Z','var')
        indxx=find(real(Z(1,1,:))==empty);% we assume that if the real part is non-existent, so is the imaginary part
        indxy=find(real(Z(1,2,:))==empty);
        indyx=find(real(Z(2,1,:))==empty);
        indyy=find(real(Z(2,2,:))==empty);
        Z(1,1,indxx)=NaN+1i*NaN;
        Zvar(1,1,indxx)=NaN+1i*NaN;
        Z(1,2,indxy)=NaN+1i*NaN;
        Zvar(1,2,indxy)=NaN+1i*NaN;
        Z(2,1,indyx)=NaN+1i*NaN;
        Zvar(2,1,indyx)=NaN+1i*NaN;
        Z(2,2,indyy)=NaN+1i*NaN;
        Zvar(2,2,indyy)=NaN+1i*NaN;
        if all(all(all(real(Z(:,:,:))==0))) %Now look to see if all the values are zero, in which case we assume there is no data
            Z(:,:,:)=NaN(2,2,nfreq)+1i*NaN(2,2,nfreq);
            Zvar(:,:,:)=NaN(2,2,nfreq)+1i*NaN(2,2,nfreq);
        end
    end
    if exist('Rhoa','var')
        indxx=find(Rhoa(1,1,:)==empty);
        indxy=find(Rhoa(1,2,:)==empty);
        indyx=find(Rhoa(2,1,:)==empty);
        indyy=find(Rhoa(2,2,:)==empty);
        Rhoa(1,1,indxx)=NaN;
        Rhoa(1,2,indxy)=NaN;
        Rhoa(2,1,indyx)=NaN;
        Rhoa(2,2,indyy)=NaN;
        Rhoaerr(1,1,indxx)=NaN;
        Rhoaerr(1,2,indxy)=NaN;
        Rhoaerr(2,1,indyx)=NaN;
        Rhoaerr(2,2,indyy)=NaN;
        if all(all(all(Rhoa(:,:,:)==0)))
            Rhoa(:,:,:)=NaN(2,2,nfreq);
            Rhoaerr(:,:,:)=NaN(2,2,nfreq);
        end
    end
    if exist('Phs','var')
        indxx=find(Phs(1,1,:)==empty);
        indxy=find(Phs(1,2,:)==empty);
        indyx=find(Phs(2,1,:)==empty);
        indyy=find(Phs(2,2,:)==empty);
        Phs(1,1,indxx)=NaN;
        Phs(1,2,indxy)=NaN;
        Phs(2,1,indyx)=NaN;
        Phs(2,2,indyy)=NaN;
        Phserr(1,1,indxx)=NaN;
        Phserr(1,2,indxy)=NaN;
        Phserr(2,1,indyx)=NaN;
        Phserr(2,2,indyy)=NaN;
        if all(all(all(Phs(:,:,:)==0)))
            Phs(:,:,:)=NaN(2,2,nfreq);
            Phserr(:,:,:)=NaN(2,2,nfreq);
        end
    end
    if exist('T','var')
        indx=find(real(T(1,:))==empty);
        indy=find(real(T(2,:))==empty);
        T(1,indx)=NaN+1i*NaN;
        T(2,indy)=NaN+1i*NaN;
        Tvar(1,indx)=NaN+1i*NaN;
        Tvar(2,indy)=NaN+1i*NaN;
        if all(all(T(:,:)==0))
            T(:,:)=NaN(2,nfreq);
            Tvar(:,:)=NaN(2,nfreq);
        end
    end
end

[d] = make_nan_data;

d.responses = {'ZXX','ZXY','ZYX','ZYY','TX','TY'}; 

%If there is no impedance data, then set rotation to the same as rhorot
if all(all(all(isnan(Z))))
    zrot=rhorot;
end

%If there is no resistivity/phase data then set rotation to the same as
%zrot
if all(all(all(isnan(Rhoa)))) && all(all(all(isnan(Phs)))) 
    rhorot=zrot;
end

%If there is no tipper, set tipper rotation to zero and remove from
%components
if all(all(all(isnan(T))))
    trot=zeros(size(T));
    d.responses(end-1:end) = [];
end

%If the EDI contains only tipper, remove from the components:
if all(all(all(isnan(Rhoa)))) && all(all(all(isnan(Phs)))) && all(all(all(isnan(Z))))
    d.responses(1:4) = [];
end

trot(trot==10^32) = 0;

d.site = {filename(1:end-4)};

%Do some checks to make sure rotations are all ok--------------------------
if ~exist('trot','var')
  trot = zrot;   % Dangerous option that can be used to fix incorrect EDI files 
end

if ~all(zrot == zrot(1))
    disp(['WARNING: ',d.site{1},' has different impedance rotations at different periods'])
end

if ~all(rhorot == rhorot(1))
    disp(['WARNING: ',d.site{1},' has different apparent resistivity rotations at different periods'])
end 

if ~all(trot == trot(1))
    disp(['WARNING: ',d.site{1},' has different tipper rotations at different periods'])
end 

if ~all(rhorot == zrot)
    disp(['WARNING: ',d.site{1},' has different apparent resistivity rotations from the impedance rotations! Impedance rotations are taken as default'])
end 


d.zrot = zrot;
d.trot = trot;    

%Done checks---------------------------------------------------------------


%Build data structure, d, in the standard format
d.f = f;
d.T = 1./d.f;
d.nf = length(d.T);
d.ns = 1;
d.nr = 4;

mu = 4*pi*10^-7;

%Convert EDI variances to standard errors (2020/07)
%   EDI files contain variances but MTcode (and ModEM) always use standard errors as of 2020/07
Zvar = sqrt(Zvar);
Tvar = sqrt(Tvar);

d.Z = reshape(permute((Z),[2 1 3]),d.nr,d.nf).';
d.Zerr = (reshape(permute(Zvar,[2 1 3]),d.nr,d.nf)'+1i*reshape(permute(Zvar,[2 1 3]),d.nr,d.nf)');

%Potential bug: not sure if/when tipper needs to be conjugate
d.tip = T.';
d.tiperr = real(Tvar.')+1i*real(Tvar.');

d.rho = reshape(permute(Rhoa,[2 1 3]),d.nr,d.nf)';
d.rhoerr = reshape(permute(Rhoaerr,[2 1 3]),d.nr,d.nf)';
d.pha = reshape(permute(Phs,[2 1 3]),d.nr,d.nf)';
d.phaerr = reshape(permute(Phserr,[2 1 3]),d.nr,d.nf)';

d.loc = coords';

d = set_map_projection(d);

d.x = 0; d.y = 0; d.z = d.loc(1,3);
d.origin = d.loc;

d.name = filename;
d.niter = filename(1:end-4);

% Calculate missing impedances, apparent resistivities and phases
if all(all(all(isnan(d.Z))))
    [d.Z,d.Zerr]=calc_Z(d.rho,d.rhoerr,d.pha,d.phaerr,d.T);
end

if all(all(all(isnan(d.rho)))) && all(all(all(isnan(d.pha)))) 
    [d.rho,d.pha,d.rhoerr,d.phaerr]=calc_rho_pha(d.Z,d.Zerr,d.T);
end

%Convert from E/B field units (mv/km / nT) to E/H SI units (Volts / Amp)
d.Z = d.Z*(mu*1000);
d.Zerr = d.Zerr*(mu*1000); %The standard errors are converted rather than the variances b/c standard errors have same units as impedance
[d.rho,d.pha,d.rhoerr,d.phaerr]=calc_rho_pha(d.Z,d.Zerr,d.T); % calculate from the converted impedance

%Check the the calculated apparent resistivities.
%If the off-diagonals are very large (e.g. 10^5) then it is possible that 
%the EDI was in field units and give a warning.
units_code = 1;
if nanmedian(d.rho(:,2:3,:))<1 %Possible that the EDI impedances were already in field units
    %This if statement checks if the median rho value is <1 Ohm m. It is assumed that if your entire
    %EDI has a median value <1 Ohm m, then something strange is going on. Either the units are
    %incorrect in the EDI, or you have some strange geology/static shift which should probably be
    %noted.
    disp(['WARNING: Site ', filename(1:end-4),' HAS VERY LOW RESISTIVITES (MEDIAN <1 OHM M). Check: (1) Does EDI file impedance have correct units? (2) Was site was processed correctly? (3) Large static shifts present?'])
    units_code = 0;
elseif nanmedian(d.rho(:,2:3,:))>10^6 %Possible that something else is wrong with your EDI involving units, processing or static shift
    disp(['WARNING: Site ', filename(1:end-4),' HAS VERY HIGH RESISTIVITES (MEDIAN >10^6 OHM M). Check: (1) Does EDI file impedance have correct units? (2) Was site was processed correctly? (3) Large static shifts present?'])
    units_code = 2;
end

end% END MAIN--------------------------------------------------------------




function [nfreq] = find_nfreq(tline)
[~,nfreq]=strtok(tline,'//');
ind=strfind(nfreq,'/');
nfreq(ind)=[];
nfreq=str2double(nfreq);

end
