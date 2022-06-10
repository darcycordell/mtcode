function write_data_edi(d,idx)
%
% This function is identical to write_edi_imp but takes a standardized data
% structure as input and re-arranged the data in this structure to be
% compatible with write_edi_imp formats
%
% Usage: write_data_edi(d,is)
%
% Inputs: "d" is a standard data structure
%       "is" is the OPTIONAL index of the site that you want to save
%           If no "is" is provided, then the script loops over all sites in
%           "d"
%
% If you want to write out all the data as EDI files then you need:
%
%       for is = 1:d.ns
%           write_data_edi(d,is)
%       end
%
% Note: the impedances in the data structure should *always* be in SI
% units. The resulting EDI file will have EDIs in field units (mV/km / nT).
% This is the same as the units of EDI files output from WinGLink
%
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

mu = 4*pi*10^-7;

%Convert impedances and impedance standard errors to field units.
d.Z = d.Z/(mu*1000);
d.Zerr = d.Zerr/(mu*1000);

zrot = d.zrot;
rhorot = d.zrot;
trot = d.trot;

if ~exist('idx','var')
    idx = 1:d.ns;
end

for is = idx
filename = [d.site{is},'.edi'];

[~, name, ~] = fileparts(filename);

Z(1,1,:,1) = d.Z(:,1,is); Z(1,2,:,1) = d.Z(:,2,is); Z(2,1,:,1) = d.Z(:,3,is); Z(2,2,:,1) = d.Z(:,4,is);
rho(1,1,:,1) = d.rho(:,1,is); rho(1,2,:,1) = d.rho(:,2,is); rho(2,1,:,1) = d.rho(:,3,is); rho(2,2,:,1) = d.rho(:,4,is);
pha(1,1,:,1) = d.pha(:,1,is); pha(1,2,:,1) = d.pha(:,2,is); pha(2,1,:,1) = d.pha(:,3,is); pha(2,2,:,1) = d.pha(:,4,is);
rhoerr(1,1,:,1) = d.rhoerr(:,1,is); rhoerr(1,2,:,1) = d.rhoerr(:,2,is); rhoerr(2,1,:,1) = d.rhoerr(:,3,is); rhoerr(2,2,:,1) = d.rhoerr(:,4,is);
phaerr(1,1,:,1) = d.phaerr(:,1,is); phaerr(1,2,:,1) = d.phaerr(:,2,is); phaerr(2,1,:,1) = d.phaerr(:,3,is); phaerr(2,2,:,1) = d.phaerr(:,4,is);
Zerr(1,1,:,1) = d.Zerr(:,1,is); Zerr(1,2,:,1) = d.Zerr(:,2,is); Zerr(2,1,:,1) = d.Zerr(:,3,is); Zerr(2,2,:,1) = d.Zerr(:,4,is);
tip(1,:,1) = d.tip(:,1,is); tip(2,:,1) = d.tip(:,2,is); 
tiperr(1,:,1) = d.tiperr(:,1,is); tiperr(2,:,1) = d.tiperr(:,2,is);

Zerr = real(Zerr).^2; %Convert standard errors to variance
tiperr = real(tiperr).^2; %Convert standard errors to variance

f = d.f;

if all(all(all(isnan(d.Z))))% if Z don't exist, then don't write them to file
    Z_flag=false;
else
    Z_flag=true;
end

if all(all(all(isnan(d.rho))))% if rho and phase don't exist, then don't write them to file
    rho_flag=false;
else
    rho_flag=true;
end

if all(all(all(isnan(d.tip))))% if tipper doesn't exist, then don't write it to file
    tip_flag=false;
else
    tip_flag=true;
end

% determine if zrot is different for each frequency, or if it is a single
% number, make it so that there is one angle per period
if length(zrot) ==1
    zrot=ones(length(f)).*zrot;
end
if length(rhorot) ==1
    rhorot=ones(length(f)).*rhorot;
end
if length(trot) ==1
    trot=ones(length(f)).*trot;
end


% change all NaN values to the empty value of 1e32
Z(isnan(Z))=1.0e32+(1i*1.0e32);
Zerr(isnan(Zerr))=1.0e32;
if tip_flag
    tip(isnan(tip))=1.0e32+(1i*1.0e32);
    tiperr(isnan(tiperr))=1.0e32;
end
if rho_flag
    rho(isnan(rho))=1.0e32;
    pha(isnan(pha))=1.0e32;
    rhoerr(isnan(rhoerr))=1.0e32;
    phaerr(isnan(phaerr))=1.0e32;
end


% Determine coordinates
lat = d.loc(is,1);
long = d.loc(is,2);
elev = d.loc(is,3);

lat=[fix(lat) abs(fix((lat-fix(lat))*60)) abs(((lat-fix(lat))*60-fix((lat-fix(lat))*60))*60)];
lat=[num2str(lat(1),'%02g') ':' num2str(lat(2),'%02g') ':' num2str(lat(3),'%06.3f')];

long=[fix(long) abs(fix((long-fix(long))*60)) abs(((long-fix(long))*60-fix((long-fix(long))*60))*60)];
long=[num2str(long(1),'%02g') ':' num2str(long(2),'%02g') ':' num2str(long(3),'%06.3f')];

elev=num2str(elev);


% Writing edi file
fid=fopen(filename,'wt');
fprintf(fid,'>HEAD\n');
fprintf(fid,'\n');
fprintf(fid,['DATAID="', strtok(name,'.') ,'"\n']);
fprintf(fid,'ACQBY=""\n');
fprintf(fid,'FILEBY=""\n');
fprintf(fid,'ACQDATE=""\n');
fprintf(fid,'FILEDATE=""\n');
fprintf(fid,'COUNTRY=""\n');
fprintf(fid,['LOC="' ,strtok(name,'.'), '"\n']);
fprintf(fid,['LAT=' lat '\n']);
fprintf(fid,['LONG=' long '\n']);
fprintf(fid,['ELEV=' elev '\n']);
fprintf(fid,'EMPTY=1.0e+32\n');
fprintf(fid,'\n');
fprintf(fid,'>INFO\n');
fprintf(fid,'\n');
fprintf(fid,'>=DEFINEMEAS\n');
fprintf(fid,'\n');
fprintf(fid,['REFLOC="' ,strtok(name,'.'), '"\n']);
fprintf(fid,['REFLAT=' lat '\n']);
fprintf(fid,['REFLONG=' long '\n']);
fprintf(fid,['REFELEV=' elev '\n']);
fprintf(fid,'\n');
fprintf(fid,'>HMEAS ID=1.0 CHTYPE=HX X=0 Y=0 AZM=0\n');
fprintf(fid,'>HMEAS ID=2.0 CHTYPE=HY X=0 Y=0 AZM=90\n');
fprintf(fid,'>HMEAS ID=3.0 CHTYPE=HZ X=0 Y=0 AZM=0\n');
fprintf(fid,'>EMEAS ID=4.0 CHTYPE=EX X=0 Y=0 X2=0 Y2=0\n');
fprintf(fid,'>EMEAS ID=5.0 CHTYPE=EY X=0 Y=0 X2=0 Y2=0\n');
fprintf(fid,'\n');
fprintf(fid,'>=MTSECT\n');
fprintf(fid,'\n');
fprintf(fid,['SECTID="' ,strtok(name,'.'), '"\n']);
fprintf(fid,['NFREQ=' num2str(length(f)) '\n']);
fprintf(fid,'\n');
if issorted(f); order='INC'; else order='DEC'; end
fprintf(fid,'%s\n','>!****FREQUENCIES****!');
fprintf(fid,['>FREQ NFREQ=' num2str(length(f)) ' ORDER=' order ' //' num2str(length(f)) '\n']);
for j=1:length(f)
    fprintf(fid,' %12.6e ',f(j));
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');

if Z_flag
fprintf(fid,'%s\n','>!****IMPEDANCE ROTATION ANGLES****!');
fprintf(fid,['>ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if zrot(j) >= 0
        fprintf(fid,' %12.6e ',zrot(j));
    else
        fprintf(fid,'%12.6e ',zrot(j));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%s\n','>!****IMPEDANCES****!');
fprintf(fid,['>ZXXR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if real(Z(1,1,j)) >= 0
        fprintf(fid,' %12.6e ',real(Z(1,1,j)));
    else
        fprintf(fid,'%12.6e ',real(Z(1,1,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZXXI ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if imag(Z(1,1,j)) >= 0
        fprintf(fid,' %12.6e ',imag(Z(1,1,j)));
    else
        fprintf(fid,'%12.6e ',imag(Z(1,1,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZXX.VAR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if Zerr(1,1,j)>= 0
        fprintf(fid,' %12.6e ',Zerr(1,1,j));
    else
        fprintf(fid,'%12.6e ',Zerr(1,1,j));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZXYR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if real(Z(1,2,j)) >= 0
        fprintf(fid,' %12.6e ',real(Z(1,2,j)));
    else
        fprintf(fid,'%12.6e ',real(Z(1,2,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZXYI ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if imag(Z(1,2,j)) >= 0
        fprintf(fid,' %12.6e ',imag(Z(1,2,j)));
    else
        fprintf(fid,'%12.6e ',imag(Z(1,2,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZXY.VAR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if Zerr(1,2,j) >= 0
        fprintf(fid,' %12.6e ',Zerr(1,2,j));
    else
        fprintf(fid,'%12.6e ',Zerr(1,2,j));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYXR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if real(Z(2,1,j)) >= 0
        fprintf(fid,' %12.6e ',real(Z(2,1,j)));
    else
        fprintf(fid,'%12.6e ',real(Z(2,1,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYXI ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if imag(Z(2,1,j)) >= 0
        fprintf(fid,' %12.6e ',imag(Z(2,1,j)));
    else
        fprintf(fid,'%12.6e ',imag(Z(2,1,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYX.VAR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if Zerr(2,1,j) >= 0
        fprintf(fid,' %12.6e ',Zerr(2,1,j));
    else
        fprintf(fid,'%12.6e ',Zerr(2,1,j));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYYR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if real(Z(2,2,j)) >= 0
        fprintf(fid,' %12.6e ',real(Z(2,2,j)));
    else
        fprintf(fid,'%12.6e ',real(Z(2,2,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYYI ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if imag(Z(2,2,j)) >= 0
        fprintf(fid,' %12.6e ',imag(Z(2,2,j)));
    else
        fprintf(fid,'%12.6e ',imag(Z(2,2,j)));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,['>ZYY.VAR ROT=ZROT //' num2str(length(f)) '\n']);
for j=1:length(f)
    if Zerr(2,2,j) >= 0
        fprintf(fid,' %12.6e ',Zerr(2,2,j));
    else
        fprintf(fid,'%12.6e ',Zerr(2,2,j));
    end
    if mod(j,5)==0 && j~=length(f)
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
end



if rho_flag
    fprintf(fid,['>RHOROT //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if zrot(j) >= 0
            fprintf(fid,' %12.6e ',rhorot(j));
        else
            fprintf(fid,'%12.6e ',rhorot(j));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    if rho_flag
        fprintf(fid,'%s\n','>!****APPARENT RESISTIVITIES AND PHASES****!');
        fprintf(fid,['>RHOXY ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rho(1,2,j));
            else
                fprintf(fid,'%12.6e ',rho(1,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOXY.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rhoerr(1,2,j));
            else
                fprintf(fid,'%12.6e ',rhoerr(1,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOYX ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rho(2,1,j));
            else
                fprintf(fid,'%12.6e ',rho(2,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOYX.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rhoerr(2,1,j));
            else
                fprintf(fid,'%12.6e ',rhoerr(2,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOXX ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rho(1,1,j));
            else
                fprintf(fid,'%12.6e ',rho(1,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOXX.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rhoerr(1,1,j));
            else
                fprintf(fid,'%12.6e ',rhoerr(1,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOYY ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rho(2,2,j));
            else
                fprintf(fid,'%12.6e ',rho(2,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>RHOYY.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',rhoerr(2,2,j));
            else
                fprintf(fid,'%12.6e ',rhoerr(2,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSXY ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',pha(1,2,j));
            else
                fprintf(fid,'%12.6e ',pha(1,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSXY.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',phaerr(1,2,j));
            else
                fprintf(fid,'%12.6e ',phaerr(1,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSYX ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',pha(2,1,j));
            else
                fprintf(fid,'%12.6e ',pha(2,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSYX.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',phaerr(2,1,j));
            else
                fprintf(fid,'%12.6e ',phaerr(2,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSXX ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',pha(1,1,j));
            else
                fprintf(fid,'%12.6e ',pha(1,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSXX.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',phaerr(1,1,j));
            else
                fprintf(fid,'%12.6e ',phaerr(1,1,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSYY ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',pha(2,2,j));
            else
                fprintf(fid,'%12.6e ',pha(2,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,['>PHSYY.ERR ROT=RHOROT //' num2str(length(f)) '\n']);
        for j=1:length(f)
            if rho(1,2,j) >= 0
                fprintf(fid,' %12.6e ',phaerr(2,2,j));
            else
                fprintf(fid,'%12.6e ',phaerr(2,2,j));
            end
            if mod(j,5)==0 && j~=length(f)
                fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');
    end
end




if tip_flag
    fprintf(fid,'%s\n','>!****TIPPER PARAMETERS****!');
    fprintf(fid,['>TROT.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if real(tip(1,j)) >= 0
            fprintf(fid,' %12.6e ',trot(j));  %MJU May 2012- removed first dim
        else
            fprintf(fid,'%12.6e ',trot(j));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TXR.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if real(tip(1,j)) >= 0
            fprintf(fid,' %12.6e ',real(tip(1,j)));  %MJU May 2012- removed first dim
        else
            fprintf(fid,'%12.6e ',real(tip(1,j)));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TXI.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if imag(tip(1,j)) >= 0
            fprintf(fid,' %12.6e ',imag(tip(1,j)));  %MJU May 2012- removed first dim
        else
            fprintf(fid,'%12.6e ',imag(tip(1,j)));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TXVAR.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if tiperr(2,j) >= 0
            fprintf(fid,' %12.6e ',tiperr(1,j));   %MJU May 2012 - removed first dim: DC April 2016 fix so that Tx var is being written correctly
        else
            fprintf(fid,'%12.6e ',tiperr(1,j));   %DC April 2016 - re-added 1 dim
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TYR.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if real(tip(2,j)) >= 0
            fprintf(fid,' %12.6e ',real(tip(2,j)));%MJU May 2012 - removed first dim
        else
            fprintf(fid,'%12.6e ',real(tip(2,j)));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TYI.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if imag(tip(2,j)) >= 0
            fprintf(fid,' %12.6e ',imag(tip(2,j))); %MJU May 2012 - removed first dim
        else
            fprintf(fid,'%12.6e ',imag(tip(2,j)));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,['>TYVAR.EXP //' num2str(length(f)) '\n']);
    for j=1:length(f)
        if tiperr(2,j) >= 0
            fprintf(fid,' %12.6e ',tiperr(2,j));  %MJU May 2012 - removed first dim
        else
            fprintf(fid,'%12.6e ',tiperr(2,j));
        end
        if mod(j,5)==0 && j~=length(f)
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'>END');
fclose(fid);

end

end
