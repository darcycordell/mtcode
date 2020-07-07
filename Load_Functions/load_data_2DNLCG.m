function [dobs, dpred] = load_data_2DNLCG(par,dpred_file,sitefile)
% Function which loads Randy Mackie's 2D NLCG Inversion format data files
% for both the observed data and inversion response (i.e. predicted) data.
%
% For more info on the inversion:
%   Rodi, W., & Mackie, R. L. (2001). Nonlinear conjugate gradients algorithm 
%   for 2-D magnetotelluric inversion. Geophysics, 66(1), 174–187.
%
% Apr. 2020: Error floor from par file is now applied. The computed misfit
% is within 2 decimal places of the values in the 2D inversion .rsp file.
%
% The data structures include:
% Z = impedance in SI units
% Zerr = impedance error in SI units
% rho = apparent resistivity in Ohm m
% rhoerr = apparent resistivity error in Ohm m
% pha = phase in degrees (e^(-iwt) convention)
% phaerr = phase error in degrees
% tip = tipper
% tiperr = tipper error
% T = period in seconds (increasing periods)
% f = frequency in Hertz (decreasing frequency)
% rotation = data rotation in degrees (clockwise from north)
% site = MT site names
% x, y = MT site locations in x and y model coordinates (x is N-S)
% origin = The center of the model mesh in [long, lat]
% ns, nf, and nr are the number of stations, frequencies and responses
% components = the string of components which are included.
%
%
% To be added:----------------------------------------------------------
%
% - Add station elevations (z in model coordinates)
% - Reference the model to lat-long space somehow
% - Add loc (site locations in lat, long and elevation projected onto the profile)
% - Add the rotation angle (which the data has been rotated to)


dobs.origin = [0,0]; %How to add origin to reference center of mesh in lat/long?
dpred.origin = dobs.origin;

dobs.rotation = 0; %How to get rotation information of data from inversion files?
dpred.rotation = 0;

dobs.name = par.data_files{1}(1:end-4);
dpred.name = dpred_file;

disp('Loading 2D NLCG Data Files')

fid=fopen(dpred_file,'r');

line = 1; isite = 0;

while line ~= -1
   
    line = fgetl(fid);
    
    if strfind(line,'block')
        isite = isite+1; 
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        textscan(fgetl(fid),'%s');
        line = fgetl(fid);
        
        iper = 0;
        while 1
            line = fgetl(fid);
            if strcmp(line,' ')==1
                break
            else
                iper = iper+1;
                
                % replace ***** values in .rsp file
                line = regexprep(line,'\<[*][*]+\s\>','NaN ');
                    
                read_data = cell2mat(textscan(line,'%f')); 
                
                dobs.rho(iper,[2 3],isite) = read_data([5 1]);
                dobs.pha(iper,[2 3],isite) = read_data([6 2]);
                dobs.tip(iper,2,isite) = read_data(9)+1i*read_data(10);
                
                dpred.rho(iper,[2 3],isite) = read_data([7 3]);
                dpred.pha(iper,[2 3],isite) = read_data([8 4]);
                dpred.tip(iper,2,isite) = read_data(11)+1i*read_data(12);
                
                dobs.T(iper,1) = read_data(13);
                dpred.T(iper,1) = read_data(13);
            end

            
        end
       
    end
           
end

rm = dobs.rho==100 & dobs.pha==-45;
dobs.rho(rm) = NaN;
dobs.pha(rm) = NaN;
dobs.tip(dobs.tip==0) = NaN+1i*NaN; % check this - zero can be data


dobs.rho(:,[1,4],:) = NaN;
dobs.pha(:,[1,4],:) = NaN;
dpred.rho(:,[1,4],:) = NaN;
dpred.pha(:,[1,4],:) = NaN;
dobs.tip(:,1,:) = NaN+1i*NaN; % no Tzx in 2D inversion
dpred.tip(:,1,:) = NaN+1i*NaN;

dobs.pha(:,2,:) = -dobs.pha(:,2,:);
dobs.pha(:,3,:) = -(dobs.pha(:,3,:)+180);

dpred.pha(:,2,:) = -dpred.pha(:,2,:);
dpred.pha(:,3,:) = -(dpred.pha(:,3,:)+180);

% dobs.pha(:,2,:) = dobs.pha(:,2,:)+90;
% dobs.pha(:,3,:) = dobs.pha(:,3,:)-90;
% 
% dpred.pha(:,2,:) = dpred.pha(:,2,:)+90;
% dpred.pha(:,3,:) = dpred.pha(:,3,:)-90;


sze = size(dobs.rho);

dobs.f = 1./dobs.T;
dpred.f = 1./dpred.T;

dobs.ns = sze(3); dpred.ns = dobs.ns;
dobs.nr = 4; dpred.nr = 4;
dobs.nf = sze(1); dpred.nf = dobs.nf;


%----------------------Load Site Names-----------------------------------
[dobs.site, locxy, ~] = load_sitefile(sitefile,dobs.ns);
dpred.site = dobs.site;
dobs.x = locxy(:,1);
dobs.y = locxy(:,2);
dobs.z = locxy(:,3);

dpred.x = dobs.x;
dpred.y = dobs.y;
dpred.z = dobs.z;

dpred.niter = ''; % placeholder

%----------------------Load the Errors------------------------------------
for i = 1:length(par.mode)
    if par.mode(i) == 1
        [dobs] = load_err(dobs,dpred,par,i);
    end
end

%Calculate the impedances from the apparents resistivity and phase
[dobs.Z,dobs.Zerr] = calc_Z(dobs.rho,dobs.rhoerr,dobs.pha,dobs.phaerr,dobs.T);

dpred.rhoerr = NaN(sze);
dpred.phaerr = NaN(sze);
dpred.tiperr = NaN(sze(1),2,sze(3))+1i*NaN(sze(1),2,sze(3));

[dpred.Z,~] = calc_Z(dpred.rho,dpred.rhoerr,dpred.pha,dpred.phaerr,dpred.T);

dpred.Z(:,4,:) = NaN + 1i*NaN;
dpred.Zerr = NaN(sze) + 1i*NaN(sze);

if par.mode(1) == 1
    dobs.responses{1} = 'ZYX';
end

if par.mode(2) == 1
    dobs.responses{2} = 'ZXY';
end

if par.mode(3)
    dobs.responses{3} = 'TX';
end

dpred.responses = dobs.responses;

end
%%
function [dobs]=load_err(dobs,dpred,par,mode)
%Function which reads the observed data files and reads the observed error
%in those files. Error floor is applied.

fid=fopen(par.data_files{mode},'r');

if mode == 1 || mode == 2
    err_block = nan(dobs.nf,5,dobs.ns);
elseif mode == 3
    err_block = nan(dobs.nf,4,dobs.ns);
end

ns=str2double(fgetl(fid));

line = str2num(fgetl(fid));

nf = line(1);
j=1;
while 1
    while j<=ns
        i=1;
        while i<=nf
            line=fgetl(fid);
            line(regexp(line,'[(,)]'))=[];
            tmp = str2num(line); %#ok<*ST2NM>
            if tmp(1) ~= nf  
                ind = nearestpoint(tmp(1),dobs.T);
                err_block(ind,:,j) = tmp;
              
                i=i+1;
            end

        end
        j=j+1;
    end
    
    if (i-1)*(j-1) == ns*nf
        break
    end
    
end

if mode == 1 % TM
    
    tm_rho_err = err_block(:,4,:);
    tm_rho_err(tm_rho_err < par.tm_errflr(1)/100) = par.tm_errflr(1)/100;
    dobs.rhoerr(:,3,:) = tm_rho_err; % these are in logarithmic absolute error
    
    tm_pha_err = err_block(:,5,:);
    tm_pha_err(tm_pha_err < (par.tm_errflr(2)*0.5/100) ) = (par.tm_errflr(2)*0.5/100) ; % 0.5 comes from the fact that error in ln(rho) is 2x error in phase (rad)
    dobs.phaerr(:,3,:) = tm_pha_err; % these are in radians absolute error
        
elseif mode == 2 % TE

    te_rho_err = err_block(:,4,:);        
    te_rho_err(te_rho_err < par.te_errflr(1)/100) = par.te_errflr(1)/100;    
    dobs.rhoerr(:,2,:) = te_rho_err; % these are in logarithmic absolute error
    
    te_pha_err = err_block(:,5,:);
    te_pha_err(te_pha_err < (par.te_errflr(2)*0.5/100) ) = (par.te_errflr(2)*0.5/100) ; % 0.5 comes from the fact that error in ln(rho) is 2x error in phase (rad)
    dobs.phaerr(:,2,:) = te_pha_err; % these are in radians absolute error
    
elseif mode == 3 % tip
    
    tip_err = err_block(:,4,:);
    tip_err(tip_err < par.tip_errflr) = par.tip_errflr;
    dobs.tiperr(:,2,:) = tip_err+1i*tip_err;
      
end

dobs.rhoerr(:,[1 4],:) = NaN;
dobs.phaerr(:,[1 4],:) = NaN;
dobs.tiperr(:,1,:) = NaN+1i*NaN; % no Tzx in 2D inversion

dobs.rhoerr(dobs.rhoerr==10^6)=NaN;
dobs.phaerr(dobs.phaerr==10^6)=NaN;
dobs.tiperr(dobs.tiperr==10^6+1i*10^6) = NaN+1i*NaN;
dobs.tiperr(dobs.tiperr==0) = 10^6+1i*10^6;

%The rho error which is listed in the Mackie data files is error in ln(rho).
%This is an absolute error. The standard data structure uses relative error
%of apparent resistivity in Ohm m. The relative errors calculated here are 
% the errors that would result in the same misfit in ln(rho) space.
% Therefore, these are neither the errors from the data files nor the
% errors that were used in the inversion.

%The phaerr which is listed in the Mackie data file is in radians. Convert
%to degrees here

if mode == 1 % TM
    dobs.rhoerr(:,3,:) = tm_rho_err .* abs(dobs.rho(:,3,:)-dpred.rho(:,3,:)) ./ abs(log(dobs.rho(:,3,:))-log(dpred.rho(:,3,:))); % TM only
    dobs.phaerr(:,3,:) = (180/pi)*dobs.phaerr(:,3,:);
elseif mode == 2
    dobs.rhoerr(:,2,:) = te_rho_err .* abs(dobs.rho(:,2,:)-dpred.rho(:,2,:)) ./ abs(log(dobs.rho(:,2,:))-log(dpred.rho(:,2,:))); % TE only
    dobs.phaerr(:,2,:) = (180/pi)*dobs.phaerr(:,2,:);
end

if numel(dobs.tiperr)==1 % not sure if this is needed
    sze = size(dobs.rho);
    dobs.tiperr = NaN(sze(1),2,sze(3))+1i*NaN(sze(1),2,sze(3));
end

fclose(fid);

end %END load_err
