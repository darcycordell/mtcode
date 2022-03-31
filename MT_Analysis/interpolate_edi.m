function d = interpolate_edi(filename)
%
% Function which loads EDIs and interpolates those EDIs onto a common set
% of frequencies. This function was previously called "load_edis_v7". Most
% of the code is unchanged.
%
% The primary difference is that this function now outputs a matfile which
% contains a "d" MT data structure instead of outputting ZZZ, ZZZ_HZ, EEE,
% EEE_HZ, etc.
%
% Usage: d = interpolate_edi(filename)
%       
%   *Optional Input* = filename to save data structure mat file
%
% Option to bin data into frequency bins or interpolate data based on an
% input frequency text file. It is most common to use a frequency text
% file.
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

max_perc_diff=20; %largest percent difference for making new period array which is uniform for all stations. Easier to bin data with smaller value?

choice1 = menu('Choose frequencies to use','Bin periods from file','Interpolate periods from file','Bin periods determined from the data','Interpolate over periods determined from the data','Exit');
if choice1 == 1
    disp('Binning with user provided frequencies')
    pause(0.5);
    [freq_file,freq_path] = uigetfile({'*.txt'},'pick  file containing frequencies');
    ff = load( [freq_path freq_file], '-ascii');
    Tnew = 1./ff';
    Tnew = sort(Tnew); %Periods need to be ascending
    interp_flag=0;
elseif choice1==2
    disp('Interpolating with user provided frequencies, handy if there exists many data gaps, then this will fill data by interpolating')
    pause(0.5);
    [freq_file,freq_path] = uigetfile({'*.txt'},'pick  file containing frequencies');
    ff = load( [freq_path freq_file], '-ascii');
    Tnew = 1./ff';
    Tnew = sort(Tnew); %Periods need to be ascending
    interp_flag=1;
elseif choice1==3
    disp('Bin over a period set determined by the data')

    interp_flag=0;
elseif choice1 == 4 
    disp('Interpolate over a period set determined by the data')
    interp_flag = 1;
    
    % Note that by uncommenting the folowing line you can hardwire the
    % periods to a pre-defined set
    %TT=logspace(-4, 5, 72);  %IMPORTANT equally sampled periods (opton used in load_edis_v2 and v3
    % USE THE FOLLOWING VALUES WHEN WORKING WITH MOSTLY NIMS DATA (STILL
    % INTERPOLATES!)
%     TT= 1./[1.781260e+00  ,1.343740e+00  ,1.000000e+00  ,7.500020e-01  ,5.624990e-01  ,4.062500e-01,...
%         3.125000e-01,  2.343750e-01,  1.718750e-01,  1.328130e-01,  1.015630e-01,  7.812501e-02 ,...
%         5.859369e-02,  4.296870e-02,  3.203130e-02,  2.421880e-02,  1.796880e-02,  1.328120e-02 ,...
%         1.015620e-02,  7.617188e-03,  5.664059e-03,  4.296870e-03,  3.320311e-03,  2.539059e-03 ,...
%         1.757810e-03,  1.318360e-03,  9.765620e-04,  7.324219e-04,  5.371091e-04,  4.150390e-04 ,...
%         3.173829e-04,  2.441410e-04,  1.831050e-04,  1.342770e-04,  1.007080e-04,  7.629388e-05 ,...
%         5.493162e-05,  3.967291e-05 ];
else
    return
end


%-------------------GET EDI FILE NAMES IN CURRENT DIRECTORY---------------
directory=dir;
fff=length(directory);
count=1;
for iio=1:fff
    [~, ~, ext] = fileparts(directory(iio).name);
    if strcmp(ext,'.edi')
        edst{count}=directory(iio).name;
        count=count+1;
    else
        disp([directory(iio).name,' is not an edi file']);
    end 
end

nsta=length(edst); %nsta is the number of edi files in the directory, fff is the maximum number of characters in the edi filename


%-------------------------LOAD RAW EDI DATA--------------------------------
% start loop over sites
T_all = [];
for iyu=1:nsta
    % read edi file
    disp([mfilename,' : Reading ',edst{iyu}])
    [d(iyu),units_code(iyu)] = load_data_edi(edst{iyu});
    
    %All frequencies present in all EDIs
    T_all = unique([T_all; d(iyu).T]);
       
    for ifreq = 1:d(iyu).nf
        % rotate impedances to geographic north
        d(iyu).Z(ifreq,:,1) = rotate_Z(d(iyu).Z(ifreq,:,1),-d(iyu).zrot(ifreq));
        %rotate tipper to geographic north
        d(iyu).tip(ifreq,:,1) = rotate_tip(d(iyu).tip(ifreq,:,1),-d(iyu).trot(ifreq));
        % Remember that resistivity and phase are not rotated!
        [d(iyu).rho(ifreq,:), d(iyu).pha(ifreq,:), d(iyu).rhoerr(ifreq,:), d(iyu).phaerr(ifreq,:)] = calc_rho_pha(d(iyu).Z(ifreq,:,1),d(iyu).Zerr(ifreq,:,1),d(iyu).T(ifreq));
        
        d(iyu).zrot(ifreq) = 0;
        d(iyu).trot(ifreq) = 0;
    end

end


%---------------------PERFORM QC CHECKS------------------------------------
indmask = units_code==0;
ind = find(units_code==0);

if sum(indmask)>0
    disp('****************************************************************')
    disp([num2str(sum(indmask)),' site(s) appear to have VERY LOW resistivity values. This may be due to:'])
    disp('(1) a unit conversion issue; (2) a processing issue or; (3) static shifts/geology')
    fprintf('\n')
    disp('If (1), then use convert_edi_units.m to correct the units')
    fprintf('\n')
    disp('Please double check the following EDI files:')
    for i = 1:length(ind)
        disp([d(ind(i)).site{1},'.edi'])
    end
    disp('****************************************************************')
end

indmask = units_code==2;
ind = find(units_code==2);

if sum(indmask)>0
    disp('****************************************************************')
    disp([num2str(sum(indmask)),' site(s) appear to have VERY HIGH resistivity values. This may be due to:'])
    disp('(1) a unit conversion issue; (2) a processing issue or; (3) static shifts/geology')
    fprintf('\n')
    disp('Please double check the following EDI files:')
    for i = 1:length(ind)
        disp([d(ind(i)).site{1},'.edi'])
    end
    disp('****************************************************************')
end

%%

if nansum([d(:).loc])==0
    disp('WARNING: It appears that you do not have any latitude and longitude information for your survey. You will need to add it manually')
end


%%
%----------------GET NEW FREQUENCY/PERIOD DATA SET-------------------------

if ~exist('Tnew','var') % make a period array based on the data, and the period tolerance defined at the top of the function
    Tnew=find_new_periods(d,max_perc_diff);
    Tnew = sort(Tnew); %Periods need to be ascending
end

%Determine if the "Tnew" frequencies to be interpolated/binned onto goes
%beyond the range of frequencies present in the EDIs.
Tnew(Tnew>max(T_all))=[];
Tnew(Tnew<min(T_all))=[];


%--------------------INTERPOLATE / BIN------------------------------------
for iyu = 1:nsta
    
    if interp_flag==0  %Binning:
        [dnew(iyu)] = bin_impedance(d(iyu),Tnew);
    else % Interpolation (spline)
        [dnew(iyu)] = interpolate_impedance(d(iyu),Tnew);
    end
   
end


%------------------SET UP DATA STRUCTURE----------------------------------
d_orig = d;
clear d;

d.T = dnew(1).T';
d.f = 1./d.T;
d.nf = length(d.T);
d.ns = nsta;
d.nr = 4;
d.responses = d_orig(1).responses;
        

for is = 1:nsta
    
    d.Z(:,:,is) = dnew(is).Z;
    d.Zerr(:,:,is) = dnew(is).Zerr;
    d.tip(:,:,is) = dnew(is).tip;
    d.tiperr(:,:,is) = dnew(is).tiperr;
    
    d.site{is,1} = d_orig(is).site{1};
    d.loc(is,:) = d_orig(is).loc;
    
    d.zrot(:,is) = zeros(d.nf,1);
    d.trot(:,is) = zeros(d.nf,1);
    
end

[d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);

d = set_map_projection(d);
d.origin = [(max(d.loc(:,1))-min(d.loc(:,1)))/2+min(d.loc(:,1)),(max(d.loc(:,2))-min(d.loc(:,2)))/2+min(d.loc(:,2)),0];
d.niter = '';

if exist('filename','var')
    d.name = filename;
else
    d.name = ['Data_',edst{1}(1:3)];
end


save(d.name,'d');

choice=menu('Plot the interpolation results?','Yes','No');
if choice == 1
    % plot data:
    is = 1;
    while 1
        
        set_figure_size(1);
        d.Zerr = nan(size(d.Zerr));
        plot_impedance(d_orig(is),1);
        plot_impedance(d,is);

        
        next_menu = menu('','Next','Previous','Done');
    
        if next_menu == 1
            is = is+1;
            close all
            if is>d.ns
                is = d.ns;
            end
        elseif next_menu == 2
            is = is-1;
            if is<1
                is = 1;
            end
            close all
        else
            close all
            break
        end
        
        
    end
end

end %end of main

%==========================================================================
function [TT] = find_new_periods(d,per_tol)
%==========================================================================

nsta=length(d);

% First step is to find a set of periods which matches the data as good as
% possible.
%put all the periods from every station into one long array
count=1;
for iyu=1:nsta
    nf=d(iyu).nf;
    freq_array(count:count+nf-1)=d(iyu).f;
    count=count+nf;
end
freq_array=unique(freq_array);% remove duplicates
differ=diff(freq_array); % find the differences between each
for j=1:length(differ)
    if abs(differ(j)/freq_array(j))*100 <= per_tol % if the difference is less than the tolerance level
        new_freq=mean([freq_array(j),freq_array(j+1)]);
        freq_array(j+1)=new_freq; % make them the same mean value
        freq_array(j)=new_freq;
        differ=diff(freq_array); %recalculate the difference now that elements have changed
    end
end
freq_array=unique(freq_array); % this is the period array which we will use from now on

TT=1./freq_array;

TT=sort(TT);% sort in ascending order


end

%==========================================================================
function [dnew] = interpolate_impedance(d,Tnew)
%=========================================================================
% this function interpolates the impedance z (and tipper z_hz) from the
% irregularily sampled 'period', to the regularily sampled array 'TT'
% same goes for the errors

%Added ability for the script to exclude large gaps in the data from the interpolation.
%If the EDI has been edited and there is a gap greater than one decade
%(e.g. a gap with log10(T2) - log10(T1) > 1), then this gap is excluded
%from the new set of periods (Tnew)
T_orig = Tnew;

dT = diff(log10(Tnew)); %Find difference in logarithmic periods
nindx = nearestpoint(1,cumsum(dT)); %Find the difference in indices for one logarithmic period


%%
nindx = 10;
Z = nan(length(T_orig),4,d.ns);
Z = Z + 1i*Z;
Zerr = nan(size(Z));
Zerr = Zerr + 1i*Zerr;
tip = nan(length(T_orig),2,d.ns);
tip = tip + 1i*tip;
tiperr = nan(size(tip));
tiperr = tiperr + 1i*tiperr;

%Do the interpolation
for is = 1:d.ns
    for ic=1:4 %this loop interpolates the impedance(z) and error(dz) data onto the 'TT' array (for all 4 impedance elements, and 2 tipper elements)
        ind=find(~isnan(real(d.Z(:,ic,is)))); %Find the nan indices
        Tnew = T_orig; %Reset the new period vector to the original
        indgaps = find(diff(ind)>nindx); %Find any gaps in the data > 1 decade
        if ~isempty(indgaps)
            %Find the indices in the new period vector which correspond to
            %the gaps
            ind1 = nearestpoint(d.T(ind(indgaps)),Tnew);
            ind2 = nearestpoint(d.T(ind(indgaps+1)),Tnew);
            
            indout = [];
            for i = 1:length(ind1) %If there is more than one gap, then loop over gaps
                indout = [indout ind1(i):ind2(i)]; %Indices in gaps
            end
            indin = setdiff(1:length(Tnew),indout); %Indices not in gaps
            
            Tnew(indout) = []; %Remove the gaps from the new period vector (i.e. the interpolation query points)
        else
            indin = 1:length(Tnew);
        end
        
        if ~isempty(ind) && ~all(real(d.Z(:,ic,is))==0)
            dnew.Z(:,ic,is)=10.^interp1(log10(d.T(ind)),log10(real(d.Z(ind,ic,is))),log10(Tnew),'pchip',NaN)+1i*10.^interp1(log10(d.T(ind)),log10(imag(d.Z(ind,ic,is))),log10(Tnew),'pchip',NaN);
            clear ind
            ind=find(~isnan(real(d.Zerr(:,ic,is))));
            dnew.Zerr(:,ic,is)=10.^interp1(log10(d.T(ind)),log10(real(d.Zerr(ind,ic,is))),log10(Tnew),'linear',NaN); %error is not complex
            dnew.Zerr(:,ic,is) = dnew.Zerr(:,ic,is)+1i*dnew.Zerr(:,ic,is);
            clear ind
        elseif isempty(ind) && ~all(real(d.Z(:,ic,is))==0)% if ind is empty, it means there is no data in z, therefore we need to define ZZint and EEint as NaN
            dnew.Z(:,ic,is)=NaN(length(Tnew),1)+1i.*NaN(length(Tnew),1);
            dnew.Zerr(:,ic,is)=NaN(length(Tnew),1)+1i.*NaN(length(Tnew),1);
        elseif all(real(d.Z(:,ic,is))==0) %if all the impedances in this element are zero (can happen with synthetic data)
            dnew.Z(:,ic,is)=zeros(length(Tnew),1);
            dnew.Zerr(:,ic,is)=zeros(length(Tnew),1);
        end

        Z(indin,ic,is) = dnew.Z(:,ic); %Add the interpolated data into the holding variable at the correct indices
        Zerr(indin,ic,is) = dnew.Zerr(:,ic);
        
    end
end

for is = 1:d.ns
    for ic = 1:2
        ind=find(~isnan(real(d.tip(:,ic,is))));       
        Tnew = T_orig;
        indgaps = find(diff(ind)>nindx);
        if ~isempty(indgaps)
            ind1 = nearestpoint(d.T(ind(indgaps)),Tnew);
            ind2 = nearestpoint(d.T(ind(indgaps+1)),Tnew);
            
            indout = [];
            for i = 1:length(ind1)
                indout = [indout ind1(i):ind2(i)];
            end
            indin = setdiff(1:length(Tnew),indout);
            
            Tnew(indout) = [];
        else
            indin =1:length(Tnew);
        end
        
        if ~isempty(ind) && ~all(real(d.tip(:,ic,is))==0)
            dnew.tip(:,ic,is)=interp1(log10(d.T(ind)),real(d.tip(ind,ic,is)),log10(Tnew),'pchip',NaN)+1i*interp1(log10(d.T(ind)),imag(d.tip(ind,ic,is)),log10(Tnew),'pchip',NaN);
            clear ind
            ind=find(~isnan(real(d.tiperr(:,ic,is))));
            dnew.tiperr(:,ic,is)=interp1(log10(d.T(ind)),real(d.tiperr(ind,ic,is)),log10(Tnew),'linear',NaN); %error is not complex
            dnew.tiperr(:,ic,is) = dnew.tiperr(:,ic,is)+1i*dnew.tiperr(:,ic,is);
            clear ind
        elseif isempty(ind) && ~all(real(d.tip(:,ic,is))==0) % if ind is empty, it means there is no data in z_hz, therefore we need to define ZZint_hz and EEint_hz as NaN
            dnew.tip(:,ic,is)=NaN(length(Tnew),1)+1i.*NaN(length(Tnew),1);
            dnew.tiperr(:,ic,is)=NaN(length(Tnew),1)+1i.*NaN(length(Tnew),1);
        elseif all(real(d.tip(:,ic,is))==0)
            dnew.tip(:,ic,is)=zeros(length(Tnew),1);
            dnew.tiperr(:,ic,is)=zeros(length(Tnew),1);
        end

        tip(indin,ic,is) = dnew.tip(:,ic); %Add the interpolated data into the holding variable at the correct indices
        tiperr(indin,ic,is) = dnew.tiperr(:,ic);

    end
end


    dnew.T = T_orig; %Return the period vector to its original size

    %Replace variables in the data structure
    dnew.Z = Z;
    dnew.Zerr = Zerr;
    dnew.tip = tip;
    dnew.tiperr = tiperr;
    

end% end of function interpolate_impedance



%==========================================================================
function [dnew] = bin_impedance(d, Tnew)
%==========================================================================
% this function bins impedance(z), tipper(z_hz), and their errors(dz, dz_hz)
% from the measured 'period' onto the new TT array
% ALTERNATIVELY - one could use the interp1 function using the 'nearest'
% method (which is a nearest neighbour interpolant). This way was chosen
% simply because there is more control
dnew.T=Tnew;
perc_diff = zeros(d.nf,1);
for is = 1:d.ns
    dnew.Z(1:length(Tnew),1:4,is)=NaN; %start with a matrix of NaN's
    dnew.Zerr(1:length(Tnew),1:4,is)=NaN;
    dnew.tip(1:length(Tnew),1:2,is)=NaN;
    dnew.tiperr(1:length(Tnew),1:2,is)=NaN;
    for ifreq=1:d.nf
        ind_move=find( abs(d.T(ifreq)-Tnew(:)) == min(abs(d.T(ifreq)-Tnew(:))) ); %find the indice of the minimum difference between measured period, and TT
        perc_diff(ifreq)= (abs(Tnew(ind_move)-d.T(ifreq))/d.T(ifreq)) *100;
        if ~isnan(dnew.Z(ind_move,:,is)) % if another period has already been moved to this TT, then ZZint will no longer be NaN at ind_move ;
            error(['Two different measured periods are trying to be binned into the same period array, Tnew= ',num2str(Tnew(ind_move),'%.4f'),'s. Need to adjust period tolerance in interpolate_edi'])
        else    
            dnew.Z(ind_move,:)=d.Z(ifreq,:,is);
            dnew.Zerr(ind_move,:)=d.Zerr(ifreq,:,is);
            dnew.tip(ind_move,:)=d.tip(ifreq,:,is);
            dnew.tiperr(ind_move,:)=d.tiperr(ifreq,:,is);
        end
    end

    disp([num2str(max(perc_diff),'%.1f'),'% - largest difference in frequency arrays when binning the data'])

end

end % end of function bin_impedance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

