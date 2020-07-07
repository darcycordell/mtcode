function mj_pcolor_v11
% Greg Nieuwenhuis - Jan 2013
%   - no longer requires a seperate station file, instead it makes its own
%   station file named after the first edi file read in (edifile.txt)


% Greg Nieuwenhuis - April 2012, mj_pcolor_v10
%       - station file (*.txt) is required containing 4 colummns:
%           1. name of station (must be the same as the .dcmp filename
%           2. longitude
%           3. latitude
%           4. elevation
%       - coordinates are forced into positive degrees east, not sure if
%           this works for all datasets
%       - got rid of coordinates in default file (retrieves them from the
%           location.txt file)
%       - BE WARNED - values of 0 in any of the parameters (including strike)
%           are removed and made into NaN values
%       - channeling angle (and strike) are rotated into the 0->90 degree
%           quadrant - is this correct?
%       - default file is now an .m file, not a .mat file
%       - everything of consequence is done in functions, not main program
% Dependant functions (external to Matlab):
%   - geo2utm, utm2geo (Dennis Rippe)
%   - m_map suite of functions ('http://www.eos.ubc.ca/~rich/map.html', updated Dec 2011)
%   - functions written for mtplot (Martyn Unsworth)
%       - calc_MT
%       - read_edi_dr (version which also outputs rotation angles)
%       - rot_z_rippe (version which rotates tipper as well)
%       - write_edi_rippe (version which writes rotation angles)
%       - calc_Z
%   - functions written specifically for mj_pcolor (located within this file):
%       - read_station_file
%       - make_default
%       - read_dcmp
%       - re_order_and_plot_data
%       - plot_pcolor
%       - plot_profile
%       - plot_pcolor_maps
%       - plot_strike_map
%       - plot_rose
%       - decomp_compare_multi
%------------------------------------------------------------------------
% Greg Nieuwenhuis, Feb 2011, mj_pcolor_v8, maps plotted in geographical
% coordinates, UTM coordinates converted using geo2utm (written by Dennis
% Rippe), frequency matrix indexed properly
%   - map_flag controls whether maps are plotted in kilometers or
%   geographical coords (careful changing this - it only works for some
%   functions)
%   - uses m_map toolbox when plotting in geographical coords
%   - pcolor maps of average shear, channelling, and twist can be
%   generated, with *.xyz output for plotting in GMT

%Ted Bertrand, mj_pcolor_v4, April 2008
%-Uses custom colormaps for pcolor plots. Updated figure titles
%-Implemented code to flip and rotate rose diagrams, 0deg is at top, +ve cw
%-Reading station coordinates and filenames is now independent of # of
% characters, only require station names to be all same length

%Ersan Turkoglu, mj_pcolor_v3, September 2007
%-Rose diagrams were somehow made using preferred strike before. Now rose
%diagrams are independent of preferred strike direction.

%Ersan Turkoglu, mj_pcolor_v2, June 2007,
%-Projected maps, can handle with different data sets, misfit map,
%defaults file

%Ersan Turkoglu, mj_pcolor_v1,December 2005

%Martyn Unsworth & Wolfgang Soyer, mj_pcolor, 2003
%---------------------------------------------------------------------------

clear all;
close all;

per_tol=5; % percentage tolerance allowed between period elements
% (if the difference is smaller than this, both periods
%  are combined)
profile_folder='./profiles/';%folder in which profile figures are saved
pcolor_folder='./pcolor/';
map_folder='./maps/';
rose_folder='./rose/';
hist_folder='./histograms/';

% read in station coordinates from station file, make a station file in
% case it is interesting for other applications, also check if dcmp file
% associated with edi file exists
menu_stations = menu('Which stations to plot?','Choose station file','Read all dcmp files in folder');
if menu_stations == 1
    [slocation,stn]=read_station_file;
else
    [slocation,stn]=make_station_file;
end


%looks for defaults_mjp.m file
if exist('defaults_mjp.m','file') == 2 % if default file exists
    defaults_mjp; % read in default file
    disp(['Default file read, Preferred strike direction set to: ',num2str(p_s)]) %#ok<NODEF>
    disp(['Project Name is: ',project_name])
    
else % if default file does not exist
    % make a new default file
    [p_s,plot_degrees, strike_scale, plot_names, nsect, pclims_strike, ...
        pclims_shear, pclims_channel, pclims_twist, pclims_misfit, pclims_skew, ...
        tlim,rlim, plim,project_name]...
        = make_default;
    
end

buffer=sqrt( (max(slocation(:,1))-min(slocation(:,1)))^2 + ( max(slocation(:,1))-min(slocation(:,1)) )^2 )*.1; % buffer is 10% of length of profile/grid

% limits are defined based on min and max coordinates of data with a buffer
limits=[min(slocation(:,1))-buffer max(slocation(:,1))+buffer...
    min(slocation(:,2))-buffer max(slocation(:,2))+buffer];


% if limits longitudes are not positive degrees east, make them
if limits(1) < 0
    limits(1)=limits(1)+360;
end
if limits(2) < 0
    limits(2)=limits(2)+360;
end

% these define the origin used to convert lat/long -> kilometers:
cent_long=limits(1)+((limits(2)-limits(1))/2);% the definitions of these are
orig_lat=limits(3)+((limits(4)-limits(3))/2);%  found in 'geo2utm.m'
%  currently defined as center
%  of station coords

% projection definition for m_map:
m_proj('mercator','long',limits(1:2),'lat',limits(3:4));%m_gshhs_l('patch',[1 1 1]);


if ~plot_degrees % if plots are in kilometers
    [limits(1), limits(3)]=geo2utm(limits(1), limits(3), cent_long, orig_lat); % change limits to kilomters
    [limits(2), limits(4)]=geo2utm(limits(2), limits(4), cent_long, orig_lat);
    limits=limits./1000; %meters to kilometers
    limits(1)=limits(1)-500; %UTM conversion includes 500 kilometer offset
    limits(2)=limits(2)-500;
    
    for is=1:length(slocation(:,1)) %change all coordinates to kilometers
        [slocation(is,1), slocation(is,2)]=geo2utm(slocation(is,1), slocation(is,2), cent_long, orig_lat);
        slocation(is,:)=slocation(is,:)./1000;
        slocation(is,1)=slocation(is,1)-500;
    end
end

% change preferred strike direction into radians
pref_strike = p_s*pi/180;

% Read dcmp files
[long, lat, per, param, nsites] = read_dcmp(slocation, stn);

% Output summary of strike, twist, shear and aniso to text file    % MJU 2016-03025
  str=['GBsummary.dat'];
  fid1=fopen(str,'w+');
  for is =1:nsites
    fprintf(fid1,'%9.3f %9.3f %9.3f %9.3f \n',param(1,is,1),param(4,is,1),param(2,is,1),param(11,is,1));  
  end
  fclose(fid1);
  
  % Plot the receovered impedances. Uses function similar to MTplot.  MJU 2016-3-25
  plot_decomp_data(17,long,lat,per,param,nsites)  
  


% reorder and project the dataset along preferred strike, bin periods, and
% plot the result
[per, param, xdist, long, lat, x, y, stn ] = re_order_and_plot_data(per, param, ...
    long, lat, stn, per_tol, plot_degrees, cent_long, orig_lat, pref_strike);

param(find(param==0))=NaN;  %Be careful! a strike of 0 will happen if the strike program did not work properly

%===========================================
% Plot
%===========================================
choice1 = 1;
while choice1 ~= 7
    menu1 = {'pcolor','profile','map','rose','Decomp Compare','histogram','exit'};
    menu2 = {'strike','shear','channelling','twist','misfit','skew'};
    
    choice1 = menu(' kind of plot ',menu1);
    if choice1 ~= 5 && choice1 ~= 7 % not for decomp compare or exit
        choice2 = menu(' parameter ',menu2);
        param_name=menu2{choice2};
        switch choice2
            case 5
                choice2 = 9; %because misfit is the 9th parameter in .dcmp file
                menu2{9}='misfit';
            case 6
                choice2 = 10; %because skew is the 10th parameter in .dcmp file
                menu2{10}='skew';
        end
        
    elseif choice1 == 5 % decomp compare
        choice2 = 7;
    else % exit
        choice2 = 1;
    end
    clear parplot %make sure parplot is refreshed for each new plot
    parplot = squeeze(param(choice2,:,:)); 
    nper = size(parplot,2);
    
    [misfit_tot]=calc_misfit_average(squeeze(param(9,:,:)));
    
    switch choice2
        case 1  % strike
            parplot = rem(rem(parplot,360)+360,90); %puts angle in the right quadrant
            pclims=pclims_strike;
            pylabel='degrees';
            cmap=[1 0 0;1 0.1111 0;1 0.2222 0;1 0.3333 0;1 0.4444 0;1 0.5556 0;1 0.6667 0;1 0.7778 0;1 0.8889 0;1 1 0;0.8571 1 0;0.7143 1 0;0.5714 1 0;0.4286 1 0;0.2857 1 0;0.1429 1 0;0 1 0;0.1429 1 0;0.2857 1 0;0.4286 1 0;0.5714 1 0;0.7143 1 0;0.8571 1 0;1 1 0;1 0.8889 0;1 0.7778 0;1 0.6667 0;1 0.5556 0;1 0.4444 0;1 0.3333 0;1 0.2222 0;1 0.1111 0;1 0 0];
            cmap_misfit=[0 1 0;0.0625 1 0;0.1250 1 0;0.1875 1 0;0.2500 1 0;0.3125 1 0;0.3750 1 0;0.4375 1 0;0.5000 1 0;0.5625 1 0;0.6250 1 0;0.6875 1 0;0.7500 1 0;0.8125 1 0;0.8750 1 0;0.9375 1 0;1 1 0;1 0.9375 0;1 0.8750 0;1 0.8125 0;1 0.7500 0;1 0.6875 0;1 0.6250 0;1 0.5625 0;1 0.5000 0;1 0.4375 0;1 0.3750 0;1 0.3125 0;1 0.2500 0;1 0.1875 0;1 0.1250 0;1 0.0625 0;1 0 0];
            %cmap=hsv(32);
        case 2 % Shear
            pclims=pclims_shear;
            pylabel='degrees';
            cmap=[1 0 0;1 0.1111 0;1 0.2222 0;1 0.3333 0;1 0.4444 0;1 0.5556 0;1 0.6667 0;1 0.7778 0;1 0.8889 0;1 1 0;0.8571 1 0;0.7143 1 0;0.5714 1 0;0.4286 1 0;0.2857 1 0;0.1429 1 0;0 1 0;0.1429 1 0;0.2857 1 0;0.4286 1 0;0.5714 1 0;0.7143 1 0;0.8571 1 0;1 1 0;1 0.8889 0;1 0.7778 0;1 0.6667 0;1 0.5556 0;1 0.4444 0;1 0.3333 0;1 0.2222 0;1 0.1111 0;1 0 0];
            %cmap=jet(32);
        case 3 % Channeling
            parplot = rem(rem(parplot,360)+360,90); %puts angle in the right quadrant
            pclims= pclims_channel;
            pylabel='degrees';
            cmap=[0 1 0;0.0625 1 0;0.1250 1 0;0.1875 1 0;0.2500 1 0;0.3125 1 0;0.3750 1 0;0.4375 1 0;0.5000 1 0;0.5625 1 0;0.6250 1 0;0.6875 1 0;0.7500 1 0;0.8125 1 0;0.8750 1 0;0.9375 1 0;1 1 0;1 0.9375 0;1 0.8750 0;1 0.8125 0;1 0.7500 0;1 0.6875 0;1 0.6250 0;1 0.5625 0;1 0.5000 0;1 0.4375 0;1 0.3750 0;1 0.3125 0;1 0.2500 0;1 0.1875 0;1 0.1250 0;1 0.0625 0;1 0 0];
            %cmap=jet(33)
        case 4 % Twist
            pclims= pclims_twist;
            pylabel='degrees';
            cmap=[1 0 0;1 0.125 0;1 0.25 0;1 0.375 0;1 0.5 0;1 0.625 0;1 0.75 0;1 0.875 0;1 1 0;0.875 1 0;0.75 1 0;0.625 1 0;0.5 1 0;0.375 1 0;0.25 1 0;0.125 1 0;0 1 0;0.125 1 0;0.25 1 0;0.375 1 0;0.5 1 0;0.625 1 0;0.75 1 0;0.875 1 0;1 1 0;1 0.875 0;1 0.75 0;1 0.625 0;1 0.5 0;1 0.375 0;1 0.25 0;1 0.125 0;1 0 0];
            %cmap=[1 0 0;1 0.1111 0;1 0.2222 0;1 0.3333 0;1 0.4444 0;1 0.5556 0;1 0.6667 0;1 0.7778 0;1 0.8889 0;1 1 0;0.8333 1 0;0.6667 1 0;0.5000 1 0;0.3333 1 0;0.1667 1 0;0 1 0;0.1667 1 0;0.3333 1 0;0.5000 1 0;0.6667 1 0;0.8333 1 0;1 1 0;1 0.9091 0;1 0.8182 0;1 0.7273 0;1 0.6364 0;1 0.5455 0;1 0.4545 0;1 0.3636 0;1 0.2727 0;1 0.1818 0;1 0.0909 0;1 0 0];
            %cmap=jet(32);
        case 9 % misfit
            pclims=pclims_misfit;
            pylabel=' ';
            cmap=[0 1 0;0.0625 1 0;0.1250 1 0;0.1875 1 0;0.2500 1 0;0.3125 1 0;0.3750 1 0;0.4375 1 0;0.5000 1 0;0.5625 1 0;0.6250 1 0;0.6875 1 0;0.7500 1 0;0.8125 1 0;0.8750 1 0;0.9375 1 0;1 1 0;1 0.9375 0;1 0.8750 0;1 0.8125 0;1 0.7500 0;1 0.6875 0;1 0.6250 0;1 0.5625 0;1 0.5000 0;1 0.4375 0;1 0.3750 0;1 0.3125 0;1 0.2500 0;1 0.1875 0;1 0.1250 0;1 0.0625 0;1 0 0];
            %cmap=jet(33);
        case 10 % skew
            parplot = tan(parplot*pi/180);
            pclims=pclims_skew;
            pylabel=' ';
            cmap=[0 1 0;0.1250 1 0;0.2500 1 0;0.3750 1 0;0.5000 1 0;0.6250 1 0;0.7500 1 0;0.8750 1 0;1 1 0;1 0.9286 0;1 0.8571 0;1 0.7857 0;1 0.7143 0;1 0.6429 0;1 0.5714 0;1 0.5000 0;1 0.4286 0;1 0.3571 0;1 0.2857 0;1 0.2143 0;1 0.1429 0;1 0.0714 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];
            %cmap=jet(33);
    end
    if choice1 ==1
        % =================> PCOLOR
        figure
        set(gcf,'Name',[pwd,' - pcolor']);
        
        plot_pcolor(xdist,parplot, per, cmap, pclims, stn);
        
        title([project_name,' ',char(menu2(choice2)),' Misfit= ',num2str(misfit_tot,'%.2f')]);
        set(gcf, 'paperPositionMode', 'auto');
        
        print_figures(pcolor_folder, param_name,project_name, min(min(per)), max(max(per)))
        
    elseif choice1 == 2 || choice1 == 3 || choice1 == 4
        prompt = {'choose period'};
        titles  = 'Choose period (range)';
        def = {[num2str(1),':',num2str(nper)]};
        ipers_string = char(inputdlg(prompt,titles,1,def));
        ipers=str2num(ipers_string);
        
        if (choice1==2)
            % =================> PROFILE
            figure
            plot_profile(parplot, ipers, xdist)
            
            xlabel(['km along preferred strike: ',num2str(pref_strike*180/pi)]);
            ylabel(pylabel);
            maxipers=max(per(1,ipers));
            minipers=min(per(1,ipers));
            title([param_name,' - ',project_name,' - T=',num2str(minipers),'s-',num2str(maxipers),'s']);
            set(gca,'XLim',[min(xdist)-(buffer*111);max(xdist)+(buffer*111)],'YLim',pclims);
            
            print_figures(profile_folder, param_name,project_name, minipers, maxipers)
            
            
        elseif (choice1==3 & choice2~=1)
            %==================> Contour Maps of shear, channeling, Skew, and twist
            
            figure
            [save_dat] = plot_pcolor_maps(parplot, per, ipers, stn, limits, long, lat, plot_degrees, plot_names);
            % save_dat is only used when saving xyz files for gmt plotting
            % (see the 'print figures' section below)
            
            maxipers=max(per(1,ipers));
            minipers=min(per(1,ipers));
            title([char(menu2(choice2)),' - ',project_name,' - T=',num2str(minipers),'s-',num2str(maxipers),'s']);
            caxis(pclims);
            colormap(cmap); colorbar;
            
            print_figures(map_folder, param_name,project_name, minipers, maxipers)
            
        elseif (choice1==3 & choice2==1)
            % =================> MAP of strike angles
            
            strikemap_choice = menu('How to map the strike angles?',{'red and blue crosses','single vector colored with misfit'});
            
            % this new function plots the strike vectors color coded to
            % misfit, with a length proportional to phase split
            figure_handle=figure;
            m_grid('box','fancy','tickdir','in'); hold on
            
            %plot_geoboundaries(figure_handle,'AB_geo.txt')
            
            if strikemap_choice==1
                % plot_strike_map is to plot as crosses:
                plot_strike_map(parplot, ipers, p_s, plot_degrees, strike_scale, x, y, cent_long, orig_lat);
                maptype='_cross';
            else
                % plot_strike_map_color is to plot as single lines with length
                % corresponding to phase difference, color corresponding to
                % misfit, 90° ambiguity removed based on preferred strike
                % orientation
                plot_strike_map_color(figure_handle,parplot, ipers, p_s, plot_degrees, strike_scale/10, x, y, cent_long, orig_lat,squeeze(param(9,:,:)),pclims_misfit,cmap_misfit,squeeze(param(12,:,:)));
                maptype='_colored';
            end
            hold off;
            
            maxipers=max(per(1,ipers));
            minipers=min(per(1,ipers));
            title(['Regional strike - ',project_name,' - T=',num2str(minipers),'s-',num2str(maxipers),'s']);
            
            print_figures(map_folder, [param_name,maptype],project_name, minipers, maxipers)
            
        elseif (choice1==4 & (choice2==1|choice2==3))
            % =================> ROSE   with either strike or
            % channeling
            
            prompt = {'choose site range'};
            stitle  = 'Choose site range';
            def = {['1:',num2str(nsites)]};
            isites = char(inputdlg(prompt,stitle,1,def));
            strike=parplot*pi/180;
            %%%%%%%%%%%%%%%%%%%%%
            parplot=squeeze((strike(:,ipers)));
            parplot = parplot(str2num(isites),:);
            
            figure
            plot_rose(parplot, nsect);
            
            maxipers=max(per(1,ipers));
            minipers=min(per(1,ipers));
            title([param_name,' ',project_name,' - T=',num2str(minipers),'s-',num2str(maxipers),'s']);
            
            print_figures(rose_folder, param_name,project_name, minipers, maxipers)
            
        end
    elseif choice1 == 5 
        % ===================> DECOMP COMPARE

        decomp_compare_multi(per, param,stn, per_tol, tlim, rlim, plim,pclims_misfit)
        
    elseif choice1 == 6
        % ====================> HISTOGRAM
        figure
        
        plot_hist(parplot);
        
        title([project_name,' ',char(menu2(choice2)),' Misfit= ',num2str(misfit_tot,'%.2f')]);
        xlabel('Angle (degrees)')
        ylabel('Number of Stations')
        
        print_figures(hist_folder, param_name,project_name, min(min(per)), max(max(per)))
        
    end
end


end % END MAIN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [long, lat,per, param, nsites] = read_dcmp(slocation, stn)
% Extract decomp information from .dcmp files
% Ersan Turgoklu(?) with edits by Greg Nieuwenhuis - April 2012
% IN:
%   slocation - [number of stations x 3] - column 1 - long
%                                                 2 - lat
%                                                 3 - elevation
%   stn - cell array of station names - [number of stations x 1]
% OUT:
%   long - [1 x number of stations]
%   lat - [1 x number of stations]
%   per - [number of stations x number of periods] - rows=stations,
%                                                   columns=periods
%   param - [12 x number of stations x number of periods] - 12 elements are
%                                       12 parameters read from .dcmp file
%               1 - regional strike
%               2 - shear
%               3 - channelling
%               4 - twist
%               5 - app rho a
%               6 - app rho b
%               7 - phase a
%               8 - phase b
%               9 - rms
%               10 - skew
%               11 - anis
%               12 - phadif
%   nsites - [1x1] - number of stations
%--------------------------------------------------------------------



%===========================================
% Site loop
%===========================================
for stn_ind=1:length(slocation(:,1))
    %===========================================
    % DCMP - Files
    %===========================================
    fname = [stn{stn_ind},'.dcmp'];
    fid=fopen(fname,'r');
    long(stn_ind)=slocation(stn_ind,1);
    lat(stn_ind)=slocation(stn_ind,2);
    %===========================================
    % Extract bands
    %===========================================
    line=1;
    while line ~= -1;
        line=fgetl(fid);
        k=findstr(line,'Per');% line 9 of *.dcmp (header line for band, min-per, max-per, start and end)
        if ~isempty(k)
            while line ~= -1
                line=fgetl(fid); % line 10 of *.dcmp
                if (length(line)>2)
                    ibound = sscanf(line(3:6),'%d'); % find the number of bounds (usually ~ line 10 of *.dcmp file)
                    bound(ibound,1) = sscanf(line(8:20),'%f'); % this is the minimum period used for the decomposition
                    bound(ibound,2) = sscanf(line(21:36),'%f'); % maximum period used for the decomposition
                else
                    line = -1;
                end
            end
        end
    end
    nbounds = ibound; % total number of bands used (line 10 of *.dcmp)
    %          Azimuth    Shear   Channeling    Twist   misfit skew
    sparam = {'regional','shear','channelling','twist','rho a','rho b','phase a','phase b','rms','skew','anis','phadif'};
    %===========================================
    % Extract parameters
    %===========================================
    for iextract=1:12
        line=1;
        while line ~= -1;
            line=fgetl(fid);
            if findstr(line,char(sparam{iextract}))
                nfreq(stn_ind)=sscanf(line,'%f',1); % number of frequencies in *.dcmp file
                for ifreq=1:nfreq(stn_ind)
                    line=fgetl(fid);
                    temp2=sscanf(line,'%f'); % contains an array of each number along that particular line
                    per_temp(ifreq)=temp2(1); % period at that line
                    per(stn_ind,ifreq)=temp2(1); % the period matrix
                    param(iextract,stn_ind,ifreq)=temp2(2);% one of the 12 parameters
                    
                end
                line=-1;
            end
        end
    end
    fclose(fid);
end   % sites loop

nsites=stn_ind;


end %END read_dcmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function plot_pcolor(xdist,parplot, per, cmap, pclims, stn)

nsites=length(per(:,1));
nper=length(per(1,:));

tx=xdist;
dx=diff(tx);

for ind_s=1:nsites %interpolate over NaN values so that pcolor image does not have holes
    not_nan=find(~isnan(parplot(ind_s,:)));
    if length(not_nan) > 1
        row=interp1(per(ind_s,not_nan),parplot(ind_s,not_nan),per(ind_s,:)); % linear interpolation over NaN points
        parplot(ind_s,:)=row;
    else
        display([stn{ind_s},' contains less than 2 datapoints, cannot plot on pcolor image'])
    end
end

% extend arrays
parplot(nsites+1,:)=parplot(nsites,:);
parplot(:,nper+1)=parplot(:,nper);
per_plot=per(1,:);
per_plot(nper+1)=per(1,nper);
tx(nsites+1)=tx(nsites)+dx(nsites-1)/2;

%shift tx to the left so stations are at the center of the
%columns:
tx(1)=tx(1)-dx(1)/2;
for i=1:length(dx)-1
    if dx(i+1) < dx(i)
        tx(i+1)=tx(i+1)-dx(i+1)/2;
    else
        tx(i+1)=tx(i+1)-dx(i)/2;
    end
end
tx(nsites)=tx(nsites)-dx(nsites-1)/2;

pcolor(tx,log10(per_plot),parplot');
hold on
plot(xdist',log10(per_plot(1))*ones(nsites,1),'kv');
set(gca,'YDir','reverse');
shading flat
axis tight
colormap(cmap)
caxis(pclims);
colorbar('vert')
xlabel('km'); ylabel('log period');


end %END plot_pcolor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_profile(parplot, ipers, xdist)

parplot_temp=parplot(:,ipers);
clear parplot
parplot=parplot_temp;
clear parplot_temp
for is=1:length(parplot(:,1)) % loop over the number of stations
    ind_nan=find(~isnan(parplot(is,:)));
    if ~isempty(ind_nan) % if no data existed within ipers at this station, make the average NaN
        parplot_temp(is,1)=mean(parplot(is,ind_nan),2);
    else
        parplot_temp(is,1)=NaN;
    end
    parplot_avg(is)=mean(parplot_temp(is,:));
end
%parplot=squeeze(mean(parplot(:,str2num(ipers)),2));
plot(xdist,parplot_avg,'+');hold on
ind=find(~isnan(parplot_avg));
avg=mean(parplot_avg(ind));
plot(xdist,ones(1,length(xdist)).*avg,'b--');hold off


end %end plot_profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_s,plot_degrees, strike_scale, nsect, plot_names, pclims_strike, pclims_shear, pclims_channel, pclims_twist, pclims_misfit, pclims_skew, tlim, rlim, plim,project_name]...
    = make_default
% function to make a default .m file containing all the plotting parameters
% for mj_pcolor


prompt = {'Preferred strike direction:',...
    'Plot maps in degrees (enter 1) or kilometers (enter 0)',...
    'Scale for strike vectors on maps',...
    'Number of sections in the rose diagrams',...
    'Plot names on pcolor maps (enter 1) or not (enter 0)',...
    'Limits for strike plots',...
    'Limits for shear plots',...
    'Limits for channeling plots',...
    'Limits for twist plots',...
    'Limits for misfit plots',...
    'Limits for skew plots',...
    'Period Limits for soundings',...
    'Apparent Resistivity Limits for soundings',...
    'Phase limits for soundings',...
    'Project Name for Figure titles',...
    };
dlg_title = 'Enter default values';
num_lines = 1;
def = {'30','1','2','32','1','[0,90]','[-45,45]','[0,90]','[-60,60]','[0,5]','[0,1]','[10^-0.3 10^4.3]','[10^-0 10^3]','[-10 100]','Default Project Name'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
p_s=str2num(char(answer{1}));
plot_degrees=str2num(char(answer{2}));
strike_scale=str2num(char(answer{3}));
nsect=str2num(char(answer{4}));
plot_names=str2num(char(answer{5}));
pclims_strike=str2num(char(answer{6}));
pclims_shear=str2num(char(answer{7}));
pclims_channel=str2num(char(answer{8}));
pclims_twist=str2num(char(answer{9}));
pclims_misfit=str2num(char(answer{10}));
pclims_skew=str2num(char(answer{11}));
tlim=str2num(char(answer{12}));
rlim=str2num(char(answer{13}));
plim=str2num(char(answer{14}));
project_name=char(answer{15});

fid=fopen('defaults_mjp.m','w'); % now save the variables into an m-file
fprintf(fid,'%s \t \t %s \n',['p_s=',num2str(p_s),';'],['%preferred strike direction (degrees)']);
fprintf(fid,'%s \t \t %s \n',['plot_degrees=',num2str(plot_degrees),';'],['%true to plot in degrees, false to plot in kilometers']);
fprintf(fid,'%s \t \t %s \n',['strike_scale=',num2str(strike_scale),';'],['% scale of strike vectors']);
fprintf(fid,'%s \t \t %s \n',['nsect=',num2str(nsect),';'],['%number of sections in the rose diagram']);
fprintf(fid,'%s \t \t %s \n',['plot_names=',num2str(plot_names),';'],['%Plot names on pcolor maps (enter 1) or not (enter 0)']);
fprintf(fid,'%s \t \t %s \n',['pclims_strike=[',num2str(pclims_strike),'];'],['%Limits for strike plots']);
fprintf(fid,'%s \t \t %s \n',['pclims_shear=[',num2str(pclims_shear),'];'],['%Limits for shear plots']);
fprintf(fid,'%s \t \t %s \n',['pclims_channel=[',num2str(pclims_channel),'];'],['%Limits for channeling plots']);
fprintf(fid,'%s \t \t %s \n',['pclims_twist=[',num2str(pclims_twist),'];'],['%Limits for twist plots']);
fprintf(fid,'%s \t \t %s \n',['pclims_misfit=[',num2str(pclims_misfit),'];'],['%Limits for misfit plots']);
fprintf(fid,'%s \t \t %s \n',['pclims_skew=[',num2str(pclims_skew),'];'],['%Limits for skew plots']);
fprintf(fid,'%s \t \t %s \n',['tlim=[',num2str(tlim),'];'],['%Period limits for sounding curves']);
fprintf(fid,'%s \t \t %s \n',['rlim=[',num2str(rlim),'];'],['%Apparent Resistivity limits for sounding curves']);
fprintf(fid,'%s \t \t %s \n',['plim=[',num2str(plim),'];'],['%Phase limits for sounding curves']);
fprintf(fid,'%s \t \t %s \n',['project_name=''',project_name,''';'],['%Project Name for Figure titles']);
fclose(fid);

end % END make_default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [per_out, param_out, xdist, long_out, lat_out, x_out, y_out, stn_out ] = re_order_and_plot_data(per, param, long, lat, stn, per_tol, plot_degrees, cent_long, orig_lat, pref_strike)
% function which:
%   - looks into the periods of each station and bins them to
%       within 10% of one another.
%   - re-orders each station along a profile which is projected ACROSS
%       strike
%   - finally plots the binned and ordered period array so that you can see
%       where data exists and where it does not
%   - zeros in the param_out matrix represent data points which do not
%       exist
% IN:
%   per - [number of stations x number of periods] - period range of data
%   param - [12 x number of stations x number of periods] - parameters from
%       the .dcmp files (strike, shear, etc)
%   long - [1 x number of stations] - longitude of stations
%   lat - [1 x number of stations] - latitude of stations
%   stn - cell array - [number of stations x 1] - names of stations
%   per_tol - [1 x 1] - tolerance level between periods for binning
%   plot_degrees - [1 x 1] - 1= plot in degrees, 0=plot in kilometers
%   cent_long - [1 x 1] - central longitude for UTM conversion into km
%   orig_lat - [1 x 1] - origin latitude for UTM conversion into km
%   pref_strike - [1 x 1] - preferred striek direction in radians
% OUT:
%   per_out - [number of stations x number of periods] - period range of data
%   param_out - [12 x number of stations x number of periods] - parameters from
%       the .dcmp files (strike, shear, etc)
%   xdist - [1 x number of stations] - distance across preferred strike
%       with respect to westernmost station
%   long_out - [1 x number of stations] - longitude of stations
%   lat_out - [1 x number of stations] - latitude of stations
%   x_out - [1 x number of stations] - longitude of stations in kilometers
%   y_out - [1 x number of stations] - latitude of stations in kilometers
%   stn_out - cell array - [number of stations x 1] - names of stations
%--------------------------------------------------------------------

nsites=length(long); %number of stations

% reorder and bin the period array so that each column contains similar
% periods
ind_zero=find(per~=0); % find where zeroes do not exist
per_array=per(ind_zero); % make one array with all the period numbers inside
per_array=unique(per_array);
differ=diff(per_array);
for j=1:length(differ)
    if abs(differ(j)/per_array(j))*100 <= per_tol % if the difference is less than the tolerance level
        per_array(j+1)=per_array(j); % make them the same
        differ=diff(per_array); %recalculate the difference now that elements have changed
    end
end
per_array=unique(per_array); % this is the period array which we will use from now on

per_new=zeros(nsites,length(per_array));
param_new=zeros(12,nsites,length(per_array));


num_per=length(param(1,1,:));

%re-order the per and param matrices:
for ind_s=1:nsites %loop over stations of the old array
    for ind_f=1:num_per%loop over periods of the old array
        if per(ind_s,ind_f) ~= 0 %only enter this loop if there is data at this period
            ind_move=find((abs(per_array-per(ind_s,ind_f))./per_array).*100 <= per_tol); %find the indice of the new per_aray within the tolerance allowed
            if length(ind_move) > 1
                % find the smallest difference between the values found within
                % the tolerance
                ind_small=find( abs(per_array(ind_move)-per(ind_s,ind_f))== min([abs(per_array(ind_move)-per(ind_s,ind_f))]) );
                temp=ind_move(ind_small);
                clear ind_move ind_small
                ind_move=temp;
            end
            if isempty(ind_move) % if no indicie was found, find the smallest difference in period and put the values there
                ind_move=find((abs(per_array-per(ind_s,ind_f))./per_array) == min(abs(per_array-per(ind_s,ind_f))./per_array));
                display(['Warning - ',num2str(per(ind_s,ind_f),'%.2f'),' Hz was put into the bin: ',num2str(per_array(ind_move),'%.2f'),' Hz, (difference of ',num2str((abs(per_array(ind_move)-per(ind_s,ind_f))/per_array(ind_move))*100,'%.1f'),'%)']);
            end
            per_new(ind_s,ind_move)=per(ind_s,ind_f);
            param_new(:,ind_s,ind_move)=param(:,ind_s,ind_f);
        end
    end
end

clear per param

if plot_degrees %if mapping units are in degrees
    [x,y]=geo2utm(long, lat, cent_long, orig_lat); % in meters
    x=(x./1000)-500; % meters to kilometers, and subtract 500km following UTM conversion
    y=y./1000;
else %otherwise if mapping units are in kilometers
    x=long;
    y=lat;
end

% project stations on line of preferred strike (sp)
beta = atan2(y,x); % angle of each plane with respect to measurement coordinate system
r = sqrt(x.^2 + y.^2); % distance of each point from the origin
xdist = r.*cos(pref_strike+beta); % adding pref_strike is like rotating the coordinate frame into the pref-strike frame, and r*cos projects the point onto the x axis
xdist=xdist-min(xdist);

% sort the stations so that the profile is in the correct order
[xdist,ind_sort]=sort(xdist);% ind_sort contains the ordered indices

for is=1:nsites
    per_out(is,:)=per_array;% make per_out to be all period values duplicated for each station
end
per_new=per_new(ind_sort,:);
param_out=param_new(:,ind_sort,:);%these are the newly sorted arrays
lat_out=lat(ind_sort); % zeros in the param matrix represent no data available
long_out=long(ind_sort);
x_out=x(ind_sort);
y_out=y(ind_sort);
stn_order=cell(length(x),1);
for j=1:length(x)
    stn_order{j}=stn{ind_sort(j)};
end
stn_out=cell(length(x),1);
stn_out=stn_order;

% now plot the period array so that you can see where data exists, and
% where it doesn't:
figure(72)
per_plot=per_array;
per_plot(length(per_plot)+1)=per_plot(length(per_plot)); % extend the array
ind_data=find(per_new == 0);
plot_matrix=ones(size(per_new));
plot_matrix(ind_data)=0;
plot_matrix(length(plot_matrix(:,1))+1,:)=plot_matrix(length(plot_matrix(:,1)),:);
plot_matrix(:,length(plot_matrix(1,:))+1)=plot_matrix(:,length(plot_matrix(1,:)));
pcolor([1:1:nsites+1], log10(per_plot),plot_matrix')
axis ij
cmap=colormap('gray');
caxis([0 1])
colormap(flipud(cmap));
title('Period Matrix; Black=data, White=no data')
xlabel('Stations')
ylabel('log(Period)')

end % END re_order_and_plot_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [save_dat] = plot_pcolor_maps(parplot, per, ipers, stn, limits, long, lat, plot_degrees, plot_names)
% function to plot pcolor maps of parameters

parplot_temp=parplot(:,ipers);
clear parplot
parplot=parplot_temp;
num_stn=length(parplot(:,1));
NaN_count=0;
count=0;
for j=1:num_stn
    notnan=find(~isnan(parplot(j,:))); % index which are not NaN, and within period range specified (ipers)
    if isempty(notnan)
        mean_parplot(j,1)=NaN; % if a station has no skew (all NaN)
        NaN_count=NaN_count+1;
    else
        count=count+1;
        mean_parplot(j,1)=mean(parplot(j,notnan),2);
        good_row(count)=j; % keep track of indices which have at least one calculated angle
    end
end
display([num2str(NaN_count),' stations without data, they are not plotted'])
clear parplot
parplot=mean_parplot;
%parplot=squeeze(mean(parplot(:,str2num(ipers)),2));

X=linspace(limits(1),limits(2),256);
Y=linspace(limits(3),limits(4),256)';
Z=griddata(long(good_row),lat(good_row),parplot(good_row),X,Y);
if plot_degrees % if plots are in degrees
    m_grid('box','fancy','tickdir','in'); hold on
    h=m_pcolor(X,Y,Z);hold on;
    set(h,'edgecolor','none')
    m_plot(long(good_row),lat(good_row),'k.','MarkerSize',12);
    if plot_names %to plot station names on maps
        offset=(max(long)-min(long))/50;
        for j=1:length(long)
            m_text(long(j)+offset,lat(j)+offset,stn{j},'FontSize',7);hold on;
        end
    end
else   %if coordinates are in kilometers, or mapping units are in kilometers
    pcolor(X,Y,Z);hold on
    shading flat
    plot(long(good_row),lat(good_row),'k.','MarkerSize',12);
end
save_dat=[long(good_row)',lat(good_row)',parplot(good_row)];

end %END plot_pcolor_maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_strike_map(strike, ipers, p_s, plot_degrees, strike_scale, x, y, cent_long, orig_lat)
% function to plot the regional strike directions on a map
% IN:
%   - strike [number of stations x number of periods] - contains strike
%       angles
%   - ipers - string containing the periods to be plotted (i.e. '1:53')
%   - p_s - [1 x 1] - preferred strike direction in degrees
%   - plot_degrees - 1 = plot map in degrees, 0=plot map in kilometers
%   - strike_scale [1 x 1] - scale of strike cross' on map (~1-10)
%   - x [1 x number of stations] - longitude coords in kilometers
%   - y [1 x number of stations] - latitude coords in kilometers
%   - cent_long - central longitude for UTM conversion (must be the same as
%       the value used to convert x and y to kilometers)
%   - orig_lat - origin latitude for UTM conversion (see above)

strike=strike.*(pi/180);%convert strike angles to radians

% rotate to align strikes to a preferred direction
sp = p_s*pi/180;%preferred strike in radians

aindex=find(strike-sp >= pi/4);
bindex=find(strike-sp < -pi/4);
if (sum(aindex)>0)
    strike(aindex)=strike(aindex)-pi/2;
end
if (sum(bindex)>0)
    strike(bindex)=strike(bindex)+pi/2;
end

%strike(:,str2num(ipers))*180/pi %displays the strike angle to screen
[num_stn,num_freq]=size(strike);
NaN_count=0;
count=0;
for is=1:num_stn
    notnan=find(~isnan(strike(is,ipers))); % index which are not NaN, and within period range specified (ipers)
    if isempty(notnan)
        %mean_strike(j,1)=pref_strike; % if a station has no strike direction (all NaN), replace with prefeered strike
        NaN_count=NaN_count+1;
    else
        count=count+1;
        mean_strike(is,1)=mean(strike(is,ipers(notnan)),2);
        good_row(count)=is; % keep track of indices which have at least one calculated strike direction
    end
    clear notnan
end
display([num2str(NaN_count),' stations without strike data'])

%strike=squeeze(mean(strike(:,str2num(ipers)),2)); %delete after +
clear strike
strike=mean_strike;
%                 create and plot bars

for is=1:length(strike)
    x_strike(is)=strike_scale.*cos(strike(is));
    y_strike(is)=strike_scale.*sin(strike(is));
    xs(:,is) = [x(is)-x_strike(is); x(is)+x_strike(is)];
    ys(:,is) = [y(is)+y_strike(is); y(is)-y_strike(is)];
    xs_pi2(:,is) = [x(is)-y_strike(is); x(is)+y_strike(is)];
    ys_pi2(:,is) = [y(is)-x_strike(is); y(is)+x_strike(is)];
end
if plot_degrees %if mapping units are in degrees
    m_grid('box','fancy','tickdir','in'); hold on
    [xs,ys]=utm2geo((xs.*1000)+500000, ys.*1000, cent_long, orig_lat);
    [xs_pi2,ys_pi2]=utm2geo((xs_pi2.*1000)+500000, ys_pi2.*1000, cent_long, orig_lat);
    m_plot(xs(:,good_row),ys(:,good_row),'r-'); hold on
    %m_plot(y,x,'k.','MarkerSize',10); hold on;
    m_plot(xs_pi2(:,good_row),ys_pi2(:,good_row),'b-')
    %xlabel('Longitude (\circE)');ylabel('Latitude (\circN)');
else
    plot(xs(:,good_row),ys(:,good_row),'r-'); hold on
    plot(xs_pi2(:,good_row),ys_pi2(:,good_row),'b-')
    axis('equal')
    xlabel('km'); ylabel('km');
end
hold on



end % END plot_strike_map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_strike_map_color(figure_handle,strike, ipers, p_s, plot_degrees, strike_scale, x, y, cent_long, orig_lat, misfit,pclims,cmap_misfit,phase_diff)
% function to plot the regional strike directions on a map
% IN:
%   - strike [number of stations x number of periods] - contains strike
%       angles
%   - ipers - string containing the periods to be plotted (i.e. '1:53')
%   - p_s - [1 x 1] - preferred strike direction in degrees
%   - plot_degrees - 1 = plot map in degrees, 0=plot map in kilometers
%   - strike_scale [1 x 1] - scale of strike cross' on map (~1-10)
%   - x [1 x number of stations] - longitude coords in kilometers
%   - y [1 x number of stations] - latitude coords in kilometers
%   - cent_long - central longitude for UTM conversion (must be the same as
%       the value used to convert x and y to kilometers)
%   - orig_lat - origin latitude for UTM conversion (see above)

strike=strike.*(pi/180);%convert strike angles to radians

% rotate to align strikes to a preferred direction
sp = p_s*pi/180;%preferred strike in radians

aindex=find(strike-sp >= pi/4);
bindex=find(strike-sp < -pi/4);
if (sum(aindex)>0)
    strike(aindex)=strike(aindex)-pi/2;
end
if (sum(bindex)>0)
    strike(bindex)=strike(bindex)+pi/2;
end

%strike(:,str2num(ipers))*180/pi %displays the strike angle to screen
[num_stn,num_freq]=size(strike);
NaN_count=0;
count=0;
for is=1:num_stn
    notnan=find(~isnan(strike(is,ipers))); % index which are not NaN, and within period range specified (ipers)
    if isempty(notnan)
        %mean_strike(j,1)=pref_strike; % if a station has no strike direction (all NaN), replace with prefeered strike
        NaN_count=NaN_count+1;
    else
        count=count+1;
        mean_strike(is,1)=mean(strike(is,ipers(notnan)),2);
        mean_misfit(is,1)=mean(misfit(is,ipers(notnan)),2);
        mean_phase_diff(is,1)=mean(phase_diff(is,ipers(notnan)),2);
        good_row(count)=is; % keep track of indices which have at least one calculated strike direction
    end
    clear notnan
end
display([num2str(NaN_count),' stations without strike data'])

%strike=squeeze(mean(strike(:,str2num(ipers)),2)); %delete after +
clear strike misfit
strike=mean_strike;
misfit=mean_misfit;
%                 create and plot bars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is for coordinates to plot scale vector
indx=find(max(x)==x);
indy=find(min(y)==y);
x_scale_loc=x(indx)-(111*0.4);
y_scale_loc=y(indy)-(111*0.4);
x_scale_tmp=strike_scale*cos(0)*50;
%y_scale_tmp=strike_scale*sin(0)*50;
x_scale=[x_scale_loc-x_scale_tmp;x_scale_loc+x_scale_tmp];
y_scale=[y_scale_loc;y_scale_loc];
[long_scale,lat_scale]=utm2geo((x_scale.*1000)+500000, y_scale.*1000, cent_long, orig_lat);
%%%%%%%%%%%%%%%%%%%%%%%%%

for is=1:length(strike)
    x_strike(is)=strike_scale.*mean_phase_diff(is,1).*cos(strike(is));
    y_strike(is)=strike_scale.*mean_phase_diff(is,1).*sin(strike(is));
    xs(:,is) = [x(is)-x_strike(is); x(is)+x_strike(is)];
    ys(:,is) = [y(is)+y_strike(is); y(is)-y_strike(is)];
    xs_pi2(:,is) = [x(is)-y_strike(is); x(is)+y_strike(is)];
    ys_pi2(:,is) = [y(is)-x_strike(is); y(is)+x_strike(is)];
    
    %now find the rgb values to plot each vector color with misfit value
    [rgb(is,:)] = plot_color(pclims(1), pclims(2), misfit(is), cmap_misfit);
end

figure(figure_handle)
if plot_degrees %if mapping units are in degrees
    [xs,ys]=utm2geo((xs.*1000)+500000, ys.*1000, cent_long, orig_lat);
    [xs_pi2,ys_pi2]=utm2geo((xs_pi2.*1000)+500000, ys_pi2.*1000, cent_long, orig_lat);
    [long,lat]=utm2geo((x.*1000)+500000, y.*1000, cent_long, orig_lat);
    %m_plot(xs(:,good_row),ys(:,good_row),'r-'); hold on
    m_plot(long,lat,'k.','MarkerSize',5); hold on;
    for is=1:length(good_row)
        m_plot(xs_pi2(:,good_row(is)),ys_pi2(:,good_row(is)),'color',rgb(good_row(is),:),'LineWidth',2)
    end
    caxis(pclims);
    h=colorbar;
    set(get(h,'ylabel'),'string','Misfit')
    %xlabel('Longitude (\circE)');ylabel('Latitude (\circN)');
    
    m_plot(long_scale,lat_scale,'k-','LineWidth',1)
    m_text(long_scale(1),lat_scale(1)+0.1,'50° Phase Difference','FontSize',6,'FontName','Arial')
else
    plot(xs(:,good_row),ys(:,good_row),'r-'); hold on
    plot(xs_pi2(:,good_row),ys_pi2(:,good_row),'b-')
    axis('equal')
    xlabel('km'); ylabel('km');
end
hold on



end % END plot_strike_map_color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_rose(parplot, nsect)
% function to plot a rose diagram of channeling angle or strike
% IN:
%   - parplot - [number of stations x number of periods] - contains strike
%       or channeling angle to plot
%   - nsect - number of sections in the rose diagram

[ns,nf]=size(parplot);
%look for unique values at each station, so that we are not just counting 
%   every strike angle which has been duplicated because it is in the same period band
count=1;
for is=1:ns
    tmp=unique(parplot(is,:));
    indnan=find(isnan(tmp));
    tmp(indnan)=[];
    if ~isempty(tmp)
        strike_angle(count,:)= tmp;
        count=count+1;
    else
        disp(['station number ',num2str(is),' contains no strike data'])
    end
end
clear parplot
parplot=strike_angle;
[ns,nf]=size(parplot);
avg_strike=mean(mean(parplot))*(180/pi);
disp([' The Average Strike angle is ',num2str(avg_strike,'%2f'),' degrees'])

% now reshape the matrix and add points so that rose vectors go out from
% opposite quadrants
parplot(ns+1:2*ns,:)=parplot+pi;
parplot=reshape(parplot,1,[]);
ind_nan=find(isnan(parplot));
parplot(ind_nan)=[];
rose(parplot+pi/2,nsect)
set(findobj(gca,'Type','line'),'Color','r');
hold on
rose(parplot,nsect);
set(gca,'XDir','reverse')  %TB - Flips through vertical axis
view(90,90);               %TB - Rotates 90deg clockwise so 0 is at top



end %END plot_rose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function decomp_compare_multi(decomp_per, param, stn,per_tol, tlim, rlim, plim,pclims_misfit)
% - reads the edi file corresponding to each station, and comapres to the
%       decomposed apprent resistivity and phase
% - edi data is rotated into the regional strike direction, determined by
%       averaging the strike over all periods
% IN:
%   decomp_per [number of stations x number of periods] periods in dcmp
%       files (equally spaced matrix
%   param [12 x number of stations x number of periods] - parameters from
%       dcmp files
%   stn [number of stations x 1] cell array of site names
%   tlim - period limits for sounding curves
%   rlim - app res limits for curves
%   plim - phase limits for curves
% OUT:
%   all figures are printed to the folder ./decomp_plots/
%------------------------------------------------------------------------

[nstn,nper]=size(decomp_per);% initialize decomp matrices
decomp_res=zeros(2,2,nstn,nper);
decomp_phs=zeros(2,2,nstn,nper);

strike=squeeze(param(1,:,:));
misfit=squeeze(param(9,:,:));
decomp_res(1,2,:,:)=(param(5,:,:));%extract relevant data from param matrix
decomp_res(2,1,:,:)=(param(6,:,:));
decomp_phs(1,2,:,:)=(param(7,:,:));
decomp_phs(2,1,:,:)=(param(8,:,:));

% find the average strike angle at each station (and find average misfit, 
% shear, and twist as well)
strike_avg=zeros(nstn,1);
shear_avg=zeros(nstn,1);
twist_avg=zeros(nstn,1);
misfit_avg=zeros(nstn,1);
for is=1:nstn
    count=0;
    for iper=1:nper
        if ~isnan(strike(is,iper)) % only avergae data which exists
            strike_avg(is)=strike_avg(is)+strike(is,iper);
            shear_avg(is)=shear_avg(is)+param(2,is,iper);
            twist_avg(is)=twist_avg(is)+param(4,is,iper);
            misfit_avg(is)=param(9,is,iper)+misfit_avg(is); 
            count=count+1;
        end
    end
    strike_avg(is)=strike_avg(is)/count;
    shear_avg(is)=shear_avg(is)/count;
    twist_avg(is)=twist_avg(is)/count;
    misfit_avg(is)=misfit_avg(is)/count;

end

h=figure;
% begin loop to choose station
choice=4;
while choice ~= 3
    if choice ~=4 %choice will only ==4 when the loop is initially entered
        choice=menu('What Next?','Next Station','Previous Station','Exit');
    end
    switch choice
        case 1 %next station
            if is >= nstn
                disp('You have reached the final station')
                is=is-1;
            end
            is=is+1;
        case 2 %previous station
            if is < 1
                disp('You are at the first station, you cannot plot the previous station!')
                is=1;
            end
            is=is-1;
        case 3% exit
            break
        case 4 %force the next station to plot when loop is first entered
            choice =1;
            is=1;
    end
    if ~isnan(strike_avg(is)) % if the decomp file lacks any strike information,
                              %  it means that the strike program failed for that 
                              %  station
        % read edi file
        if exist([stn{is},'.edi'])
            
            % read the edi file using read_edi_rippe
            [Z,Zvar,Tip,Tvar,Rhoa,Rhoaerr,Phs,Phserr,f_edi,rot_edi,coords]= read_edi_imp([stn{is},'.edi']);
            T_edi=1./f_edi;
            nfreq_edi=length(f_edi);
            disp([stn{is},'.edi opened']);
            
            % rotate edi data into strike coordinates from dcmp file
            if rot_edi ~=0 %if the edi rotation is not initially zero
                
                %rotate the edi impedance data to zero
                [Z,Zvar]=rot_z_dr(Z,Zvar, nfreq_edi,-rot_edi);
                
                rot_edi=0;
            end
            
            %rotate the edi impedance data into the same coordinate system as the decomp data
            Z_calc=(4*pi*10^-4)*Z; 
            Zvar_calc=(4*pi*10^-4)*Zvar;
            [Z,Zvar]=rot_z_dr(Z_calc,Zvar_calc, nfreq_edi,strike_avg(is));
                        
            % calculate app res and phase from rotated impedance
            [Rhoa,Rhoaerr,Phs,Phserr]=calc_MT(Z,Zvar,f_edi);
            
            % ----------------------calculate impedance for edi writing
            
            % now find which period values corresponds between decomp matrix and edi file
            for iper=1:length(T_edi(:))
                ind_decomp=find((abs(decomp_per(is,:)-T_edi(iper))./decomp_per(is,:)).*100 <= per_tol);
                if length(ind_decomp) > 1
                    % find the smallest difference in ind_edi to assign
                    % period values
                    ind_small=find( abs(decomp_per(is,ind_decomp)-T_edi(iper))== min([abs(decomp_per(is,ind_decomp)-T_edi(iper))]) );
                    temp=ind_decomp(ind_small);
                    clear ind_edi ind_small
                    ind_decomp=temp;
                end
                if isempty(ind_decomp)
                    T_ind(iper)=NaN;
                else
                    T_ind(iper)=ind_decomp; %T_ind is the same length as T_edi, and each indice points to the corresponding period in decomp_per
                end
            end
            
            
            
            %-------------- Plot the results
            figure(h)
            set(gcf, 'Units','normalized')
            window=[0.25 0.2 0.5 0.75];
            set(gcf,'OuterPosition',window);
            
            
            subplot(8,1,[1:3])
            loglog(T_edi,squeeze(Rhoa(1,1,:)),'r:') ; hold on;
            loglog(T_edi,squeeze(Rhoa(1,2,:)),'r-') ;
            loglog(T_edi,squeeze(Rhoa(2,2,:)),'b:') ;
            loglog(T_edi,squeeze(Rhoa(2,1,:)),'b-') ;
            
            loglog(decomp_per(is,:),squeeze(decomp_res(1,2,is,:)),'r.--')
            loglog(decomp_per(is,:),squeeze(decomp_res(2,1,is,:)),'b.--')
            
            hold off
            axis([tlim(1) tlim(2) rlim(1) rlim(2)])
            title([stn{is}, '; Strike = ', num2str(strike_avg(is),'%.0f'),'^o','; Shear = ',num2str(shear_avg(is),'%.1f'),'; Twist = ',num2str(twist_avg(is),'%.1f'),'; RMS = ',num2str(misfit_avg(is),'%.1f')])
            ylabel ('\rho_a (\Omega m)')
            
            subhandle=subplot(8,1,[4:6]);
            window=get(subhandle,'Position');
            set(subhandle,'Position',[window(1) window(2) window(3) window(4)*0.8])% set the second subplot to be 80% of the height it would normally be (to fit the title in)
            semilogx(T_edi,squeeze(Phs(1,2,:)),'r-') ; hold on;
            semilogx(T_edi,squeeze(Phs(2,1,:)),'b-') ;
            
            semilogx(decomp_per(is,:),squeeze(decomp_phs(1,2,is,:)),'r.--') ;
            semilogx(decomp_per(is,:),squeeze(decomp_phs(2,1,is,:)),'b.--') ;
            
            axis([tlim(1) tlim(2) plim(1) plim(2)]); hold off;
            title('Red = xy; blue = yx; solid line = original; dotted line = decomposition, dashed line = xx and yy original data')
            ylabel ('phase (deg)')
            
            ind_notnan=find(~isnan(misfit(is,:)));
            subhandle=subplot(8,1,8);
            window=get(subhandle,'Position');
            set(subhandle,'Position',[window(1) window(2) window(3) window(4)*1.5])% set the second subplot to be 150% of the height it would normally be
            semilogx(decomp_per(is,ind_notnan), [1:1:length(decomp_per(is,ind_notnan))], 'white') % Use semilogx to set up axes but plot the line invisible
            bar(decomp_per(is,ind_notnan),misfit(is,ind_notnan),'hist');
            set(gca,'XScale','log')
            axis([tlim(1) tlim(2) pclims_misfit(1) pclims_misfit(2)])
            title('Misfit')
            xlabel ('Period (s)')
            
            %--------------------- END plotting
            set(gcf, 'Units','default')
            
            print_figures('./decomp_plots/', 'decomp_compare',stn{is}, min(decomp_per(is,:)), max(decomp_per(is,:)))
            
            clear Z Zvar Tip Tvar Rhoa Phs f_edi rot_edi coords
        else
            disp([stn{is},'.edi Does not exist'])
        end
    else
        clf
        title(['Strike determination for ',stn{is},' failed'])
    end
end



end % END decomp_compare_multi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function print_figures(save_folder, param_name,project_name, minipers, maxipers)
% opens a new directory and saves the figures inside

if ~exist(save_folder,'dir')%check if the folder already exists
    mkdir(save_folder);% make new folder for plots
end
oldfolder=cd(save_folder); %change directory to the new folder for printing
%print('-depsc',[char(menu2(choice2)),'-map.eps']);
print('-djpeg100',[param_name,'-',project_name,'-',num2str(minipers),'s-',num2str(maxipers),'s.jpg']);
print('-depsc','-painters',[param_name,'-',project_name,'-',num2str(minipers),'s-',num2str(maxipers),'s.eps']);
%saveas(gcf, [char(menu2(choice2)),'-map'], 'fig');
cd(oldfolder);%return to the old folder for the next station

end % END print_figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [misfit_tot]=calc_misfit_average(misfit_temp)
% calculates the average misfit while ignoring NaN values
[nx,ny]=size(misfit_temp);
misfit_tot=0;
count=0;
for ix=1:nx
    for iy=1:ny
        if ~isnan(misfit_temp(ix,iy))
            misfit_tot=misfit_tot+misfit_temp(ix,iy); %
            count=count+1;
        end
    end
end
misfit_tot=misfit_tot/count; %find average
end %END calc_misfit_average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rgb] = plot_color(minimum, maximum, num, cmap)
% function takes a number, the minimum, maximum, colormap, and returns the
% rgb values for that number, within the colormap, scaled to the min and
% max specified.
% The output rgb can be used directly in a plot function, for example:
   %plot(x,y,'color',rgb);

%cmap - [: x 3] - example jet()
%rgb - [1 x 3]

values=linspace(minimum,maximum,length(cmap(:,1))); %this is the array of numbers which matches the array of colors in cmap

ind=find(min(abs(values-num))==abs(values-num)); %this is the indice of the closest number to the input number

rgb=cmap(ind,:); %rgb is the RGB value associated with the input number in cmap, scaled to the minimum and maximum


end %plot_color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%=====================================================

function plot_hist(parplot)

for is=1:length(parplot(:,1)) % loop over the number of stations
    ind_nan=find(~isnan(parplot(is,:)));
    if ~isempty(ind_nan) % if no data existed within ipers at this station, make the average NaN
        parplot_temp(is,1)=mean(parplot(is,ind_nan),2);
    else
        parplot_temp(is,1)=NaN;
    end
    parplot_avg(is)=mean(parplot_temp(is,:));
end

%bins=[-60:5:60];
histfit(parplot_avg,50);
v=axis;
axis([v(1) v(2) 0 50])


end


%===================================================================
function [slocation,stn] = make_station_file
% This function will read in all the edi files in the working directory and
% make a station location file with the filename "fname"

% Reading in file names
files=dir('*.edi');
% Looping over all files
for i=1:length(files)
    % Determining site name
    station_in=files(i).name;
    % Reading in and writing edi file
    fid=fopen(station_in,'rt');
    count=0;
    while count < 3
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if strncmp(tline,'LAT=',4)
            lat_str = sscanf(tline(5:length(tline)),'%s');
            if isempty(strfind(lat_str,':')) %if the latitude is in decimal degree format
                files(i).lat=str2double(lat_str);
            else % if the latitude is in deg:min:sec format
                [deg_str,remain]=strtok(lat_str,':');
                deg=str2double(deg_str);
                if deg < 0 %if the degree is negative, then you must multiply the final answer by -1
                    deg=deg.*-1;
                    mult=-1;
                else %otherwise multiply it by +1
                    mult=1;
                end
                remain(1)=[];
                [min_str,remain]=strtok(remain,':');
                mini=str2double(min_str);
                remain(1)=[];
                sec=str2double(remain);
                files(i).lat=dms2deg(deg,mini,sec)*mult;
            end
            
            count=count+1;
        elseif strncmp(tline,'LONG=',5)
            long_str = sscanf(tline(6:length(tline)),'%s');
            if isempty(strfind(long_str,':')) %if the longitude is in decimal degree format
                files(i).long=str2double(long_str);
            else %if teh longitude is in deg:min:sec format
                [deg_str,remain]=strtok(long_str,':');
                deg=str2double(deg_str);
                if deg < 0 %if the degree is negative, then you must multiply the final answer by -1
                    deg=deg.*-1;
                    mult=-1;
                else %otherwise multiply it by +1
                    mult=1;
                end
                remain(1)=[];
                [min_str,remain]=strtok(remain,':');
                mini=str2double(min_str);
                remain(1)=[];
                sec=str2double(remain);
                files(i).long=dms2deg(deg,mini,sec)*mult;
            end
            count=count+1;
            
        elseif strncmp(tline,'ELEV=',5)
            elev_str = sscanf(tline(6:length(tline)),'%s');
            files(i).elev=str2double(elev_str);
            
            count=count+1;
        end
    end
    fclose(fid);
end% end reading edi files
%===============================================

%===============================================
%organize sloaction and stn variables for output from the function
for is=1:length(files)
    files(is)
    slocation(is,1)=files(is).long;
    slocation(is,2)=files(is).lat;
    slocation(is,3)=files(is).elev;
    stn{is}=strtok(files(is).name,'.');
end

%find the stations which exist as .dcmp files, and disregard all others
count=1;
for j=1:length(files)
    fname=[stn{j},'.dcmp'];
    if exist(fname,'file')~=2
        %         slocation(j,1)=NaN; % Use NaN as a marker for stations which do not exist
        disp(['site ',fname,' not found'])
    else
        stn_temp{count}=stn{j};
        slocation_temp(count,:)=slocation(j,:);
        count=count+1;
    end
end
clear slocation stn
slocation=slocation_temp; % long, lat & elev of station coords
stn=stn_temp; % cell array of station names (indices match slocation)


% make sure all longitude coords are in positive degrees east
for j=1:length(slocation(:,1))
    if slocation(j,1) < 0
        slocation(j,1)=slocation(j,1)+360;
    end
end
%=====================================

%=============================================
%write station file in case it is of use elsewhere:
fname=[strtok(files(1).name,'.'),'.txt'];
fid=fopen(fname,'w');
for is=1:length(files)
        
    fprintf(fid,'%s\t%3.8f\t%3.8f\t%4.3f\n',strtok(files(is).name,'.'),files(is).long,files(is).lat,files(is).elev);
    
end
fclose(fid);  

end % END make_station_file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [slocation,stn]=read_station_file

pause(0.1)
[fname,fpath,tmp]=uigetfile('*.txt');

if tmp ~= 1
    error('Station file not opened')
end


fid=fopen([fpath,fname]);

line=fgetl(fid);
count=1;
while line ~= -1
    stn{count}=sscanf(line,'%s',1);
    line(1:length(stn{count}))=[];
    slocation(count,:)=sscanf(line,'%f %f %f',3);
    
    count=count+1;
    line=fgetl(fid);
end

%find the stations which exist as .dcmp files, and disregard all others
count=1;
for j=1:length(slocation(:,1))
    fname=[stn{j},'.dcmp'];
    if exist(fname,'file')~=2
        %         slocation(j,1)=NaN; % Use NaN as a marker for stations which do not exist
        disp(['site ',fname,' not found'])
    else
        stn_temp{count}=stn{j};
        slocation_temp(count,:)=slocation(j,:);
        count=count+1;
    end
end
clear slocation stn
slocation=slocation_temp; % long, lat & elev of station coords
stn=stn_temp; % cell array of station names (indices match slocation)


% make sure all longitude coords are in positive degrees east
for j=1:length(slocation(:,1))
    if slocation(j,1) < 0
        slocation(j,1)=slocation(j,1)+360;
    end
end
%=====================================

end % end read_station_file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%