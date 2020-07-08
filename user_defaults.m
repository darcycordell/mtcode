function [u] = user_defaults
% for MTcode version 2020-06-10
%------------------User input structures-----------------------------------
%
% These are variables which are likely survey dependent and need to be
% changed by the user. I've grouped them all together into a "u" (for
% "user") structure. These could also made into an external text file or
% parameter file or something that each user has their own version of.
%
% Since these options are survey dependent, it may be useful to copy this
% file to the current directory for each survey and then edit it there.
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)


%GENERAL OPTIONS-----------------------------------------------------------

u.output_figure = {'png','eps'}; %Select file format to save figures.
                   % 'none' no figure is saved
                   % 'jpg' JPG is saved
                   % 'png' PNG is saved
                   % 'eps' EPS vector image is saved.
                   % Enter a cell array to print multiple types of files. For example, {'png','eps'} will print both png and eps files.
                   %EPS files are quite large so it is recommended to only
                   %use 'eps' option when producing figures you want to use
                   %in presentations or publications.
u.geofile = 'none'; %Default is 'none', put your geofile text file here if you want it.
u.projection = 'mercator'; %Projection used for map plotting
    %Options include: 'mercator', 'lambert conformal conic', 'Albers Equal-Area Conic'
        %'Stereographic', 'UTM', etc. See m_map
u.plot_topo = false; %Automatically plot topography on maps, isosurface plots, etc. 
                %Topo is not plotted on plots that already have color maps
                %(e.g. model slice plots, interpolated rho/pha plots, etc.)
                %Topo is automatically downloaded using download_srtm.m
u.topo_file = 'none'; %If u.plot_topo = true but no topo_file is specified, then the SRTM topography is automatically downloaded
                %using download_srtm.m. This can make plotting very slow.
                %It is preferable to download the required SRTM files only once using download_srtm. This
                %function automatically saves a file called "elevation.mat".
u.topo_transparency = 0.5; %Set transparency for topography plots (1 = opaque; 0 = transparent)
u.elev_colim = [0 3000]; %Set topography colorbar axis when plotting topography on a map

%MODEL PLOTTING INFO------------------------------------------------------------

%General Model Inputs

u.plot_input_model = false; %Plot the inversion starting model? True or false.
u.num = 'default'; %Default is 'default' for iteration. If your inversion has iteration information then it will put the
                    %inversion iteration number here.

%Colormap
% enter the name of the colormap to use in pcolor plots
% 'mtcode_default' is the default colormap used by MTcode, it is similar to matlab's 'jet'. this is defined below
% 'mtcode_red_blue' is a red-to-blue colormap
% can also type the name of any embedded MATLAB colormap, e.g. 'default', 'parula', 'jet', 'gray' ...
u.cmap_name = 'mtcode_default';

u.colim=[0 3];%resistivity colorbar limits for color plots (logarithmic)
u.clabel = 'linear'; % type of colorbar labels for resistivity plot. for example, log = 0,1,2,3 and linear = 1, 10, 100, 1000. 
                    %note that in both cases the colorbar is logarithmic; only the labels are changed

%Contours
u.plot_contours = true; %"true" to plot contours, "false" to not
u.contours = [-6 1 6]; %contours to plot on log10(resistivity) slices
u.contour_text = 'off'; %"on" to plot contour labels, "off" to not

%Axis limits
u.zmin = -4; %top of cross-section plots (km)
u.zmax = 15; %maximum depth to plot cross-sections (km)
%u.xylims = [-100 100 -500 500]; %x (NS) and y (EW) model coordinate limits to plot in kilometers
u.xylims = 'default';
        %Example: [-15 15 -20 20]; First two coordinates are NS
        %'default' plots default non-padding cells
u.maplims = 'default'; %latitude and longitude limits to plot on map slices
        %Example: [-36 -35 -70 -71]; %First two coordinates are latitude
        %'default' plots default non-padding cells OR stations area (when
        %plotting data only)
u.lim_2D = [1 10 0 3]; % plotting limits for 2D in km; [ymin ymax zmin zmax];
u.ve = 1; %Vertical exaggeration (ve>1 = more exaggeration)
u.tol = 1; %maximum distance (km) that a site can be away from a profile to be projected onto that profile.

%Diagonal sections can be plotted two ways. "Smooth" diagonal sections
%averages model cells along the profile while "line" plots the actual model
%line distances. Sometimes the "line" method makes the diagonal slice look
%irregular. Sometimes "smooth" method makes the diagonal slice too few
%cells. "Interp" interpolated the "line" method so that the profile always
%looks smooth (this makes for nice figures but also means you are not
%showing the "raw" inversion model but rather an interpolated version).
u.diagonal_section_mode = 'interp'; %'smooth' or 'line' or 'interp'

%When plotting cross-sections
u.zoff = 0; % km to move a station upward in a vertical section. this prevents stations from overlapping the model

%DATA PLOTTING INFO--------------------------------------------------------

%Axis limits
u.Tlim = [0.0005 10000]; %Period limits in seconds
u.rholim = [10^0 10^3]; %Apparent resistivity limits for off-diagonals
%u.rholimdiag = u.rholim;
u.rholimdiag = [10^-3 10^3]; %Apparent resistivity limits for diagonals
u.phalim = [0 90]; %Phase limits for data plotting
u.phalimdiag = [-180 180];
u.tiplim = [-0.5 0.5]; %Tipper limits for data plotting (e.g. pseudo sections)
u.Zlim = [10^-5 10^2]; %Impedance limits for data plotting
u.residual_lim = [-5 5]; % limits for normalized residuals, impedance and tipper
u.nskip = 2; %Number of periods to skip throughout
u.sskip = 1; %Number of stations to skip throughout; 1 is plot every station

%RMS Map Scaling
u.rmslim = [0 2]; %Root mean square error limits for misfit plotting
u.rmsscale = 50; %Scale factor for plotting rms misfit as circles on a map (see plot_misfit_rms_map.m)

%Induction Vector Scaling
u.iv_scale = 2; %Induction vector scaling (e.g. 2) for plotting IVs on map
u.iv_convention = -1; %Induction vector convention. Parkinson = 1; Weise = -1

%Data interpolation (for map view of data)
u.interp_method = 'cubic'; %Interpolation method for gridding. Cubic is usually fine.
u.dx = 0.005; %Interpolation gridding for north-south direction in latitude degrees
u.dy = 0.005; %Interpolation gridding for east-west direction in longitude degrees

%Phase Tensor and Polar Diagram Scaling
u.phase_tensor_polar_scale = 1000; %Scaling for phase tensors (usually between 1000 and 10000 depending on survey size?)
u.phase_tensor_ellipse_fill = 'beta'; % quantity to use to fill phase tensor ellipses. can be 'beta' or 'phimin'
u.phase_tensor_beta_colim = [0 5]; %Phase tensor color bar limits (beta skew degrees)
u.phase_tensor_phimin_colim = [0 60]; %Phase tensor color bar limits (phi min degrees)
u.pt_pseudo_scale = 10; %Scaling for phase tensor pseudo sections (usually 1 is ok as default; >1 makes them bigger)
u.pseudo_tol_km = 100; %If a station is less than this distance from a pseudo-section line, it will not be projected on that line
                    %Units is kilometers. For example if it is 2 km, then
                    %all stations that are >2 km from the line will not be
                    %projected onto the pseudo-section
u.rose_histogram = 30; %Maximum # of stations for rose diagram histogram
u.inset_size = 0.4; %Inset size for rose diagrams on phase tensor and/or induction vector map
u.plot_inset_pt = true; % Plot rose diagram inset on phase tensor map. true or false
u.plot_inset_iv = true; % Plot rose diagram inset on induction vector map. true or false
u.inset_loc_pt = [0.1 0.6];%'NW'; % Inset location for rose diagram on phase tensor map. Can be a string: NW, NE, SW, SE, 
                       % or an [x y] location where x and y are normalized units between 0 and 1.
u.inset_loc_iv = [0.6 0.6];%'NE'; % Inset location for rose diagram on induction vector map

% Pseudo-section profile settings
u.profile_azimuth = 0; % If you choose to plot the default profile, it will intersect the center of the station grid with this azimuth. default is 0 degrees, i.e. an E-W profile
u.profile_stations = 'default'; % must be 'default' or a text file containing a list of station indices. default is all stations selected to plot on profile.
u.rotate_data_to_azimuth = 0; % set to 1 to rotate data, or 0 to not rotate data.
                              % You may want to allow data rotation when using MTplot and exploring the data at different azimuths.
                              % You should disable data rotation when
                              % viewing inversion results so that you are always viewing data in the inversion reference frame.


u.tplot = [0.1 10 100]; %Periods to plot for plot_polar_distorted.m (must be length 3)

% K-S test option
u.significance_level = 0.05; % This is the significance level (alpha) for the K-S test. If the p-value is less than alpha, the result is statistically significant.
                             % default is 0.05

%DATA EDITING OPTIONS------------------------------------------------------

%D+ smoothing percentage (% error floor of |Z|)
u.dplus_percent = 10; %Enter "0" to use the real data error rather than applying an error floor
u.edit_mode = 'Full Tensor'; %Options: 'Full Tensor' or 'Individual'.
    %When editing data points by clicking, the default is to remove
    %all components (xx,xy,yx,yy) at the clicked site/period. However,
    %under some special circumstances, you may want to edit individual
    %components one-by-one. This has not been well-tested and should be
    %used with caution.

%----------------------END USER INPUTS-------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CHECK USER INPUTS FOR VALIDITY--------------------------------------------

%Check if geofile exists and if not, set geofile to 'none'
fid = fopen(u.geofile);
if fid== -1 && ~strcmp(u.geofile,'none') %If no geo file is found on path, then set geofile to none and it will not be plotted
    disp('Geoboundary File cannot be found on your path or was not specified and so will not be plotted')
    u.geofile = 'none';
end

if fid~=-1
    fclose(fid); %Make sure you close the file after making this check other you can have problems in other scripts when opening files
end

%Check if x,y and z limits are valid
if ~strcmp(u.xylims,'default')
    if length(u.xylims)~=4
        error('Your specified model x and y limits are not valid')
    end
   
    %Sort the x and y lims from smallest to largest
    xlims = sort([u.xylims(1) u.xylims(2)]);
    ylims = sort([u.xylims(3) u.xylims(4)]);
    
    if abs(u.xylims(1)-u.xylims(2))==0 || abs(u.xylims(3)-u.xylims(4))==0
        error('Your specified model x and y limits give a range of 0 km which is invalid')
    end
    
    u.xylims = [xlims ylims];
end

if u.zmin>u.zmax
    error('Your z limits are invalid. zmin must be less than zmax')
end

%Check if lat and long limits are valid
if ~strcmp(u.maplims,'default')
    if length(u.maplims)~=4
        error('Your specified latitude and longitude limits are not valid')
    end
   
    %Sort the x and y lims from smallest to largest
    latlims = sort([u.maplims(1) u.maplims(2)]);
    lonlims = sort([u.maplims(3) u.maplims(4)]);
    
    if abs(u.maplims(1)-u.maplims(2))==0 || abs(u.maplims(3)-u.maplims(4))==0
        error('Your specified model x and y limits give a range of 0 km which is invalid')
    end
    
    u.maplims = [latlims lonlims];
end

%%% Set colormap based on above input

if strcmp(u.cmap_name,'mtcode_default') % custom jet colormap defined in S3D for many years...

    u.cmap = [0.502 0 0;0.5678 0 0;0.6337 0 0;0.6996 0 0;0.7655 0 0;
        0.8314 0 0;0.8651 0.03294 0;0.8988 0.06588 0;0.9325 0.09882 0;0.9663 0.1318 0;
        1 0.1647 0;1 0.2322 0;1 0.2996 0;1 0.3671 0;1 0.4345 0;
        1 0.502 0;1 0.5678 0;1 0.6337 0;1 0.6996 0;1 0.7655 0;
        1 0.8314 0;0.9663 0.8651 0.03294;0.9325 0.8988 0.06588;0.8988 0.9325 0.09882;0.8651 0.9663 0.1318;
        0.8314 1 0.1647;0.7655 1 0.2322;0.6996 1 0.2996;0.6337 1 0.3671;0.5678 1 0.4345;
        0.502 1 0.502;0.4345 1 0.5678;0.3671 1 0.6337;0.2996 1 0.6996;0.2322 1 0.7655;
        0.1647 1 0.8314;0.1318 0.9663 0.8651;0.09882 0.9325 0.8988;0.06588 0.8988 0.9325;0.03294 0.8651 0.9663;
        0 0.8314 1;0 0.7655 1;0 0.6996 1;0 0.6337 1;0 0.5678 1;
        0 0.502 1;0 0.4345 1;0 0.3671 1;0 0.2996 1;0 0.2322 1;
        0 0.1647 1;0 0.1318 0.9663;0 0.09882 0.9325;0 0.06588 0.8988;0 0.03294 0.8651;
        0 0 0.8314;0 0 0.8108;0 0 0.7902;0 0 0.7696;0 0 0.749;
        0 0 0.7284;0 0 0.7078;0 0 0.6873;0 0 0.6667];

elseif strcmp(u.cmap_name,'mtcode_red_blue') % a custom version of red-to-blue
    
    inc = 0.05;

    blue = [0 0 1];
    red = [1 0 0];
    for ic = 1:1/inc
        blue = [blue; [inc*ic inc*ic 1] ]; 
        red = [red; [1 inc*ic inc*ic] ]; 
    end
    bled = [blue;flipud(red(1:end-1,:))];

    for it = 1:8
        bled = [bled(1,:) - [0 0 inc]; bled];
        bled = [bled; bled(end,:) - [inc 0 0] ];
    end

    u.cmap = flipud(bled);
    
else % user entered the name of a matlab embedded colormap
    u.cmap = u.cmap_name;

end


end


