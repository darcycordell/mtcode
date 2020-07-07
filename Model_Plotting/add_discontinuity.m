% This script is used for adding "exceptions" to the ModEM covariance file
%
% Exceptions split the model space into different regions. In normal ModEM
% inversions we have two regions: air and earth. In the covariance file,
% these are indicated by 0 and 1, respectively.
%
% It is possible to add additional regions to split up the model space
% above/below different discontinuities or to add known a priori geological
% information. ModEM allows for up to 8 different regions in a model (1
% through 8). The 9 is reserved for regions with fixed cells.
%
% You can set different regions to have "no smoothing" discontinuities. For
% example, if you have regions 1 and 2, you can say no smoothing between 1
% and 2 by adding an "exception" at the top of the covariance file [1 2 0]
% or [2 1 0]. 
% *****You must do this manually after saving the covariance file*****
%
% Example: 
%
% 30 40 10 %Number of cells in x, y and z directions
% 
% 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 %smoothing in the x direction for each z layer
% 
% 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3  %smoothing in the y direction for each z layer
% 
% 0.3 %smoothing in the z direction
% 
% 1 %number of times smoothing is applied
% 
% 1 %NUMBER OF EXCEPTIONS  <-----------------------------------------------
% 
% 1 2 0 %EXCEPTIONS <----------------------------------------------------
%



%% USER INPUTS
clear all
model_file = 'modem_topo_bath_91.model';
data_file = 'modem_topo_bath_91.data';
new_covariance_file = 'covariance_with_exceptions.cov';

%% LOAD DATA
m = load_model_modem(model_file);
d = load_data_modem(data_file);
[m,d] = link_model_data(m,d);

%%  LOAD SURFACE TO GET INDICES <-----------------------------------------

%This part is unique to each problem. You will have to edit this to get the
%data into the correct format

A = load('sam_slab1.0_clip.xyz'); %Here I am loading a subduction slab surface
if max(A(:,1))>180
    A(:,1) = A(:,1)-360; %Slab1.0 models sometimes default from 0 to 360
end

%Find latitude and longitudes of A matrix
lat = unique(A(:,2));
lon = unique(A(:,1));
lat = flipud(lat); %Slab1.0 model defaults decreasing latitude

%Reshape the loaded data into a 2-D matrix of depths (in km) at (lat,long)
%coordinates. This is the format you will need to have your surface in to
%continue the rest of the script
B = reshape(A(:,3),length(lon),length(lat));


%% INTERPOLATE SURFACE ONTO 2-D MODEL MESH IN LAT,LONG

[LON,LAT]=meshgrid(lon,lat);  
Bq = interp2(LON,LAT,B',m.LON,m.LAT);
%pcolor(m.LON,m.LAT,Bq) %Plot interpolated surface (optional)

%% FIND INDICES IN DIFFERENT REGIONS
%The way the script is currently written splits the model into two regions
%(above and below the interface)

ind_surf = nan(size(Bq));
%Loop through (lat,long) coordinates to find the model cells that the 2-D surface cuts through
for i = 1:m.nx
    for j = 1:m.ny
        
        ind_surf(i,j) = nearestpoint(-Bq(i,j)*1000,m.cz); %Find z model cell that surface cuts through at (i,j)

        if ~isnan(ind_surf(i,j))
            m.A(i,j,ind_surf(i,j):end) = 10^9; %Replace cells with some identifier (10^9). (For de-bugging purpose)
        end
    end
end

%Optional plotting to show regions
% plot_slice_menus(m,70,d);
% plot_diagonal_section(m,d);

%% ADD DIFFERENT REGIONS TO COVARIANCE FILE
ind_below = find(m.A==10^9); %Indices below interface
ind_air = find(isnan(m.A)); %Air indices
ind_ocean = find(m.A>0.2999 & m.A<0.3001); %If ocean exists, this finds ocean indices and assumes ocean is 0.3 Ohm m

cov = ones(m.nx,m.ny,m.nz); %Default covariance mask is 1 (for all indices above discontinuity)
cov(ind_air) = 0; %Air indices are identified by 0 in ModEM covariance file
cov(ind_ocean) = 9; %Fixed cells are identified by 9 (ocean is assumed to be fixed)
cov(ind_below) = 2; %Region #2 (below discontinuity)

a = m;
a.A = cov;

%plot_slice_menus(a,20) %View covariance file (colorbar is weird)

%% SAVE COVARIANCE FILE AND INDICES
[~,ind_above] = setdiff(1:length(m.A(:)),[ind_below; ind_air; ind_ocean]);
[ixr,iyr,izr] = ind2sub(size(a.A),ind_above);
indices = [ixr iyr izr];
save('indices_above.mat','indices');

[ixr,iyr,izr] = ind2sub(size(a.A),ind_below);
indices = [ixr iyr izr];
save('indices_below.mat','indices');

write_modem_covariance(new_covariance_file,cov,0.3*ones(m.nz,1),0.3*ones(m.nz,1),0.3,1);


