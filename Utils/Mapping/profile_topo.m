function [yz] = profile_topo(strtpnt,endpnt,grdfile)
% Function which plots a topography profile in 2-D
%
% Inputs:
% strtpnt is the starting point for your profile in [lat,long]
% endpnt is the end point of your profile in [lat,long]
% grdfile is the SRTM grid file you want to use 'grdfile.grd'
%   The grdfile can be downloaded from e.g. https://maps.ngdc.noaa.gov/viewers/wcs-client/

%-------------------------EXAMPLE INPUTS----------------------------------------
% strtpnt = [-36.074, -70.062];
% endpnt = [-34.601,-76.21];
% grdfile = 'regional_chile.grd'; %Downloaded from https://maps.ngdc.noaa.gov/viewers/wcs-client/
%------------------------------------------------------------------------
%%
%load grid file using function from 
%    https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2
[Y,X,Z] = grdread2(grdfile);
Z = double(Z);
%find indices for start and end points using function from
%    https://www.mathworks.com/matlabcentral/fileexchange/8939-nearestpoint-x-y-m-
indsx = nearestpoint(strtpnt(1),X);
indsy = nearestpoint(strtpnt(2),Y);
index = nearestpoint(endpnt(1),X);
indey = nearestpoint(endpnt(2),Y);
%
%%
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%
%Step 1: Find the distance between the two points in EW and NS directions
%in lat/long
lat2 = Y(indey);
lon2 = X(index);
lat1 = Y(indsy);
lon1 = X(indsx);

ew_vec=X; ns_vec=Y;
nns=length(ns_vec); new=length(ew_vec);

dlat=lat2-lat1;
dlon=lon2-lon1;


%Step 2: Find the slope of the line and intercept of the line between the
%two points

s=(lat2-lat1)/(lon2-lon1); %slope
b=lat1-s*lon1; %intercept


%Step 3: Find the intersections between the line and the grid edges in NS
%and EW directions

%EW first (latint: dealing with known longitudes)
int=[];
for k=1:new;
    latint=s*ew_vec(k)+b; %slope formula to find intersections with known x
    if ew_vec(k)>=min([lon1 lon2]) && ew_vec(k)<=max([lon1 lon2])
        int=[int; [ew_vec(k) latint]];
    end
end

%NS second (lonint: dealing with known latitudes)
for k=1:nns;
    lonint=(ns_vec(k)-b)/s; %slope formula to find intersections with known y
    if ns_vec(k)>=min([lat1 lat2]) && ns_vec(k)<=max([lat1 lat2])
        int=[int; [lonint ns_vec(k)]];
    end
end

%add beginning and end points to intercept matrix and sort
int=[int;[lon1 lat1];[lon2 lat2]];
int=unique(int,'rows'); %remove any duplicates

%Step 4: Find the distance and midpoint of each line segment through each grid cell
%   Finds distance in terms of "lat-lon" (meaningless measure) AND distance
%   in metres

%Get lengths and midpoints of line segments
L=zeros(length(int(:,1))-1,1);L_mid=zeros(length(L),2);
for k=1:length(int)-1; %loops over all grid intersections
    L_set=sqrt((int(k+1,1)-int(k,1))^2+(int(k+1,2)-int(k,2))^2);
    L(k)=L_set;
    L_mid_set=(int(k+1,:)+int(k,:))./2;
    L_mid(k,:)=L_mid_set;
    
    %Conversion from two lat/lon coordinates to distance in metres
    %Taken from http://andrew.hedges.name/experiments/haversine/
    latint1= int(k,2); lonint1=int(k,1);
    latint2= int(k+1,2); lonint2=int(k+1,1);
    dlon=lonint2-lonint1;
    dlat=latint2-latint1;
    d2r=pi/180; %convert degrees to radians
    a = (sin(dlat*d2r/2))^2 + cos(latint1*d2r) * cos(latint2*d2r) * (sin(dlon*d2r/2))^2 ;
    c = 2 * atan2( sqrt(a), sqrt(1-a) ) ;
    d(k) = 6373000 * c; %R is radius of the Earth in metres
    
end
dist = cumsum(d);

%Step 5: Find which cell the midpoint belongs in
for j=1:length(L)
    x_grid(j) = find(ew_vec <= L_mid(j,1), 1, 'last');
    y_grid(j) = find(ns_vec <= L_mid(j,2), 1, 'last');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%Get topography profile
for i = 1:length(x_grid)
    z(i) = Z(x_grid(i),y_grid(i));
end


%
%--------------PLOTTING----------------------------------------------------
%
%Plot profile
d = haversine(strtpnt(1),strtpnt(2),endpnt(1),endpnt(2));
newx = 0:d*1000/500:d*1000;
dist(1) = 0;
vq = interp1(dist,z,newx,'linear');

newx = abs(newx-newx(end));
newx = newx./1000;

figure(3)
subplot(1,2,1)
plot(newx,vq)
axis([0 d -7000 6000])
xlabel('Distance Along Profile (km)')
ylabel('Elevation (m)')


%Plot map with profile on it
[xx yy] = meshgrid(Y,X);
indx = 1:floor(length(X)/500):length(X);
indy = 1:floor(length(Y)/500):length(Y);

subplot(1,2,2)
pcolor(xx(indx,indy),yy(indx,indy),Z(indx,indy)); shading flat
colorbar; colormap(gray); axis equal; hold on
axis([min(Y) max(Y) min(X) max(X)])
xlabel('Longitude'); ylabel('Latitude')
caxis([-5000 5000])

plot([lat1 lat2],[lon1 lon2],'o-')

[y,indsort] = sort(newx);
y = y-y(round(length(y)/2));
z = -vq(indsort);

yz = [1000*y' z'];



