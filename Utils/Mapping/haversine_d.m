function [lat2,lon2] = haversine_d(bearing,d,lat1,lon1)
%Function which calculates the final latitude and longitude after travelling a
%distance, d, along a bearing from a starting latitude and longitude.
%
% Usage: [lat2, lon2] = haversine_d(bearing,d,lat1,lon1)
%
%Bearing is in degrees.
%d is distance in kilometers
%Latitude and longitude are given in degrees


R = 6378.1; %Radius of the Earth
rad = pi/180;

bearing = bearing*rad;

lat1 = lat1*rad; %Current lat point converted to radians
lon1 = lon1*rad; %Current long point converted to radians

lat2 = asin(sin(lat1).*cos(d/R) + cos(lat1).*sin(d/R).*cos(bearing));
lon2 = lon1 + atan2(sin(bearing).*sin(d/R).*cos(lat1),cos(d/R)-sin(lat1).*sin(lat2));

lat2 = lat2/rad;
lon2 = lon2/rad;



