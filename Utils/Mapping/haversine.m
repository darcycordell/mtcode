function [d,theta] = haversine(lat1,lon1,lat2,lon2)
%Function which calculates the distance between two latitude and longitude
%points on the Earth's surface
%
% Usage: [d,theta] = haversine(lat1,lon1,lat2,lon2)
%
% Inputs: (lat1,lon1) is the first point in degrees
%         (lat2,lon2) is the second point in degrees
%
% Outputs:
%       d: The distance between (lat1,lon1) and (lat2,lon2) in kilometers'
%       theta: The angle between the two points in degrees

R = 6378.1; %Radius of the Earth
rad = pi/180;

dlat = lat2*rad-lat1*rad;
dlon = lon2*rad-lon1*rad;

a = sin(dlat/2).^2 + cos(lat1*rad).*cos(lat2*rad).* sin(dlon/2).^2;
c = 2*atan2(sqrt(a), sqrt(1-a));
d = R*c;

theta = atan2(sin(dlon)*cos(lat2*rad),cos(lat1*rad)*sin(lat2*rad)-sin(lat1*rad)*cos(lat2*rad)*cos(dlon))/rad;