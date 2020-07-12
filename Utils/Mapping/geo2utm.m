function [x,y] = geo2utm(long,lat,cent_long,orig_lat)
%
% Function to convert geographic lat-long coordinates to utm coordinates 
% 
% Usage: [x, y] = utm2geo(long,lat,cent_long,orig_lat)
%
% Dennis Rippe, 2010
%
% x: easting (in m) (don't forget about the 500,000)
% y: northing (in m)
% long: longitude (ddd.mmmm°E)
% lat: latitude (dd.mmmm°N)
% cent_long: central longitude
% orig_lat: origin latitude
%
% Documentation:
% http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.HTM
% http://www.linz.govt.nz/geodetic/conversion-coordinates/projection-conversions/transverse-mercator-preliminary-computations/index.aspx
% 
% Only tested for North American coordinates.

lat = lat*pi/180;
long = long*pi/180;
cent_long = cent_long*pi/180;
orig_lat = orig_lat*pi/180;

a = 6378137; % Equatorial radius
b = 6356752.3142; % Polar radius
f = (a-b)/a; % Flattening

k0 = 0.9996; % Scale factor
e = sqrt(1-b^2/a^2); % Eccentricity
e2 = e^2/(1-e^2);
n = (a-b)/(a+b);
rho = a*(1-e^2)./(1-e^2*(sin(lat)).^2).^(3/2); % Radius of curvature
nu = a./(1-e^2*(sin(lat)).^2).^(1/2); % Perpendicular radius of curvature
p = (long-cent_long);

A = a*(1-n+5/4*(n^2-n^3)+81/64*(n^4-n^5));
B = 3*a*n/2*(1-n+7/8*(n^2-n^3)+55/64*(n^4-n^5));
C = 15*a*n^2/16*(1-n+3/4*(n^2-n^3));
D = 35*a*n^3/48*(1-n+11/16*(n^2-n^3));
E = 315*a*n^4/512*(1-n);

S1 = A*lat-B*sin(2*lat)+C*sin(4*lat)-D*sin(6*lat)+E*sin(8*lat); % Meridional arc
S2 = A*orig_lat-B*sin(2*orig_lat)+C*sin(4*orig_lat)-D*sin(6*orig_lat)+E*sin(8*orig_lat);

K1 = (S1-S2)*k0;
K2 = k0*nu.*sin(2*lat)/4;
K3 = (k0*nu.*sin(lat).*(cos(lat)).^3/24).*(5-(tan(lat)).^2+9*e2*(cos(lat)).^2+4*e2^2*(cos(lat)).^4);
K4 = k0*nu.*cos(lat);
K5 = (k0*nu.*(cos(lat)).^3/6).*(1-(tan(lat)).^2+e2*(cos(lat)).^2);

y = K1+K2.*p.^2+K3.*p.^4;
x = K4.*p+K5.*p.^3+500000;