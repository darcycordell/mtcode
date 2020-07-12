function [long,lat] = utm2geo(x,y,cent_long,orig_lat)
%
% Function to convert utm coordinates to geographic lat-long coordinates
% 
% Usage: [long, lat] = utm2geo(x,y,cent_long,orig_lat)
%
% x: easting (in m) (don't forget the extra 500,000)
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


x = 500000-x;
cent_long = cent_long*pi/180;
orig_lat = orig_lat*pi/180;

a = 6378137; % Equatorial radius
b = 6356752.3142; % Polar radius
f = (a-b)/a; % Flattening

k0 = 0.9996; % Scale factor
e = sqrt(1-b^2/a^2); % Eccentricity
e2 = e^2/(1-e^2);
n = (a-b)/(a+b);

A = a*(1-n+5/4*(n^2-n^3)+81/64*(n^4-n^5));
B = 3*a*n/2*(1-n+7/8*(n^2-n^3)+55/64*(n^4-n^5));
C = 15*a*n^2/16*(1-n+3/4*(n^2-n^3));
D = 35*a*n^3/48*(1-n+11/16*(n^2-n^3));
E = 315*a*n^4/512*(1-n);

S2 = A*orig_lat-B*sin(2*orig_lat)+C*sin(4*orig_lat)-D*sin(6*orig_lat)+E*sin(8*orig_lat);

M = S2+y/k0; % Meridional arc

mu = M/(a*(1-e^2/4-3*e^4/64-5*e^6/256));
e1 = (1-(1-e^2)^(1/2))/(1+(1-e^2)^(1/2));

J1 = 3*e1/2-27*e1^3/32;
J2 = 21*e1^2/16-55*e1^4/32;
J3 = 151*e1^3/96;
J4 = 1097*e1^4/512;

fp = mu+J1*sin(2*mu)+J2*sin(4*mu)+J3*sin(6*mu)+J4*sin(8*mu); % Footprint latitude

C1 = e2*(cos(fp)).^2;
T1 = (tan(fp)).^2;
R1 = a*(1-e^2)./(1-e^2*(sin(fp)).^2).^(3/2); % Radius of curvature
N1 = a./(1-e^2*(sin(fp)).^2).^(1/2); % Perpendicular radius of curvature
D = x./(N1*k0);

Q1 = N1.*tan(fp)./R1;
Q2 = D.^2/2;
Q3 = (5+3*T1+10*C1-4*C1.^2-9*e2).*D.^4/24;
Q4 = (61+90*T1+298*C1+45*T1.^2-3*C1.^2-252*e2).*D.^6/720;
Q5 = D;
Q6 = (1+2*T1+C1).*D.^3/6;
Q7 = (5-2*C1+28*T1-3*C1.^2+8*e2+24*T1.^2).*D.^5/120;

lat = fp-Q1.*(Q2-Q3+Q4);
%long = cent_long+(Q5-Q6+Q7)/cos(fp); % Need to find out when + or - sign is needed.
long = cent_long-(Q5-Q6+Q7)./cos(fp);

lat = lat*180/pi;
long = long*180/pi;