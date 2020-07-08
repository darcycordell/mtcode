function[decimal] = dms2deg(deg,mini,sec)

% Function to convert degree, minute, second format into 
% decimal degrees format

% Greg Nieuwenhus, 2012

decimal=deg+(mini/60)+(sec/3600);

    