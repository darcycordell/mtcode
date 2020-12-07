function write_section(x,y,z,rho,d)
%
% Function which writes out a 2-D resistivity cross section as a text file
% which can be easily shared or used in Surfer or ArcGIS. Outputs files in
% (lat long z rho)
%
% Usage: write_section(x,y,z,rho,d)
%
% Inputs: x,y,z are the slice coordinates in UTM
%       rho is the resistivity.
%       d is a standard MT data structure (necessary for lat-long
%       conversions)
%
% The x,y,z,rho coordinates are automatically output from most
% "plot_cross_section" functions. For example, see "plot_cross_section.m",
% or "plot_cross_section_beneath_points.m".

disp('Writing text file with profile longitude, latitude, depth and rho')

if ~exist('text_files','dir')
    mkdir text_files
end

[lon,lat] = utm2geo(y*1000+500000,x*1000,d.origin(2),d.origin(1));

[lonmin, ind] = min(lon);
latmin = lat(ind);

rho(isnan(rho)) = 10^9;

tic;
str=['text_files\section_',num2str(lonmin),'_EW_',num2str(latmin),'_NS.dat'];
fid=fopen(str,'w+');
for iz=1:length(z)
  for i=1:length(x)
      fprintf(fid,'%10.4f %10.4f %10.4f %10.4f \n',[lon(i),lat(i),z(iz),rho(iz,i)]);
  end
end
toc
fclose(fid);



end
