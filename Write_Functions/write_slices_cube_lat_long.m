function write_slices_cube_lat_long(m)
%
% Function which writes out each slice of a model in (lat,long,z,rho) text
% file as well a larger file of the entire cube
%
% Usage: write_slices_cube_lat_long(m)
%
% Inputs: "m" is a standard model structure
%
%

prompt = {'Choose slice number(s) e.g.  1:25 or single number'};
titles  = 'Choose slice';
def = {['1:',num2str(m.nz)]};
slices = str2num(char(inputdlg(prompt,titles,1,def))); 


disp('Writing all slices as xyzd file (columns of longitude, latitude, elevation and resistivity values)');

if ~exist('text_files','dir')
    mkdir text_files
end

m.A(isnan(m.A)) = 10^9;

for iz=slices
    %z = depth(iz)+ 0.5*(depth(iz+1)-depth(iz));
% output resistivity to a text file (use in Surfer or ArcGIS)
    disp(['Layer #',num2str(iz)]);
    str=['text_files\layer_',num2str(iz,'%02.0f'),'_depth_',num2str(0.001*m.cz(iz),'%6.2f\n'),'_km.dat'];
    fid=fopen(str,'w+');
    for ix=m.npad+1:m.nx-m.npad-1
      for iy=m.npad+1:m.ny-m.npad-1
        fprintf(fid,'%10.4f %10.4f %10.4f %10.4f \n',[m.lon(iy),m.lat(ix),m.cz(iz),m.A(ix,iy,iz)]);
      end
    end
    fclose(fid);

end   
    
% output resistivity to a text file (use in Surfer or ArcGIS) - all layers in one file
% Point in centre of cell and padding excluded
disp('Writing out one large block text file with all selected layers together')
tic; tim=0;
str=['text_files\all_layers.dat'];
fid=fopen(str,'w+');
for iz=slices
  for ix=m.npad+1:m.nx-m.npad-1
    for iy=m.npad+1:m.ny-m.npad-1 
      fprintf(fid,'%10.4f %10.4f %10.4f %10.4f \n',[m.lon(iy+5),m.lat(ix+5),m.cz(iz),m.A(ix+5,iy+5,iz)]);
      if toc>tim
          disp('...')
          tim=tim+10;
      end
    end
  end
end
toc
fclose(fid);
    
    
    
    




