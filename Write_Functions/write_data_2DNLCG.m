function write_data_2DNLCG(datafile,d,mode)
% function to write a data file for use with the Rodi-Mackie NLCG 2D inversion
% Impedance units are assumed to be [V/m]/[T]
% sign convention is assumed to be -1, i.e.: exp(-i\omega t)
% orientation is assumed to be 0.0 (unrotated)
% origin is assumed to be 0,0,0

% siteName is a cell array of characters containing site names
% coord is a [ns x 3] matrix containing the site locations in x, y, and z
%       (in meters)
% data is a complex valued matrix of size [nr x nf x ns], where nr is
% either 4,8, or 12
% (diagonals, full tensor, full tensor and tipper)

if mode == 1 %TM MODE
    j = 3;
elseif mode == 2 %TE MODE
    j = 2;
else
    return
end

d.pha(:,2,:) = d.pha(:,2,:)-90;
d.pha(:,3,:) = d.pha(:,3,:)+90;

d.phaerr = d.phaerr*(pi/180);
d.rhoerr = log(d.rhoerr);

d.rho(isnan(d.rho))= 100;
d.pha(isnan(d.pha)) = -45;
d.rhoerr(isnan(d.rhoerr)) = 10^6;
d.phaerr(isnan(d.phaerr)) = 10^6;

fid1=fopen(datafile,'w+');
fprintf(fid1,'%3.0f',d.ns);
fprintf(fid1,'\n');

for k=1:d.ns
   fprintf(fid1,'%3.0f',d.nf);
   fprintf(fid1,' %5.4f',1);
   fprintf(fid1,'\n');
   for i=1:d.nf
      fprintf(fid1,'%11.5f',d.T(i));
      fprintf(fid1,'  (');
      fprintf(fid1,' %10.3f',d.rho(i,j,k));
      fprintf(fid1,'  ,');
      fprintf(fid1,' %7.2f',d.pha(i,j,k));
      fprintf(fid1,'  )');
      fprintf(fid1,'  %6.4f',d.rhoerr(i,j,k));
      fprintf(fid1,'  %6.4f',d.phaerr(i,j,k));
      fprintf(fid1,'\n');
   end
end

fclose(fid1);