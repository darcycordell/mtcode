function export_nlcg_mdl(model,z,y,indsit)

display ('IN : export_mdl')

% Output in format linked to geographical units %  MJU 2016-04-09
%  y - centres of model columns in 2-D model

ny = length(y);  nz = length(z);    ns = length(indsit);

% This paramter file contains profile specific information
% including station locations, azimuth of exported model, and
%along strike offset in case stations do not coincide with model
% i.e. when strike of profile is not the same geoelectric strike from GB/PT 

% % % if ~exist('cd I:\mju\field-projects\ABC\data-analysis\nlcg6\rippe-final\ABC-N\alpha3tau10SS\TETMHZ','file'); 
% % %    display ('export_nlcg_mdl_param.m missing')
% % %    display ('Returning to main menu')
% % %    return
% % % else 
% % %     export_nlcg_mdl_param
% % % end    
    
export_nlcg_mdl_param

[~,name,~] =fileparts(filename_stn);

% As a check, load station locations from text file for plotting
data = load (filename_stn,'-ascii') ;
stn_lat = data(:,1);    stn_long = data(:,2); 
nstat = length(stn_lat);

% Tie local to global co-ordinates at first station
lat1 = stn_lat(1);   long1 = stn_long(1);
% Convert this point to UTM
[EW1,NS1] = geo2utm(long1,lat1,long1,lat1);


figure(98)

% (1)  Calculate locations of model columns in UTM, relative to first station.
%      At this point the co-ordinate system has profile extending NS
  for iy=1:ny
    NS_orig(iy) = 1000*y(iy)-1000*y(indsit(1))+NS1;    EW_orig(iy) = EW1;
  end

  plot(EW_orig,NS_orig,'ko')
  hold on
  
% (2) Define point about which rotation will occur (in UTM)
  rot_cent_EW = EW_orig(indsit(1));     rot_cent_NS = NS_orig(indsit(1));

% (3) Convert to relative co-ordinate system for rotation
  EW_rel = EW_orig-rot_cent_EW;         NS_rel = NS_orig-rot_cent_NS;
  
  plot(EW_rel,NS_rel,'k+')
 
% (4) Rotate clockwise (positive)
  c = cosd(profile_azimuth);      s = sind(profile_azimuth);     R = [ c, -s ; s, c];
  for iy=1:ny
    loc = [NS_rel(iy),EW_rel(iy)];         loc = R*loc';
    NS_rot(iy) = loc(1);                  EW_rot(iy)=loc(2);
  end
  
  plot(EW_rot,NS_rot,'r.')
  
  

% (5) Return to absolute UTMs
  NS = NS_rot + rot_cent_NS;   
  EW = EW_rot + rot_cent_EW;

  plot(EW,NS,'g.-')
  
% (6) Move along strike if needed (some issues here)
   NS = NS+move_dist*cosd(move_dir); 
   EW = EW+move_dist*sind(move_dir); 

   plot(EW,NS,'b.-')
   
   print(98,'-djpeg',[name,'_UTM.jpg'])
   
   close(98)
   
% (7) Convert back to geographic
  for iy=1:ny  
    [long(iy),lat(iy)] = utm2geo(EW(iy),NS(iy),long1, lat1);
  end
    
% Plot figure to check if along strike offset is needed
figure(99)

for iy = indsit(1)-1: indsit(ns)+1;
    plot(long(iy),lat(iy),'b.'); hold on;
end
plot(long(indsit(1)),lat(indsit(1)),'ro')
%plot(long(indsit(nstat)),lat(indsit(nstat)),'ro')
plot(stn_long,stn_lat,'r*')
title(strrep(name,'_','\_'))
axis('equal')
print(99,'-djpeg',[name,'.jpg'])
close(99)

% Replace NaN values with upper most non-NaN in column
for iy =1 :ny
    for iz = nz:-1:1
       if isnan(model(iz,iy)) == 1;
           model(iz,iy) = model(ilast,iy);
       else
           ilast = iz;
       end
    end
end

% write resistivity model to file                
str=[name,'_nlcg_model.dat'];   
fid1=fopen(str,'w+');                                                                                       
                                                           
                                                           
for iy=indsit(1)-1: indsit(ns)+1;
  for iz = 1:nz-nzskip
    depth = z(iz)+0.5*(z(iz+1)-z(iz));   % Mid-point in row
    fprintf(fid1,'%9.3f %9.3f %9.3f %10.3f\n',[long(iy);lat(iy);depth; 10^(model(iz,iy))]);
  end                                                
end
fclose(fid1);                                      


% write resistivity model to file beneath stations           
str=[name,'_nlcg_model_stns.dat'];   
fid2=fopen(str,'w+');                                                                                       
                                                           
                                                           
for iy=1 : nstat
  for iz = 1:nz-nzskip
    depth = z(iz)+0.5*(z(iz+1)-z(iz));   % Mid-point in row
    fprintf(fid2,'%9.3f %9.3f %9.3f %10.3f\n',[long(indsit(iy));lat(indsit(iy));depth; 10^(model(iz,indsit(iy)))]);
  end                                                
end
fclose(fid2);         



 
end