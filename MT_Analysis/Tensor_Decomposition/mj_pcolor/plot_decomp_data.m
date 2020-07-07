function [p1] = plot_decomp_data(fig_number,long, lat, per, param, ns)

% Plots results of the tensor decomposition
% Closely based onplot_rho_pha_pseudo.m 

% Plots pseudosection of MT data in co-ordinate frame defined by profile
% direction (variable = profile_dir). 
% Jan 2016 :  Data are rotated to an azimuth = profile_dir -90

T = per(1,:);

rad = 180./pi;       nfreq = length(T);   mu = 4*pi*1e-7;      p1 = 1;
profile_dir = 90;
tlim = [0.003, 3000];  rlim = [1 1000];     plim = [0 90];
survey_name = 'Decomp test';
rot_angle = profile_dir-90;  % Angle to which data is rotated


% Process co-ords
x = long; y=lat;
x_mean = mean(x);  y_mean = mean(y); % km
x = cos(y_mean/rad)*111*(x-x_mean);y = 111*(y-y_mean); % km

for ifreq=1:nfreq
  logf(ifreq) = -log10(T(ifreq));
end


  c = cosd(rot_angle);      s = sind(rot_angle);     R = [ c, -s ; s, c];

  for is=1:ns
    loc = [x(is),y(is)];    loc = R*loc';
    x_rot(is) = loc(1);     y_rot(is)=loc(2);
  end
  
 % Add loop to put stations in order on monotonically increasing x_rot
    [x_sort,index] = sort(x_rot,'ascend');

  for is =1:ns    
    for ifreq=1:nfreq
      rho_xy(ifreq,is)= param(5,is,ifreq);  rho_yx(ifreq,is)= param(5,is,ifreq);
      pha_xy(ifreq,is)= param(7,is,ifreq);  pha_yx(ifreq,is)= param(8,is,ifreq);
    end
  end

  figure(fig_number)  
  %========================================================================
  subplot(2,2,1)
  pcolor(x_sort,logf,-log10(rho_xy));
  hold on
  shading 'flat'
  plot(x_sort,logf(1)+0.25,'kv')
  hold off;
  title(['\rho_{xy} '])
  ylabel ('freq (Hz)')
  axis([x_sort(1),x_sort(ns), -log10(tlim(2)) ,-log10(tlim(1))])
  caxis([-log10(rlim(2)),-log10(rlim(1))])

  
  %========================================================================
  subplot(2,2,3)
  pcolor(x_sort,logf,pha_xy);
  hold on;
  shading 'flat'
  plot(x_sort,logf(1)+0.25,'kv')
  hold off;
  title(['\Phi_{xy}'])
  ylabel ('freq (Hz)')
  axis([x_sort(1),x_sort(ns), -log10(tlim(2)) ,-log10(tlim(1))])
  xlabel ('distance (km)')
  caxis([plim(1),plim(2)])
  
  %========================================================================
  subplot(2,2,2)
  pcolor(x_sort,logf,-log10(rho_yx));
  hold on;
  caxis([-log10(rlim(2)),-log10(rlim(1))])
  shading 'flat'
  plot(x_sort,logf(1)+0.25,'kv')
  hold off;
  title(['\rho_{yx}'])
  ylabel ('freq (Hz)')
  axis([x_sort(1),x_sort(ns), -log10(tlim(2)) ,-log10(tlim(1))])
  
  %========================================================================
  subplot(2,2,4)
  pcolor(x_sort,logf,pha_yx);
  hold on;
  shading 'flat'
  plot(x_sort,logf(1)+0.25,'kv')
  hold off;
  title(['\Phi_{yx}'])
  xlabel ('distance (km)')
  ylabel ('freq (Hz)')
  axis([x_sort(1),x_sort(ns),  -log10(tlim(2)) ,-log10(tlim(1))])
  caxis([plim(1),plim(2)])
  xlabel ('distance (km)')
  %========================================================================

  top_title([survey_name,' : \theta =', num2str(rot_angle),'^o : Profile = ', num2str(rot_angle+90),'^o']);
  
  
  eval(['print -djpeg100 pseudo_rho_pha_',num2str(rot_angle),'_deg.jpg']);

  % Output TE apparent resistivity to file
  str=[char(survey_name),'_app_res_TE_mode.dat'];
  fid1=fopen(str,'w+');
  fprintf(fid1,'               ');
    for is=1:ns;        fprintf(fid1,'%9.3f',x_sort(is));   end
    fprintf(fid1,'\n');
    fprintf(fid1,'\n');
    for ifreq=1:nfreq
      fprintf(fid1,'%9.3f      ',T(ifreq));
    for is=1:ns; fprintf(fid1,'%9.3f',rho_xy(ifreq,is)); end
    fprintf(fid1,'\n');
  end
  fclose(fid1);
  
  
    % Output TE phase to file
  str=[char(survey_name),'_phase_TE_mode.dat'];
  fid1=fopen(str,'w+');
  fprintf(fid1,'               ');
    for is=1:ns;        fprintf(fid1,'%9.3f',x_sort(is));   end
    fprintf(fid1,'\n');
    fprintf(fid1,'\n');
    for ifreq=1:nfreq
      fprintf(fid1,'%9.3f      ',T(ifreq));
    for is=1:ns; fprintf(fid1,'%9.3f',pha_xy(ifreq,is)); end
    fprintf(fid1,'\n');
  end
  fclose(fid1);
  
    % Output TM apparent resistivity to file
  str=[char(survey_name),'_app_res_TM_mode.dat'];
  fid1=fopen(str,'w+');
  fprintf(fid1,'               ');
    for is=1:ns;        fprintf(fid1,'%9.3f',x_sort(is));   end
    fprintf(fid1,'\n');
    fprintf(fid1,'\n');
    for ifreq=1:nfreq
      fprintf(fid1,'%9.3f      ',T(ifreq));
    for is=1:ns; fprintf(fid1,'%9.3f',rho_yx(ifreq,is)); end
    fprintf(fid1,'\n');
  end
  fclose(fid1);
  
  
    % Output TM phase to file
  str=[char(survey_name),'_phase_TM_mode.dat'];
  fid1=fopen(str,'w+');
  fprintf(fid1,'               ');
    for is=1:ns;        fprintf(fid1,'%9.3f',x_sort(is));   end
    fprintf(fid1,'\n');
    fprintf(fid1,'\n');
    for ifreq=1:nfreq
      fprintf(fid1,'%9.3f      ',T(ifreq));
    for is=1:ns; fprintf(fid1,'%9.3f',pha_yx(ifreq,is)); end
    fprintf(fid1,'\n');
  end
  fclose(fid1);
  

end 