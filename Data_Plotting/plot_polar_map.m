function plot_polar_map(d,polar)
%%
% Function to plot polar diagrams for all stations in map view as a
% function of period
%
% Usage: plot_polar_map(d,polar)
%
% "d" is MT data structure
% "polar" is a structure of polar diagram variables
%           See calc_polar function
%

u = user_defaults;
[L] = load_geoboundary_file_list;
[d] = set_map_projection(d);

% Convert station locations from lat/long to UTM
centlong=(max(d.loc(:,2))-min(d.loc(:,2)) )/2 + min(d.loc(:,2));
centlat=(max(d.loc(:,1))-min(d.loc(:,1)) )/2 + min(d.loc(:,1));
[loc_ew,loc_ns]=geo2utm(d.loc(:,2),d.loc(:,1),centlong,centlat); 

%========================================================================
  % Plot results for all stations on map
%========================================================================  
 
for ifreq = 1:u.nskip:d.nf  % Loop over frequencies
  
  set_figure_size(1);
  m_grid('box','fancy','tickdir','in','xlabeldir','end'); 
  hold on  
  for is = 1:u.sskip:d.ns  % Loop over stations
       
    %Normalized polar diagram radius
    r = 1.1*max([squeeze(abs(polar.x(ifreq,2,is,:))); squeeze(abs(polar.y(ifreq,2,is,:)))]);
    scale = u.phase_tensor_polar_scale/r;
    
    %Convert station locations (e.g. center of ellipse) to lat long
    [long_cent,lat_cent]=utm2geo(loc_ew(is),loc_ns(is),centlong,centlat);
    m_plot(long_cent,lat_cent,'.k')
   
    % Plot off-diagonal peanut
    [polar1,polar2]=utm2geo(squeeze(loc_ew(is)+scale*polar.y(ifreq,2,is,:)),...
        squeeze(loc_ns(is)+scale*polar.x(ifreq,2,is,:)),centlong,centlat);
    m_plot(polar1,polar2,'k-')
    
    % Plot diagonal peanut
   [polar1,polar2]=utm2geo(squeeze(loc_ew(is)+scale*polar.y(ifreq,1,is,:)),...
       squeeze(loc_ns(is)+scale*polar.x(ifreq,1,is,:)),centlong,centlat);
    m_plot(polar1,polar2,'r-')
       
  end
  
  plot_topo(d,3);
  plot_geoboundaries(L);
  
  if d.T(ifreq) > 1;
    title(['T = ',num2str(d.T(ifreq)),' s'])
  else
    title(['f = ',num2str(1/d.T(ifreq)),' Hz'])  
  end
  
  xlabel('Longitude')
  ylabel('Latitude')
  
  print_figure('polar',['polar_map_',num2str(ifreq,'%03.0f'),'_',num2str(d.T(ifreq))]); %Save figure

  
  pause(0.1)
 
end

end



