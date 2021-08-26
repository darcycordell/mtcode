function plot_phase_tensor_map(d,iv)
%%
% Function which plots phase tensors in map view and colors ellipses with
% beta skew angle.
%
% Usage: plot_phase_tensor_map(d,iv)
%
% "d" is standard MT structure
% "iv" is an optional flag to plot induction vectors on the map as well. 
%       1 = plot ivs
%       0 = do not plot ivs (default)
%
% if IVs are plotted, rose diagram contains the phase tensor azimuths

[L] = load_geoboundary_file_list;
close all;

if ~exist('iv','var')
    iv = 0;
end

u = user_defaults;
[d] = set_map_projection(d);

% Convert station locations from lat/long to UTM
centlong=(max(d.loc(:,2))-min(d.loc(:,2)) )/2 + min(d.loc(:,2));
centlat=(max(d.loc(:,1))-min(d.loc(:,1)) )/2 + min(d.loc(:,1));
[loc_ew,loc_ns]=geo2utm(d.loc(:,2),d.loc(:,1),centlong,centlat); 

strike = nan(d.ns,1);
%========================================================================
  % Plot results for all stations on map
%========================================================================  
 
for ifreq = 1:u.nskip:d.nf  % Loop over frequencies
   %%
    h = set_figure_size(777); % this number is important! plot_induction_vector_map needs to know whether or not to make new figure
    m_grid('box','fancy','tickdir','in','xlabeldir','end');
    
    
    hold on  
    
    for is = 1:u.sskip:d.ns  % Loop over stations

        [p] = calc_phase_tensor(d.Z(ifreq,:,is)); %Calculate phase tensors and put in "p" structure
        strike(is) = p.strike;
        
        %Normalized radius of the ellipse
        r = 1.1*max([squeeze(abs(p.x)); squeeze(abs(p.y))]);
        scale = u.phase_tensor_polar_scale/r;
        
        %Convert ellipse locations to lat-long
        [pt1,pt2]=utm2geo(squeeze(loc_ew(is)+scale*p.y),...
            squeeze(loc_ns(is)+scale*p.x),centlong,centlat);
        
        if strcmp(u.phase_tensor_ellipse_fill,'phimin')
            m_fill(pt1,pt2,abs(p.phimin)); hold on %Plot ellipses and fill them with m_map
        elseif strcmp(u.phase_tensor_ellipse_fill,'beta')
            m_fill(pt1,pt2,abs(p.beta)); hold on
        else
            disp('Unrecognized input for u.phase_tensor_ellipse_fill. Ellipses are filled with beta skew angle values. Check your user_defaults')
            m_fill(pt1,pt2,abs(p.beta)); hold on
        end
                        
        %m_text(d.loc(is,2),d.loc(is,1)+0.0016,num2str(strike(ifreq,is),2))    % debugging 
    end
    
    m_plot(d.loc(:,2), d.loc(:,1),'k.'); hold on; %Plot station locations
    
    plot_geoboundaries(L);
    
    colormap(flipud(u.cmap));
    if strcmp(u.phase_tensor_ellipse_fill,'phimin')
        caxis(gca,u.phase_tensor_phimin_colim);
        hcb = colorbar;
        hcb.Label.String = '\Phi_{min} (degrees)';
    else
        caxis(gca,u.phase_tensor_beta_colim);
        hcb = colorbar;
        hcb.Label.String = 'Beta Skew Angle (degrees)';
    end  
    
    set(gca,'SortMethod','childorder') % make sure phase tensors and stations plot ON TOP of the geoboundaries
      
    if d.T(ifreq) > 1
        title(['T = ',num2str(d.T(ifreq)),' s'])
    else
        title(['f = ',num2str(1/d.T(ifreq)),' Hz'])  
    end

    xlabel('Longitude')
    ylabel('Latitude')
    
    if iv==1 %Plot induction vectors on top of phase tensors
        plot_induction_vector_map(d,ifreq,'k')
    end
         
    if u.plot_inset_pt %Plot inset if user wants
        
        %Plot rose diagram inset
        g = figure(100); % temporary figure to plot rose diagram
        strike_vec = [strike; strike+90; strike+180; strike+270];
        strike_vec = mod(90 - strike_vec, 360); % Reverse direction
        rose_geog(strike_vec,24,u.rose_histogram,'b');
      
        main_fig = findobj(h,'Type','axes');
        inset_fig = findobj(g,'Type','axes');
        h_inset = copyobj(inset_fig,h);
        ax=get(main_fig,'Position');
        
        if iscell(ax) % if IV inset already plotted, there will be more than one axis on the main_fig. just need to take the position of the main figure axes
            ax = ax{end};
        end
        
        close(g); 
                    
        set_inset_position(ax,h_inset,u.inset_loc_pt)              
    end
       
    if iv
        print_figure(['phase_tensor_',d.niter],['phase_tensor_iv_map_',u.phase_tensor_ellipse_fill,'_',num2str(ifreq,'%03.0f'),'_',num2str(d.T(ifreq)),'_s']); %Save figure
    else
        print_figure(['phase_tensor_',d.niter],['phase_tensor_map_',u.phase_tensor_ellipse_fill,'_',num2str(ifreq,'%03.0f'),'_',num2str(d.T(ifreq)),'_s']); %Save figure
    end

    pause(0.1)
 
end



end



