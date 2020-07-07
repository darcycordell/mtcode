function plot_induction_vector_map(d,ip,linecolor)
% Function which plots induction vectors for a single period in map view.
%
% Usage: plot_induction_vector_map(d,ip,linecolor)
%
% "d" is MT data structure
% "ip" is the index of the period to plot
% "linecolor" is a string which specifies line color (e.g. 'r','k')
%
% Note: user_defaults.m contains info about IV conventions

u = user_defaults;

if get(gcf,'number') ~= 777 % check if phase tensor map figure is plotted. if not, make new figure   
    close all;
    h = set_figure_size(888);
    if isempty(get_projection)
        d = set_map_projection(d);
    end
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    [L] = load_geoboundary_file_list;
    plot_geoboundaries(L);
    plot_topo(d,3)
elseif get(gcf,'number') == 777
    h = gcf; % need this in order to plot two insets on one figure
end

m_plot(d.loc(:,2), d.loc(:,1),'k.'); hold on; %Plot station locations

%Calculate induction vectors
%Note the negative is taken to make them point in the right convention
%(this is taken from Ersan's code in the old MTplot. Not sure why this
%is necessary).
x_vec = squeeze(-u.iv_convention*u.iv_scale*real(d.tip(ip,2,:))); % Tzy, E-W component
y_vec = squeeze(-u.iv_convention*u.iv_scale*real(d.tip(ip,1,:))); % Tzx, N-S component

%Plot vectors in map view
m_vec(1,d.loc(:,2),d.loc(:,1),x_vec,y_vec,linecolor,'shaftwidth',0.2,'headwidth',3,'headlength',4); hold on;

m_vec(1,min(d.loc(:,2)),min(d.loc(:,1)),-u.iv_convention*u.iv_scale,0,'r','shaftwidth',0.2,'headwidth',3,'headlength',4); hold on;

az = atan2d(y_vec,x_vec); % this is azimuth with geographic east as 0°, positive degrees CCW
% m_text(d.loc(:,2),d.loc(:,1),num2str(mod(az,360),'%.0f')) % debugging

if get(gcf,'number') ~= 777 % do not show title if phase tensors are already plotted    
    title(['Induction Vectors for Period: ',num2str(d.T(ip)),' s'])        
end

if u.plot_inset_iv %Plot inset if user wants

    %Plot rose diagram inset
    g = figure(100); % temporary figure to plot rose diagram
    rose_geog(az,24,u.rose_histogram,'g');

    main_fig = findobj(h,'Type','axes');
    inset_fig = findobj(g,'Type','axes');
    h_inset = copyobj(inset_fig,h);
    ax=get(main_fig,'Position');

    close(g); 
    
    set_inset_position(ax,h_inset,u.inset_loc_iv)
end

if get(gcf,'number') == 888 % if only plotting IVs, print figure 
    print_figure('iv_map',['iv_map_',num2str(ip,'%03.0f'),'_',num2str(d.T(ip)),'s']); %Save figure
    pause(0.1)
end
    


end
