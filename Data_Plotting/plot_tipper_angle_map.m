function plot_tipper_angle_map(d)
%
% Function which plots interpolated real and imaginary tipper angle in
% degrees E of N in map view using m_map
% Also plots interpolated real and imaginary tipper magnitude in map view
%
% Usage: plot_tipper_angle_map(d)
%
% "d" is an MT data structure
%

u = user_defaults;
[L] = load_geoboundary_file_list;
d = set_map_projection(d);

close all


%Set up interpolation grid
xgrid = d.lim(1):u.dx:d.lim(2);     
ygrid = d.lim(3):u.dy:d.lim(4); 
[X,Y] = meshgrid(xgrid,ygrid);

%Compute the real and imaginary tipper angle
angle_r=squeeze(atan2d(u.iv_convention*real(d.tip(:,1,:)),u.iv_convention*real(d.tip(:,2,:))));% in degrees, the -1 converts from Weise to Parkinson convention
angle_i=squeeze(atan2d(u.iv_convention*imag(d.tip(:,1,:)),u.iv_convention*imag(d.tip(:,2,:))));% in degrees, the -1 converts from Weise to Parkinson convention


angle_r(angle_r < 0) = 360 + angle_r(angle_r<0);
angle_i(angle_i < 0) = 360 + angle_i(angle_i<0);

mag_r=squeeze(sqrt(real(d.tip(:,1,:)).^2+real(d.tip(:,2,:)).^2));
mag_i=squeeze(sqrt(imag(d.tip(:,1,:)).^2+imag(d.tip(:,2,:)).^2));

%%  
for ip = 1:d.nf
  
    set_figure_size(1);
    %Create gridded real tipper angle
    grid_r = scatteredInterpolant(d.loc(:,2),d.loc(:,1), angle_r(ip,:)','linear','none');
    V_r = grid_r(X,Y);
    
    %Plot real angle
    subplot(2,2,1)
    m_pcolor(X,Y,V_r); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);
    flipud(colormap('jet')); caxis([0 360]); hcb = colorbar;
    hcb.Label.String = 'Angle (degrees)';
    title('Real Induction Vector Angle')
    
    %Create gridded imaginary tipper angle
    grid_i = scatteredInterpolant(d.loc(:,2),d.loc(:,1), angle_i(ip,:)','linear','none');
    V_i = grid_i(X,Y);
    
    %Plot imaginary angle
    subplot(2,2,2)
    m_pcolor(X,Y,V_i); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);
    flipud(colormap('jet')); caxis([0 360]); hcb = colorbar;
    hcb.Label.String = 'Angle (degrees)';
    title('Imag Induction Vector Angle')
    
    %Create gridded real tipper magnitude
    grid_magr = scatteredInterpolant(d.loc(:,2),d.loc(:,1), mag_r(ip,:)','linear','none');
    V_magr = grid_magr(X,Y);
    
    %Plot real magnitude
    subplot(2,2,3)
    m_pcolor(X,Y,V_magr); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);
    flipud(colormap('jet')); caxis([-1 1]); hcb = colorbar;
    hcb.Label.String = 'Magnitude';
    title('Real Induction Vector Magnitude')
    
    
    %Create gridded imaginary tipper magnitude
    grid_magi = scatteredInterpolant(d.loc(:,2),d.loc(:,1), mag_i(ip,:)','linear','none');
    V_magi = grid_magi(X,Y);
    
    %Plot imaginary tipper magnitude
    subplot(2,2,4)
    m_pcolor(X,Y,V_magi); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);
    flipud(colormap('jet')); caxis([-1 1]); hcb = colorbar;
    hcb.Label.String = 'Magnitude';
    title('Imag Induction Vector Magnitude')


    annotation('textbox', [0 0.9 1 0.08], ...
        'String', ['Interpolated Induction Vector Angle & Magnitude @ T = ',num2str(d.T(ip)),' s'], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

    pause(0.1)
end




end
