function plot_rho_pha_map(d,ip)
%
% Function which plots interpolated apparent resistivity and phase for a single period in map view
%
% Usage: plot_rho_pha_map(d,ip)
%
% "d" is an MT data structure
% "ip" is the index of the period to plot.

u = user_defaults;
[L] = load_geoboundary_file_list;
[d] = set_map_projection(d);

close all

if d.ns == 1
    return
end

if isnan(d.Zerr) %If all errors are NaN, then data is predicted data (from an inversion)
    datatype = 'Modelled';
else
    datatype = 'Observed';
end

%Set up interpolation grid
xgrid = d.lim(1):u.dx:d.lim(2);     
ygrid = d.lim(3):u.dy:d.lim(4); 
[X,Y] = meshgrid(xgrid,ygrid);
  
set_figure_size(1);

rhogrid = zeros(length(ygrid),length(xgrid),4);
for i = 1:4 %Loop over components
%%
    %Create gridded apparent resistivity
    rhogrid(:,:,i) = griddata(d.loc(:,2),d.loc(:,1), squeeze(log10(d.rho(ip,i,:))),X,Y,u.interp_method);
    
    subplot(2,4,i)
    m_pcolor(X,Y,squeeze(rhogrid(:,:,i))); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);

    colormap(gca,u.cmap);
    if i == 1 || i == 4
        caxis(log10(u.rholimdiag));
    else
        caxis(log10(u.rholim));
    end
    
    add_rho_colorbar(u);
    
    set(gca,'Layer','top')

    if d.nr == 2 % fix title plotting for off-diagonal only impedances
        if i == 2 || i == 3
            title([d.responses{i-1}(2:end),' App. Res.'])
        else % clear axes for diagonal components
            cla
            colorbar('off')
        end
    else      
        title([d.responses{i}(2:end),' App. Res.'])
    end
end

phagrid = zeros(length(ygrid),length(xgrid),4);
for i = 1:4 %Loop over components
    
    %Create gridded phase
    phagrid(:,:,i) = griddata(d.loc(:,2),d.loc(:,1), squeeze(d.pha(ip,i,:)),X,Y,u.interp_method);

    if i == 3 || i == 4
        phagrid(:,:,i) = phagrid(:,:,i)+180;
    end

    subplot(2,4,i+4)
    m_pcolor(X,Y,squeeze(phagrid(:,:,i))); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);

    colormap(gca,flipud(u.cmap))
    caxis(u.phalim)
    hcb = colorbar;
    hcb.Label.String = 'Phase (degrees)';
    set(gca,'Layer','top')
    
    if d.nr == 2 % fix title plotting for off-diagonal only impedances
        if i == 2 || i == 3
            title([d.responses{i-1}(2:end),' Phase'])
        else % clear axes for diagonal components
            cla
            colorbar('off')
        end
    else           
        title([d.responses{i}(2:end),' Phase'])
    end
end


annotation('textbox', [0 0.9 1 0.08], ...
    'String', ['Interpolated ',datatype,' App. Res. and Phase @ T = ',num2str(d.T(ip)),' s'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')




end
