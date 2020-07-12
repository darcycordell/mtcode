function plot_misfit_residual_map(dobs,dpred,ip,s)
% Function which plots the residuals between observed and
% predicted MT data as interpolated residuals of apparent
% resistivity and phase. 
%
% Usage: plot_misfit_residual_map(dobs,dpred,ip,s)
%
% "dobs" is the observed MT data structure
% "dpred" is the predicted MT data structurce
% "ip" is the index of the period to plot
% "s" is an OPTIONAL input of a structure of detailed statistics
%

u = user_defaults;
[L] = load_geoboundary_file_list;

if ~exist('s','var')
    s = detailed_statistics(dobs,dpred);
end

close all
drho = (dobs.rho - dpred.rho)./dobs.rhoerr; %Calcuate rho residual

%Calculate phase residual. Must be done separately for xx/xy and yy/yx
%components because yy/yx components are in the -180 to -90 quadrant
dpha = zeros(dobs.nf,4,dobs.ns);
dpha = (dobs.pha - dpred.pha)./dobs.phaerr;


%Set up interpolation grid
xgrid = dobs.lim(1):u.dx:dobs.lim(2);     
ygrid = dobs.lim(3):u.dy:dobs.lim(4); 
[X,Y] = meshgrid(xgrid,ygrid);
  
set_figure_size(1);

drhogrid = zeros(length(ygrid),length(xgrid),4);
for i = 1:4 %Loop over components

    %Create gridded apparent resistivity % differences
    drhogrid(:,:,i) = griddata(dobs.loc(:,2),dobs.loc(:,1), squeeze(drho(ip,i,:)),X,Y,u.interp_method);

    subplot(2,4,i)
    m_pcolor(X,Y,squeeze(drhogrid(:,:,i))); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(dobs.loc(:,2),dobs.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);

    colormap('gray'); caxis([-2*nanstd([drho(:); dpha(:)]) 2*nanstd([drho(:); dpha(:)])]); hcb = colorbar;
    hcb.Label.String = 'Normalized Residuals';

    if dobs.nr == 2 % fix title plotting for off-diagonal only impedances
        if i == 2 || i == 3
            title([dobs.responses{i-1}(2:end),' App. Res. R.M.S. = ',num2str(s.rms_fr(i,ip),'%.2f')])
        else % clear axes for diagonal components
            cla
            colorbar('off')
        end
    else      
        title([dobs.responses{i}(2:end),' App. Res. R.M.S. = ',num2str(s.rms_fr(i,ip),'%.2f')])
    end            
    
end

dphagrid = zeros(length(ygrid),length(xgrid),4);
for i = 1:4 %Loop over components
    
    %Create gridded apparent resistivity % differences
    dphagrid(:,:,i) = griddata(dobs.loc(:,2),dobs.loc(:,1), squeeze(dpha(ip,i,:)),X,Y,u.interp_method);

    subplot(2,4,i+4)
    m_pcolor(X,Y,squeeze(dphagrid(:,:,i))); shading flat; hold on
    m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
    m_plot(dobs.loc(:,2),dobs.loc(:,1),'k.','markersize',12);
    plot_geoboundaries(L);

    colormap('gray'); caxis([-2*nanstd([drho(:); dpha(:)]) 2*nanstd([drho(:); dpha(:)])]); hcb = colorbar;
    hcb.Label.String = 'Normalized Residuals';
    
    if dobs.nr == 2 % fix title plotting for off-diagonal only impedances
        if i == 2 || i == 3
            title([dobs.responses{i-1}(2:end),' Phase'])
        else % clear axes for diagonal components
            cla
            colorbar('off')
        end
    else           
        title([dobs.responses{i}(2:end),' Phase'])
    end
           
end

annotation('textbox', [0 0.9 1 0.08], ...
    'String', ['Interpolated Normalized Residuals @ T = ',num2str(dobs.T(ip)),' s'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')


end
