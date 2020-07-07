function plot_misfit_rms_by_period_map(dobs,dpred,s)
%
% Function which plots rms misfit for two MT datasets in map view for each
% period. The rms for each site as the given period is show as a colored
% circle.
%
% Usage: 
% plot_misfit_rms_by_period_map(dobs,dpred)
% plot_misfit_rms_by_period_map(dobs,dpred,s)
%
% "dobs" is observed MT data. Note: the errors used to normalize the
% residuals are taken from dobs.
% "dpred" is the predicted MT data. The errors in dpred are ignored.
% "s" is an OPTIONAL input that is the stats structure output from
% detailed_statistics or detailed_staistics_2D. If s is not an input, the
% misfit stats are computed from impedances by default.
%
%%
u = user_defaults;
[L] = load_geoboundary_file_list;

%Get detailed statistics with rms by frequency and rms by site by frequency
% [~,~,rms_freq,~,rms_sf,~,~] = detailed_statistics(dobs,dpred);
if ~exist('s','var')
    disp('No stats structure input. Calculating rms misfit from impedances')
    s = detailed_statistics(dobs,dpred);
end

%Choose period range
per_range={[num2str(1),':',num2str(1),':',num2str(dobs.nf)]};
answer=char(inputdlg({'Enter period range (ex. : 1:1:13)'},'Plot Periods',1,per_range));
%The periods to plot (ptp = periods to plot)
ptp=str2num(answer);  %#ok<*ST2NM>

%Check if user inputs are valid
if isempty(ptp)
    ptp = 1;
    disp('Invalid Entry: Plotting Period 001 as Default');
end

if max(ptp) > dobs.nf || min(ptp) < 1
    ptp = 1;
    disp('Invalid Entry: Plotting Period 1 as Default')
end


misfit_by_period_menu = menu('','Interpolate Misfit on Map','Plot as Colored Circles at Site Locations');

for j = ptp %Loop over periods
    
    if misfit_by_period_menu == 1 %Interpolate misfit by period in map view
        
        %Set up interpolation grid
        xgrid = dobs.lim(1):u.dx:dobs.lim(2);     
        ygrid = dobs.lim(3):u.dy:dobs.lim(4); 
        [X,Y] = meshgrid(xgrid,ygrid);

        rmsgrid = griddata(dobs.loc(:,2),dobs.loc(:,1), s.rms_sf(:,j),X,Y,u.interp_method);
        
        m_pcolor(X,Y,rmsgrid); shading flat; hold on
        m_grid('box','fancy','xlabeldir','end','tickdir','in'); hold on
        m_plot(dobs.loc(:,2),dobs.loc(:,1),'k.','markersize',12);
        plot_geoboundaries(L);

        caxis(u.rmslim);
        colormap(jet)
        hcb = colorbar;
        hcb.Label.String = 'Root Mean Square Misfit';
        m_grid('box','fancy','tickdir','in','xlabeldir','end'); hold on

        title(['Root Mean Square By Station @ T = ',num2str(dobs.T(j)),' with RMS = ',num2str(s.rms_freq(j))]);

        
    else %Plot misfit at each site in map view as colored circles

        for i = 1:dobs.ns %Loop over sites and create a circle (xp,yp) for each site and plot it at the site location

            th = 0:pi/20:2*pi;
            r = abs(dobs.lim(4)-dobs.lim(3))/u.rmsscale; %This is the radius of the circle. Usually u.rmsscale = 50 is good
            xp = r * cos(th) + dobs.loc(i,2);
            yp = r * sin(th)*0.85+ dobs.loc(i,1);

            m_fill(xp,yp,s.rms_sf(i,j)); hold on; %
            shading flat;

        end

        caxis(u.rmslim);
        colormap(jet)
        hcb = colorbar;
        hcb.Label.String = 'Root Mean Square Misfit';
        m_grid('box','fancy','tickdir','in','xlabeldir','end'); hold on
        plot_geoboundaries(L);
        title(['Root Mean Square By Station @ T = ',num2str(dobs.T(j)),' with RMS = ',num2str(s.rms_freq(j))]);
    
        
        
    end
    
    print_figure(['misfit_stats_',dpred.niter],['misfit_map_',num2str(dobs.T(j)),'s']); %Save figure
    
    close all
    


    
end