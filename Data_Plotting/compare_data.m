function compare_data(d,dpred1,dpred2)
%Function to plot an observed data set and compare two separate predicted
%data sets using apparent resistivity and phase curves
%
% Usage: compare_data(d,dpred1,dpred2);
%
% Input: d = Observed data structure
%       dpred1 = Predicted Data Set #1 structure
%       dpred2 = Predicted Data Set #2 structure
%
%
% Figures saved as PNG, EPS, JPG etc (depending on user_defaults setting).
%
%test
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)
%%
G = load_geoboundary_file_list;

is = 1; check = 0;
while 1
    
    if check == 0
        set_figure_size(1);

        %Plot station locations in map view
        subplot(2,3,[3 6])
        plot_topo(d,1);
        plot_geoboundaries_geoshow(G)
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
        axis(d.lim)
        xlabel('Longitude'); ylabel('Latitude')

        plot_rho_pha(d,is); hold on %Plot observed data
        plot_rho_pha(dpred1,is); %Plot predicted data #1
        
        %Plot predicted data #2 with a different linestyle and linewidth
        subplot(2,3,1) %Diagonal rho
        h = logerrorbar(dpred2.T,dpred2.rho(:,1,is),dpred2.rhoerr(:,1,is),'--m','-r');hold on; set(h,'LineWidth',2)
        h = logerrorbar(dpred2.T,dpred2.rho(:,4,is),dpred2.rhoerr(:,4,is),'--g','-b'); set(h,'LineWidth',2)
        manual_legend('XX Predicted Data #1','-m','YY Predicted Data #1','-g','XX Predicted Data #2','--m','YY Predicted Data #2','--g');
        
        subplot(2,3,2) %Off diagonal rho
        h = logerrorbar(dpred2.T,dpred2.rho(:,2,is),dpred2.rhoerr(:,2,is),'--r','-r');hold on; set(h,'LineWidth',2)
        h = logerrorbar(dpred2.T,dpred2.rho(:,3,is),dpred2.rhoerr(:,3,is),'--b','-b'); set(h,'LineWidth',2)
        manual_legend('XY Predicted Data #1','-r','YX Predicted Data #1','-b','XY Predicted Data #2','--r','YX Predicted Data #2','--b');
        
        subplot(2,3,4); %Diagonal phase
        semilogx(dpred2.T,dpred2.pha(:,1,is),'--m','LineWidth',2); hold on
        semilogx(dpred2.T,dpred2.pha(:,4,is)+180,'--g','LineWidth',2);
        
        subplot(2,3,5); %Off diagonal phase
        semilogx(dpred2.T,dpred2.pha(:,2,is),'--r','LineWidth',2); hold on
        semilogx(dpred2.T,dpred2.pha(:,3,is)+180,'--b','LineWidth',2);

        print_figure('compare_rho_pha',[d.site{is},'_d_',d.name,'_dpred1_',dpred1.name,'_dpred2_',dpred2.name]); %Save figure
    end
    
    site_menu = menu('','Next Site','Previous Site','Select Site from List','Plot RMS Misfit Map Comparison','Return');
    
    if site_menu == 1 %Advance through stations
        is = is+1; check = 0;
        if is>d.ns
            is = ns;
        end
    elseif site_menu == 2
        is = is-1; check = 0;
        if is==0
            is = 1;
        end
    elseif site_menu == 3 %Select station from a list
    
        table(d.site);
        check = 0;

        is = listdlg('PromptString',['Select station to Plot (',num2str(d.ns),' total) '],'ListString',d.site,'Name','Select Stations','ListSize',[300 300],'SelectionMode','single');

        if isempty(is)
            break
        end
    elseif site_menu == 4 %Compute rms misfit map
        
        s1 = detailed_statistics(d,dpred1);
        s2 = detailed_statistics(d,dpred2);

        set_figure_size(2);
        subplot(1,2,1)
        plot_misfit_rms_map(d,dpred1,s1)
        title(['Misfit Between dobs and dpred1 (',dpred1.name,'). Overall RMS = ',num2str(s1.rms)],'Interpreter','none');

        subplot(1,2,2)
        plot_misfit_rms_map(d,dpred2,s2)
        title(['Misfit Between dobs and dpred2 (',dpred2.name,'). Overall RMS = ',num2str(s2.rms)],'Interpreter','none');

        print_figure('compare_rho_pha',['rms_misfit_map_d_',d.name,'_dpred1_',dpred1.name,'_dpred2_',dpred2.name]); %Save figure
        
        check = 1;

    else
        break
    end
    
end

