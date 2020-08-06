function [d, d_orig] = edit_data_tipper(d)
%
% Function which allows user to load data from EDI or ModEM data file and
% edit the tipper data
%
% Usage: [d,d_orig] = edit_data(d)
%
% Input: "d" is a standard data structure
%
% Output: "d" is the output data structure after editing.
%       "d_orig" is the original data structure used on input
%
% Note: If you edit a data point in tipper, it will NOT remove that same
% point from the impedance. HOWEVER, if you "remove a period from entire dataset" or
% "remove a site from dataset" it WILL remove both tipper and impedances at that
% frequency and/or station
%
% Options:
%   Cycle through sites one by one, click on map to view that site, enter
%   site name to jump that specific site name
%
%   Toggle D+ fit on and off and view RMS misfit. Note: set error floor for
%   D+ in user_defaults (e.g. 10%)
%
%   Edit sites by clicking individuals points, removing entire periods or
%   removing entire sites from the data structure. While editing, you can
%   hit "U" to undo last edit or "R" to restart from beginning
%
%   Save data as EDI or ModEM format.
%
%   Each click or edit updates the file data_editing_last_click.mat
%       If the program crashes for some reason, you can return to your last
%       click by loading data_editing_last_click and then running
%       edit_data(d)
%
%   Note: If you load EDIs via interpolate_edi, you CANNOT save as ModEM
%   data file. In order to save as ModEM, you need a mesh which can be
%   constructed using M3D. The workflow is to then edit EDI data, save as
%   EDIs and then load those EDIs into M3D to build mesh, model and data
%   files for inversion
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)
%%
close all

if all(isnan(d.tip(:)))
    disp('No tipper data in your data file!')
    return
end

d_orig = d;
d = set_map_projection(d);
L = load_geoboundary_file_list;
u = user_defaults;

is = 1;

set_figure_size(1);
plot_tipper(d,is);
%Plot station locations in map view
subplot(2,2,[2 4])
plot_geoboundaries_geoshow(L);
geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
axis(d.lim);
xlabel('Longitude'); ylabel('Latitude')

h = gcf;
axesObjs = get(h,'Children');
get(axesObjs,'Children');
ax(1) = axesObjs(1); %Map
ax(2) = axesObjs(2); %Imaginary Tipper
ax(3) = axesObjs(3); %Legend
ax(4) = axesObjs(4); %Real Tipper

subplot(2,2,1)
hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

subplot(2,2,3)
hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

irun = 1;

while irun

main_menu = menu('','Next Site','Previous Site','Select From Map','Enter Site Name to View','Edit Site (Clicking)','View Periods and Remove a Period from Dataset','Remove A Site from Dataset','Save Data','Quit');

if main_menu == 1 %NEXT SITE-----------------------------------------------
    is = is+1;
    if is>d.ns
        is = d.ns;
    end
    
    set_figure_size(1);
    plot_tipper(d,is);
    %Plot station locations in map view
    subplot(2,2,[2 4])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim)
    xlabel('Longitude'); ylabel('Latitude')

    h = gcf;
    axesObjs = get(h,'Children');
    get(axesObjs,'Children');
    ax(1) = axesObjs(1);
    ax(2) = axesObjs(2);
    ax(3) = axesObjs(3);
    ax(4) = axesObjs(4);

    subplot(2,2,1)
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(2,2,3)
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

elseif main_menu == 2 %PREVIOUS SITE---------------------------------------
    is = is-1;
    if is<1
        is = 1;
    end
    
    set_figure_size(1);
    plot_tipper(d,is);
    %Plot station locations in map view
    subplot(2,2,[2 4])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim);
    xlabel('Longitude'); ylabel('Latitude')

    h = gcf;
    axesObjs = get(h,'Children');
    get(axesObjs,'Children');
    ax(1) = axesObjs(1);
    ax(2) = axesObjs(2);
    ax(3) = axesObjs(3);
    ax(4) = axesObjs(4);

    subplot(2,2,1)
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(2,2,3)
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);
    
    
elseif main_menu == 3 %CLICK ON MAP TO SELECT SITE TO PLOT-----------------
    
    [X, Y] = ginput(1);%get x,y points on either phase or rho axes       

    clickedAx=gca; 
    if clickedAx == ax(2) || clickedAx == ax(3) || clickedAx == ax(4)
        subplot(2,2,[2 4]); title('Click on the map to select a site to plot'); 
    else
        subplot(2,2,[2 4]); title('');
        
        [S] = distance(Y,X,d.loc(:,1),d.loc(:,2),6371000);
        [~,is] = min(S);
        
        set_figure_size(1);
        plot_tipper(d,is);
        %Plot station locations in map view
        subplot(2,2,[2 4])
        plot_geoboundaries_geoshow(L);
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
        axis(d.lim);
        xlabel('Longitude'); ylabel('Latitude')

        h = gcf;
        axesObjs = get(h,'Children');
        get(axesObjs,'Children');
        ax(1) = axesObjs(1);
        ax(2) = axesObjs(2);
        ax(3) = axesObjs(3);
        ax(4) = axesObjs(4);

        subplot(2,2,1)
        hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
        hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

        subplot(2,2,3)
        hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
        hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);
       
    end
    
elseif main_menu == 4 %JUMP TO A SPECIFIC SITE
    
    table(d.site)
    
    pick_site = 1;
    while pick_site

        prompt={'Site Name to Go To'};
        dlg_title='Site to View';
        def={d.site{1}}; %#ok<*CCAT1>
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        i = 1;
        while 1

            if strcmp(dinp{1},d.site{i});
                disp(['Site to Show: ',d.site{i}])
                pick_site = 0;
                break
            elseif i<d.ns
                i=i+1;
            else
                disp('Error: Site Name Not Found. Check to make sure the station name you entered is in the station list. Please enter another site.')
                break
            end

        end
    
    end
    
    is = i;
    
    set_figure_size(1);
    plot_tipper(d,is);
    %Plot station locations in map view
    subplot(2,2,[2 4])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim);
    xlabel('Longitude'); ylabel('Latitude')

    h = gcf;
    axesObjs = get(h,'Children');
    get(axesObjs,'Children');
    ax(1) = axesObjs(1);
    ax(2) = axesObjs(2);
    ax(3) = axesObjs(3);
    ax(4) = axesObjs(4);

    subplot(2,2,1)
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(2,2,3)
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);


elseif main_menu == 5 %EDIT DATA-------------------------------------------
%%
    annotation('textbox', [0 0.9 1 0.08], ...
    'String', 'RIGHT CLICK WHEN FINISHED EDITING. HIT "U" to UNDO. HIT "R" TO RESTART FRESH', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

    count = 1; points=[];
    while 1 
        
        save data_editing_last_click d
        
        [X, Y,click] = ginput(1);%get x,y points on either phase or rho axes       
        if click == 3
            break
        end

        if click == 1
            clickedAx=gca; 
            if clickedAx == ax(2) || clickedAx == ax(4)
                subplot(2,2,[2 4]); title(''); 
                click_flag = 1;
            else
                subplot(2,2,[2 4]); title('Click on either of the tipper plots only!');
                click_flag = 0;
            end

            if clickedAx==ax(4) %user clicked real tipper

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(real((d.tip(:,1,is)))-((Y))).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(real((d.tip(:,2,is)))-((Y))).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click on Tx
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click on Ty
                [~,idx] = min(val);
                id = id(idx);

                if val(1)<val(2)
                    clicked_cmp = 1;
                else
                    clicked_cmp = 2;
                end
                
                points(count) = id; count = count+1;

            elseif clickedAx==ax(2) %user clicked imag tipper

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(imag((d.tip(:,1,is)))-((Y))).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(imag((d.tip(:,2,is)))-((Y))).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click on Tx
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click on Ty
                [~,idx] = min(val);
                id = id(idx);

                if val(1)<val(2)
                    clicked_cmp = 1;
                else
                    clicked_cmp = 2;
                end  
                
                points(count) = id; count = count+1;

            end


            if click_flag
                
                if strcmp(u.edit_mode,'Individual')
                    icmp = clicked_cmp;     
                else
                    icmp = 1:2;
                end
                
                d.tip(id,:,is) = NaN+1i*NaN;
                d.tiperr(id,:,is) = NaN+1i*NaN;

                set(ax(4).Children(5),'XData',d.T,'YData',real(d.tip(:,2,is)));
                set(ax(4).Children(6),'XData',d.T,'YData',real(d.tip(:,1,is)));
                set(ax(2).Children(3),'XData',d.T,'YData',imag(d.tip(:,2,is)));
                set(ax(2).Children(4),'XData',d.T,'YData',imag(d.tip(:,1,is)));

            
            end

        end

        if click==117 || click==85 %UNDO

            if ~isempty(points)

                count = count-1;
                id = points(count);
                points(count) = [];

                d.tip(id,:,is) = d_orig.tip(id,:,is);
                d.tiperr(id,:,is) = d_orig.tiperr(id,:,is);

                set(ax(4).Children(5),'XData',d.T,'YData',real(d.tip(:,2,is)));
                set(ax(4).Children(6),'XData',d.T,'YData',real(d.tip(:,1,is)));
                set(ax(2).Children(3),'XData',d.T,'YData',imag(d.tip(:,2,is)));
                set(ax(2).Children(4),'XData',d.T,'YData',imag(d.tip(:,1,is)));

            end

        end
        
    
    end
    
elseif main_menu == 6 %PLOT HISTOGRAM OF PERIODS AND DELETE PERIODS--------
    
    per_range = squeeze(sum(~isnan(d.Z),3));
    perplot = unique(per_range','rows');
    if size(perplot,1) ~= 1
        disp('Warning: Some periods have data missing at a single component but not all components. This is sketchy!')
    end
    pfig = figure(10);
    bar(log10(d.T),perplot,0.6);
    hold on
    plot(log10(d.T),ones(numel(d.T)).*numel(d.site),'k--')
    legend('No. of Stations','Total No. of Stations')
    title('Number of stations with data')
    xlabel('Log10 Period (s)')
    ylabel('Number of Stations with Data')
    
    del_per_menu = menu('','Delete Period(s) (Clicking)','Delete Period(s) (From List)','Return');
    
    if del_per_menu == 1
        title('Number of stations with data. **RIGHT CLICK TO STOP**')
        sel = []; val = 0;
        while 1
            [X,~,click] = ginput(1);
            if click == 3
                break;
            else       
                sel = [sel nearestpoint(X,log10(d.T))];
                pfig.Children(2).Children(end).XData(sel(end))=NaN;
                val = 1;
            end              
        end
        
        sel = unique(sel);
        
        close(pfig)
              
    elseif del_per_menu == 2
        
        close(pfig);
        
        table([1:d.nf]',d.T) %#ok<NBRAK>

        [sel,val] = listdlg('PromptString',['Select periods to delete (',num2str(d.nf),' total) '],'ListString',cellstr(num2str(d.T)),'Name','Delete Periods','ListSize',[300 300]);
        
    end
    
    if del_per_menu == 1 || del_per_menu == 2
        if ~val
            disp('No periods selected. Returning to Main Menu')
        else % 2 possibilities: impedance+tipper data at period, or just tipper data
            freq = sel;
            
            if all(mv(isnan(d.Z(freq,:,:)))) % if all impedance data at selected period are NaN, delete entire period from d
                disp('No impedance data detected at selected period - deleting entire period from the dataset!')
                
                d.rho(freq,:,:) = [];
                d.rhoerr(freq,:,:) = [];
                d.pha(freq,:,:) = [];
                d.phaerr(freq,:,:) = [];
                d.Z(freq,:,:) = [];
                d.Zerr(freq,:,:) = [];
                d.tip(freq,:,:) = [];
                d.tiperr(freq,:,:) = [];
                d.T(freq) = [];
                d.f(freq) = [];
                d.zrot(freq,:) = [];
                d.trot(freq,:) = [];

                d.nf = length(d.T);
            else % non-NaN impedance data at selected period - only make the tipper data NaN
                d.tip(freq,:,:) = NaN+1i*NaN;
                d.tiperr(freq,:,:) = NaN+1i*NaN;
            end
                
        set_figure_size(1);
        plot_tipper(d,is);
        %Plot station locations in map view
        subplot(2,2,[2 4])
        plot_geoboundaries_geoshow(L);
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
        axis(d.lim);
        xlabel('Longitude'); ylabel('Latitude')

        h = gcf;
        axesObjs = get(h,'Children');
        get(axesObjs,'Children');
        ax(1) = axesObjs(1);
        ax(2) = axesObjs(2);
        ax(3) = axesObjs(3);
        ax(4) = axesObjs(4);

        subplot(2,2,1)
        hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
        hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

        subplot(2,2,3)
        hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
        hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

        save data_editing_last_click d

        end
    end
    
elseif main_menu == 7 %DELETE ENTIRE SITE----------------------------------
    %%
    table(d.site)
    
    [sel,val] = listdlg('PromptString',['Select stations to delete (',num2str(d.ns),' total) '],'ListString',d.site,'Name','Delete Stations','ListSize',[300 300]);
    
    if ~val
        disp('No stations selected. Returning to Main Menu')
        
    else % 2 possibilities: impedance+tipper data at site, or just tipper data
        if all(mv(isnan(d.Z(:,:,sel)))) % if all impedance data at selected station are NaN, delete entire station from d
            disp('No impedance data detected at selected station - deleting entire station from the dataset!')
            
            d.rho(:,:,sel) = [];
            d.rhoerr(:,:,sel) = [];
            d.pha(:,:,sel) = [];
            d.phaerr(:,:,sel) = [];
            d.Z(:,:,sel) = [];
            d.Zerr(:,:,sel) = [];
            d.tip(:,:,sel) = [];
            d.tiperr(:,:,sel) = [];
            d.site(sel) = [];
            d.loc(sel,:) = [];
            d.zrot(:,sel) = [];
            d.trot(:,sel) = [];

            d.ns = length(d.rho(1,1,:));
    
            if isfield(d,'x')
                d.x(sel) = [];
                d.y(sel) = [];
                d.z(sel) = [];
            end
        else % non-NaN impedance data at selected station - only make the tipper data NaN
            d.tip(:,:,sel) = NaN+1i*NaN;
            d.tiperr(:,:,sel) = NaN+1i*NaN;
        end
       
        save data_editing_last_click d
    end
    
elseif main_menu == 8 %SAVE DATA-------------------------------------------
    
    save_menu = menu('Save As:','EDI','ModEM','d Data Structure','WSINV (does not work)');
   
    
    if save_menu == 1
        
        for is = 1:d.ns
            d.site{is} = [d.site{is},'_edited_',date];
            write_data_edi(d,is);
        end

        d.site = d_orig.site;

    elseif save_menu == 2
        
        if isfield(d,'x')
        
            prompt={'New Data Filename'};
            dlg_title='Filename';
            def={d.name,'_edited.data'};
            num_lines=1;
            dinp = inputdlg(prompt,dlg_title,num_lines,def);
            write_data_modem(dinp{1},d)
            
        else
            
            disp('ModEM Data file not saved! The data you loaded has no mesh associated with it. Save as EDIs and then use M3D to build a mesh');
        end
        
    elseif save_menu == 3
        
        save([d.name,'_edited.mat'],'d')
        
    elseif save_menu == 4
        
    end
    
else %QUIT-----------------------------------------------------------------
    irun = 0;
    
    
end

end


