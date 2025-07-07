function [d, d_orig] = edit_data(d)
%
% Function which allows user to load data from EDI or ModEM data file and
% edit the data as well as apply D+ smoothing.
%
% Usage: [d,d_orig] = edit_data(d)
%
% Input: "d" is a standard data structure
%
% Output: "d" is the output data structure after editing.
%       "d_orig" is the original data structure used on input
%
% Note: If you edit a data point in impedance, it will NOT remove that same
% point from the tipper. HOWEVER, if you "remove a period from entire dataset" or
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
%

close all
u = user_defaults;

d_orig = d; %Hold onto original data stored in variable

d = set_map_projection(d);
L = load_geoboundary_file_list;


%---------------SET UP ERROR BAR BOUNDS FOR LATTER PLOTTING----------------
% This section is necessary for toggling the D+ error floor on/off

%Find upper and lower bounds for rho error bars in XY and YX
dYo = [d.rho(:,2:3,:)-d.rhoerr(:,2:3,:) d.rho(:,2:3,:)+d.rhoerr(:,2:3,:)];
dYo(dYo<=0) = 10^-40; %Avoid zero for log plot
dYo = dYo(:,[1 3 2 4],:); %dYo is the original error bars (before error floor is applied)
%   It is a (nf x 4 x ns) matrix with the 4 columns in the second dimension
%   set as:
%       [lower_bound on XY, upper_bound on XY, lower_bound on YX, upper_bound on YX]

dYorig = dYo; %Hold onto original, unedited error bar bounds

%Calculate error floor for entire dataset in impedance
derr = d;
for i = 1:d.ns
    for j = 2:3
        derr.Zerr(:,j,i) = apply_errorfloor(real(d.Zerr(:,j,i)),u.dplus_percent/100,abs(d.Z(:,j,i)));
    end
end
[~,~,derr.rhoerr,derr.phaerr] = calc_rho_pha(derr.Z,derr.Zerr,derr.T); %Calculate new rhoerr and phaerr
derr_orig = derr; %Hold onto original, unedited error floor bounds

%Find upper and lower error floor bounds
dY = [derr.rho(:,2:3,:)-derr.rhoerr(:,2:3,:) derr.rho(:,2:3,:)+derr.rhoerr(:,2:3,:)];
dY(dY<=0) = 10^-40;
dY = dY(:,[1 3 2 4],:); %dY has the same dimensions as dYo

dYd = dY; %Hold onto original, unedited upper and lower error floor bounds

%--------------------------------------------------------------------------

%Initialize flags: is is station counter, dplus_flag is 1 (D+ fit on) or 0 (off)
is = 1; dplus_flag = 0; don = 2; %"don" was an old variable that I don't think is necessary

set_figure_size(1);
plot_rho_pha(d,is);

%Plot station locations in map view
subplot(2,3,[3 6])
plot_geoboundaries_geoshow(L);
geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15);
axis(d.lim);
xlabel('Longitude'); ylabel('Latitude')

[ax] = get_ax; %Get axes variables of the 4 subplots
%   "ax" is a 4x1 axis handle with ax(1) = rho diagonals, ax(2) = rho
%   off-diagonals, ax(3) = pha diagonals and ax(4) = pha off-diagonals

%"hd" are the figure handles for the D+ plots to be toggled on/off
subplot(ax(2))
hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

subplot(ax(4))
hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

irun = 1; residual_figure = 0;

while irun

main_menu = menu('','Next Site','Previous Site','Select From Map','Enter Site Name to View','D+ Toggle On/Off','Set D+ Error Floor','Plot D+ Misfit','Edit Site (Clicking)','View Periods and Remove a Period from Dataset','Remove A Site from Dataset','Add Noise','Save Data','Quit');

if ~isnumeric(residual_figure) %If the D+ residuals figure is open then close it
    if ishandle(residual_figure)
        close(residual_figure);
        residual_figure = 0;
    end
end

if main_menu == 1 %NEXT SITE-----------------------------------------------
    is = is+1;
    if is>d.ns
        is = d.ns;
    end
    
    set_figure_size(1);
    plot_rho_pha(d,is);
    %Plot station locations in map view
    subplot(2,3,[3 6])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim);
    xlabel('Longitude'); ylabel('Latitude')
   
    [ax] = get_ax;   

    %Re-initialize D+ figure handles
    subplot(ax(2))
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(ax(4))
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);
    
    if dplus_flag
        %If D+ is toggled on, then calculate D+ and plot the curves
        [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
        [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

        set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
        set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
        set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
        set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
        
        %Update plotted error bars to reflect the error floor that was
        %applied. ax(4) is the pha off-diagonals
        % Because pha off-diagonals are plotted using the "errorbar"
        % built-in function, ax(4) has 4 Children plotted on the axis.
        % Children(1) is the YX D+ curve, Children(2) is the XY D+ curve,
        % Children(3) is the YX errorbar data and Children(4) is the XY
        % errorbar data.
        %
        % The errorbar data contains the YNegativeDelta and YPositiveDelta
        % variables which are equal to the error in the phase for each
        % component.
        ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

        %ax(2) contains the error bars and data for the off-diagonal rho
        %data. Because this is plotted using logerrorbar function, the
        %handles and Children are much more confusing.
        %
        %ax(2).Children is a long vector of Line handles. Children(1) is YX
        %D+ curve fit, Children(2) is XY D+ curve fit. 3 and 4 are the
        %legend labels, and the error bar lines for YX begin at Children(5)
        %through Children(d.nf+4). 
        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
            %Confusingly, the Line data are sequenced
            %in the opposite order they are plotted (first item (e.g.
            %Children(1)) is the last item that was plotted). For this
            %reason the section of Children counts forward but the sequence
            %of error bar bounds in dY counts backwards (subtracting "i"
            %counter)
        end
        
        %ax(2).Children continues with Children(d.nf+5) is the YX data.
        %Then Children(d.nf+6) through Children(2*d.nf+5) are the error bar
        %bounds
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dY(k,1:2,is);
            k = k-1; %Similar to above, the sequence is reversed for the Children versus for the dY frequencies
            %       so must count backwards
        end
        
        %Note: ax(2).Children(2*d.nf+6) is the XY data
        
        subplot(2,3,[3 6])
        title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])
        
    end

elseif main_menu == 2 %PREVIOUS SITE---------------------------------------
    is = is-1;
    if is<1
        is = 1;
    end
    
    set_figure_size(1);
    plot_rho_pha(d,is);
    %Plot station locations in map view
    subplot(2,3,[3 6])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim);
    xlabel('Longitude'); ylabel('Latitude')

    [ax] = get_ax;    

    subplot(ax(2))
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(ax(4))
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);
    
    if dplus_flag
        
        %See comments in "Next Site" menu
        [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
        [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

        set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
        set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
        set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
        set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
        
        ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
        end
        
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dY(k,1:2,is);
            k = k-1;
        end

        subplot(2,3,[3 6])
        title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])
        
    end
    
elseif main_menu == 3 %CLICK ON MAP TO SELECT SITE TO PLOT-----------------
    
    [X, Y] = ginput(1);%get x,y points on either phase or rho axes       

    clickedAx=gca; %We have to figure out which axis the user clicked on to make sure they clicked on the map
    if clickedAx==ax(1) || clickedAx == ax(2) || clickedAx == ax(3) || clickedAx == ax(4)
        subplot(2,3,[3 6]); title('Click on the map to select a site to plot'); 
    else
        subplot(2,3,[3 6]); title('');
        
        %Find the distance from the clicked point to the sites and find
        %minimum
        [S] = distance(Y,X,d.loc(:,1),d.loc(:,2),6371000);
        [~,is] = min(S);
        
        %Plot selected site
        set_figure_size(1);
        plot_rho_pha(d,is);
        %Plot station locations in map view
        subplot(2,3,[3 6])
        plot_geoboundaries_geoshow(L);
        geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
        geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
        axis(d.lim);
        xlabel('Longitude'); ylabel('Latitude')

        [ax] = get_ax;       

        subplot(ax(2))
        hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
        hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

        subplot(ax(4))
        hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
        hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

        if dplus_flag
            %See comments in "Next Menu" block of code
            [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
            [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

            set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
            set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
            set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
            set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
            
            ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
            ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
            ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
            ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

            for i = don+3:d.nf+4
                ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
            end

            k = d.nf;
            for i = don+d.nf+4:don+2*d.nf+3   
                ax(2).Children(i).YData = dY(k,1:2,is);
                k = k-1;
            end

            subplot(2,3,[3 6])
            title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])

        end
       
    end
    
elseif main_menu == 4 %JUMP TO A SPECIFIC SITE
    
    %table(d.site);

    prompt={'Site Name to Go To'};
    dlg_title='Site to View';
    def={d.site{1}}; %#ok<*CCAT1>
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    isq = find(strcmp(dinp{1},d.site));

    if isempty(isq)
        disp('Error: Site Name Not Found. Check to make sure the station name you entered is in the station list. Please enter another site.')
        isq = is;
    end
    
    
    is = isq;
    
    set_figure_size(1);
    plot_rho_pha(d,is);
    %Plot station locations in map view
    subplot(2,3,[3 6])
    plot_geoboundaries_geoshow(L);
    geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
    geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
    axis(d.lim);
    xlabel('Longitude'); ylabel('Latitude')

    [ax] = get_ax;   

    subplot(ax(2))
    hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
    hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

    subplot(ax(4))
    hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
    hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);
    
    if dplus_flag
        %See comments in "Next Menu" block of code
        [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
        [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

        set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
        set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
        set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
        set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
        
        ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
        end
        
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dY(k,1:2,is);
            k = k-1;
        end

        subplot(2,3,[3 6])
        title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])
        
    end

elseif main_menu == 5 %TOGGLE D+ ON/OFF------------------------------------
    
    if dplus_flag
        
        %If D+ is toggled off, then set D+ curves to NaN
        set(hd(1),'XData',NaN,'YData',NaN);
        set(hd(2),'XData',NaN,'YData',NaN);
        set(hd(3),'XData',NaN,'YData',NaN);
        set(hd(4),'XData',NaN,'YData',NaN);
        
        
        %Reset error bars to the original error bars with no error floor
        %applied
        ax(4).Children(1+don).YNegativeDelta = d.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = d.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = d.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = d.phaerr(:,2,is);

        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dYo(d.nf+1-(i-4),3:4,is);
        end
        
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dYo(k,1:2,is);
            k = k-1;
        end
        
        subplot(2,3,[3 6])
        title('')
        
        dplus_flag  = 0;
        
    else
    
        %If D+ is toggled on then calculate D+
        [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
        [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);
        
        %Set D+ curves figure handles to the computed D+ curves
        set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
        set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
        set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
        set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
        
        %Display plotted errorbars as calculated error floors
        % Note: Error floors are calculated for impedance but plots are
        % shown as rho/pha
        ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
        end
        
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dY(k,1:2,is);
            k = k-1;
        end
        
        subplot(2,3,[3 6])
        title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])

        dplus_flag = 1;
        
    end
    
    
elseif main_menu == 6 %CHANGE D+ ERROR FLOOR
    
    answer=char(inputdlg({'Enter D+ Error % (ex. 5)'}));
    u.dplus_percent = str2double(answer);
    
    %---------------SET UP ERROR BAR BOUNDS FOR LATTER PLOTTING----------------
    % This section is necessary for toggling the D+ error floor on/off

    %Find upper and lower bounds for rho error bars in XY and YX
    dYo = [d.rho(:,2:3,:)-d.rhoerr(:,2:3,:) d.rho(:,2:3,:)+d.rhoerr(:,2:3,:)];
    dYo(dYo<=0) = 10^-40; %Avoid zero for log plot
    dYo = dYo(:,[1 3 2 4],:); %dYo is the original error bars (before error floor is applied)
    %   It is a (nf x 4 x ns) matrix with the 4 columns in the second dimension
    %   set as:
    %       [lower_bound on XY, upper_bound on XY, lower_bound on YX, upper_bound on YX]

    dYorig = dYo; %Hold onto original, unedited error bar bounds

    %Calculate error floor for entire dataset in impedance
    derr = d;
    for i = 1:d.ns
        for j = 2:3
            derr.Zerr(:,j,i) = apply_errorfloor(real(d.Zerr(:,j,i)),u.dplus_percent/100,abs(d.Z(:,j,i)));
        end
    end
    [~,~,derr.rhoerr,derr.phaerr] = calc_rho_pha(derr.Z,derr.Zerr,derr.T); %Calculate new rhoerr and phaerr
    derr_orig = derr; %Hold onto original, unedited error floor bounds

    %Find upper and lower error floor bounds
    dY = [derr.rho(:,2:3,:)-derr.rhoerr(:,2:3,:) derr.rho(:,2:3,:)+derr.rhoerr(:,2:3,:)];
    dY(dY<=0) = 10^-40;
    dY = dY(:,[1 3 2 4],:); %dY has the same dimensions as dYo

    dYd = dY; %Hold onto original, unedited upper and lower error floor bounds

    %--------------------------------------------------------------------------
    
    if dplus_flag
        
        %If D+ is toggled on then calculate D+
        [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
        [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);
        
        %Set D+ curves figure handles to the computed D+ curves
        set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
        set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
        set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
        set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
        
        %Display plotted errorbars as calculated error floors
        % Note: Error floors are calculated for impedance but plots are
        % shown as rho/pha
        ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
        ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
        ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

        for i = don+3:d.nf+4
            ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
        end
        
        k = d.nf;
        for i = don+d.nf+4:don+2*d.nf+3   
            ax(2).Children(i).YData = dY(k,1:2,is);
            k = k-1;
        end
        
        subplot(2,3,[3 6])
        title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])

    end
    
elseif main_menu == 7 %PLOT D+ MISFIT
    

    [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
    [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

    [rhoxy_orig, phaxy_orig,rhoxyerr,phaxyerr] = calc_rho_pha(dp_xy.Zorig,dp_xy.Zerr,dp_xy.T);
    [rhoyx_orig, phayx_orig,rhoyxerr,phayxerr] = calc_rho_pha(dp_yx.Zorig,dp_yx.Zerr,dp_yx.T);

    drhoxy = (dp_xy.rho - rhoxy_orig)./rhoxyerr; drhoyx = (dp_yx.rho - rhoyx_orig)./rhoyxerr;
    dphaxy = (dp_xy.pha - phaxy_orig)./phaxyerr; dphayx = (dp_yx.pha - phayx_orig)./phayxerr;
    
    residual_figure = figure(100);
    subplot(2,1,1)
    semilogx(dp_xy.T,drhoxy,'or'); hold on;
    semilogx(dp_yx.T,drhoyx,'sb');
    xlabel('Period (s)');
    ylabel('App Rho D+ Residual')
    title(['Normalized D+ Residuals for Site ',d.site{is}])
    legend('XY','YX')
    
    subplot(2,1,2)
    semilogx(dp_xy.T,dphaxy,'or'); hold on;
    semilogx(dp_yx.T,dphayx,'sb');
    xlabel('Period (s)');
    ylabel('Phase D+ Residual')
    
elseif main_menu == 8 %EDIT DATA-------------------------------------------

    %Click to edit points, Press "U" to undo last click and press "R" to
    %reset. Note: "R" resets **ALL** edits at all sites
    annotation('textbox', [0 0.9 1 0.08], ...
    'String', 'RIGHT CLICK WHEN FINISHED EDITING. HIT "U" to UNDO.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

    count = 1; points=[]; %Points is a vector of clicked indices to keep track of "undo" commands
    
    while 1 
        
        %In case the program crashes, the d data structure is saved after
        %each click to update so you don't lose your work.
        save data_editing_last_click d
        
        origAx = gca;
            
        [X, Y,click] = ginput_editdata(1);%get x,y points on either phase or rho axes       
        
        clickedAx = gca;
        
%         if clickedAx == origAx
%             1
%         end
        
        if click == 3
           break
        end

        if click == 1 % left click (MASK DATA POINT)
            if clickedAx==ax(1) || clickedAx == ax(2) || clickedAx == ax(3) || clickedAx == ax(4)
                subplot(2,3,[3 6]); title(''); 
                click_flag = 1;
            else
                subplot(2,3,[3 6]); title('Click on either the apparent resistivity or phase plots only!');
                click_flag = 0;
            end

            if clickedAx==ax(1) %user clicked rhoxx axis

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10((d.rho(:,1,is)))-log10((Y))).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10((d.rho(:,4,is)))-log10((Y))).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click on xx
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click on yy
                [~,idx] = min(val);
                id = id(idx);
                
                if val(1)<val(2)
                    clicked_cmp = 1;
                else
                    clicked_cmp = 4;
                end

                points(count) = id; count = count+1;

            elseif clickedAx==ax(2) %user clicked rhoxy axis

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10((d.rho(:,2,is)))-log10((Y))).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10((d.rho(:,3,is)))-log10((Y))).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click on xy
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click on yx
                [~,idx] = min(val);
                id = id(idx);
                
                if val(1)<val(2)
                    clicked_cmp = 2;
                else
                    clicked_cmp = 3;
                end

                points(count) = id; count = count+1;

            elseif clickedAx==ax(3) %user clicked phaxx axis

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10(d.pha(:,1,is))-log10(Y)).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10(d.pha(:,4,is))-log10(Y)).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click
                [~,idx] = min(val);
                id = id(idx);
                
                if val(1)<val(2)
                    clicked_cmp = 1;
                else
                    clicked_cmp = 4;
                end

                points(count) = id; count = count+1;

            elseif clickedAx==ax(4) %user clicked phaxy axis

                dist1=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10(d.pha(:,2,is))-log10(Y)).^2);
                dist2=sqrt(((log10(d.T(:,1)))-log10((X))).^2+(log10(d.pha(:,3,is)+180)-log10(Y)).^2);
                [val(1), id(1)]=min(dist1);%index of closest point to mouse click
                [val(2), id(2)]=min(dist2);%index of closest point to mouse click
                [~,idx] = min(val);
                id = id(idx);
                
                if val(1)<val(2)
                    clicked_cmp = 2;
                else
                    clicked_cmp = 3;
                end

                points(count) = id; count = count+1;

            end


            if click_flag % make data point at selected period and station NaN
                
                if strcmp(u.edit_mode,'Individual')
                    icmp = clicked_cmp;     
                else
                    icmp = 1:4;
                end
                
                d.rho(id,icmp,is) = NaN;
                d.pha(id,icmp,is) = NaN;
                d.Z(id,icmp,is) = NaN+1i*NaN;
                d.rhoerr(id,icmp,is) = NaN;
                d.phaerr(id,icmp,is) = NaN;
                d.Zerr(id,icmp,is) = NaN+1i*NaN;
                
                dY(id,icmp,is) = NaN;
                dYo(id,icmp,is) = NaN;
                derr.phaerr(id,icmp,is)=NaN;
                
                if all(isnan(d.rho(:,:,is)))
                    dplus_flag = 0;
                end
                
                [ax] = set_ax_data(d, ax, is, don);
                             
                if isgraphics(ax(1)) % only true for plots with xx and yy
                    set(ax(1).Children(d.nf+3-id),'Visible','off'); % only needed for rho axes bc of errorbarloglog
                    set(ax(1).Children(2*d.nf+4-id),'Visible','off');
                end
                set(ax(2).Children(don+d.nf+3-id),'Visible','off');
                set(ax(2).Children(don+2*d.nf+4-id),'Visible','off');
                          
            end

        end

        if click==117 || click==85 % u or U key is pressed (UNDO)

            if ~isempty(points)

                %Restore the last clicked index back to its original value
                % reverse counter by one
                count = count-1;
                id = points(count); %recover last clicked point
                points(count) = []; %make points vector smaller
                
                d.rho(id,:,is) = d_orig.rho(id,:,is);
                d.pha(id,:,is) = d_orig.pha(id,:,is);
                d.Z(id,:,is) = d_orig.Z(id,:,is);
                d.rhoerr(id,:,is) = d_orig.rhoerr(id,:,is);
                d.phaerr(id,:,is) = d_orig.phaerr(id,:,is);
                d.Zerr(id,:,is) = d_orig.Zerr(id,:,is);
                
                dYo(id,:,is) = dYorig(id,:,is);
                dY(id,:,is) = dYd(id,:,is);
                
                derr.phaerr(id,:,is) = derr_orig.phaerr(id,:,is);
               
                [ax] = set_ax_data(d, ax, is, don);
                
                if isgraphics(ax(1)) % only true for plots with xx and yy
                    set(ax(1).Children(d.nf+3-id),'Visible','on');
                    set(ax(1).Children(2*d.nf+4-id),'Visible','on');
                end                
                set(ax(2).Children(d.nf+3-id+don),'Visible','on');
                set(ax(2).Children(2*d.nf+4-id+don),'Visible','on');

            end

        end
        
        if dplus_flag
            
            [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
            [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);
            
            set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
            set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
            set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
            set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);
            
            ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
            ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
            ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
            ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

            for i = don+3:d.nf+4
                ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
            end

            k = d.nf;
            for i = don+d.nf+4:don+2*d.nf+3   
                ax(2).Children(i).YData = dY(k,1:2,is);
                k = k-1;
            end

            subplot(2,3,[3 6])
            title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])



        end
    
    end
    
elseif main_menu == 9 %PLOT HISTOGRAM OF PERIODS AND DELETE PERIODS--------
    
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
        close(pfig);
              
    elseif del_per_menu == 2
        
        close(pfig);
        
        table([1:d.nf]',d.T) %#ok<NBRAK>

        [sel,val] = listdlg('PromptString',['Select periods to delete (',num2str(d.nf),' total) '],'ListString',cellstr(num2str(d.T)),'Name','Delete Periods','ListSize',[300 300]);
        
    end
    
    if del_per_menu == 1 || del_per_menu == 2
        if ~val
            disp('No periods selected. Returning to Main Menu')
        else

            freq = sel;

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

            dY(freq,:,:) = [];
            dYo(freq,:,:) = [];
            derr.phaerr(freq,:,:) = [];

            set_figure_size(1);
            plot_rho_pha(d,is);
            %Plot station locations in map view
            subplot(2,3,[3 6])
            plot_geoboundaries_geoshow(L);
            geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
            geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
            axis(d.lim);
            xlabel('Longitude'); ylabel('Latitude')

            [ax] = get_ax;  

            subplot(ax(2))
            hd(1) = loglog(NaN,NaN,'-r','LineWidth',1.5);
            hd(2) = loglog(NaN,NaN,'-b','LineWidth',1.5);

            subplot(ax(4))
            hd(3) = semilogx(NaN,NaN,'-r','LineWidth',1.5);
            hd(4) = semilogx(NaN,NaN,'-b','LineWidth',1.5);

            if dplus_flag

                [dp_xy] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T,u.dplus_percent);
                [dp_yx] = dplus(d.Z(:,3,is)*exp(1i*pi),d.Zerr(:,3,is),d.T,u.dplus_percent);

                set(hd(1),'XData',dp_xy.T,'YData',dp_xy.rho);
                set(hd(2),'XData',dp_yx.T,'YData',dp_yx.rho);
                set(hd(3),'XData',dp_xy.T,'YData',dp_xy.pha);
                set(hd(4),'XData',dp_yx.T,'YData',dp_yx.pha);

                ax(4).Children(1+don).YNegativeDelta = derr.phaerr(:,3,is);
                ax(4).Children(2+don).YNegativeDelta = derr.phaerr(:,2,is);
                ax(4).Children(1+don).YPositiveDelta = derr.phaerr(:,3,is);
                ax(4).Children(2+don).YPositiveDelta = derr.phaerr(:,2,is);

                for i = don+3:d.nf+4
                    ax(2).Children(i).YData = dY(d.nf+1-(i-4),3:4,is);
                end

                k = d.nf;
                for i = don+d.nf+4:don+2*d.nf+3   
                    ax(2).Children(i).YData = dY(k,1:2,is);
                    k = k-1;
                end


                subplot(2,3,[3 6])
                title(['XY RMS = ',num2str(dp_xy.rms),'. YX RMS = ',num2str(dp_yx.rms)])

            end

            save data_editing_last_click d
    
        end
        
    end

    
elseif main_menu == 10 %DELETE ENTIRE SITE----------------------------------
    %%
    table(d.site);
    
    [sel,val] = listdlg('PromptString',['Select stations to delete (',num2str(d.ns),' total) '],'ListString',d.site,'Name','Delete Stations','ListSize',[300 300]);
    
    if ~val
        disp('No stations selected. Returning to Main Menu')
        
    else
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
        
        dY(:,:,sel) = [];
        dYo(:,:,sel) = [];
        derr.phaerr(:,:,sel) = [];
    
        if isfield(d,'x')
            d.x(sel) = [];
            d.y(sel) = [];
            d.z(sel) = [];
        end
        
        save data_editing_last_click d
    end
    
elseif main_menu == 11 %ADD NOISE
    
    d = add_noise(d,0);
                         
elseif main_menu == 12 %SAVE DATA-------------------------------------------
    
    save_menu = menu('Save As:','EDI (All sites)','ModEM Data File','d Data Structure Matfile','Text File (Individual Station)','WSINV (does not work)');
   
    
    if save_menu == 1 %EDI
        
        for is = 1:d.ns
            d.site{is} = [d.site{is},'_edited_',date];
            write_data_edi(d,is);
        end

        d.site = d_orig.site;

    elseif save_menu == 2 %MODEM
        
        if isfield(d,'x')
        
            prompt={'New Data Filename'};
            dlg_title='Filename';
            def={[strtok(d.name,'.'),'_edited.data']};
            num_lines=1;
            dinp = inputdlg(prompt,dlg_title,num_lines,def);
            
            if ~isempty(dinp)
                write_data_modem(dinp{1},d);
            end
            
        else
            
            disp('ModEM Data file not saved! The data you loaded has no mesh associated with it. Save as EDIs and then use M3D to build a mesh');
        end
        
    elseif save_menu == 3 %MATFILE
        
        save([d.name,'_edited'],'d')
        
    elseif save_menu == 4 %TEXT FILE
              
        for ir = 1:4
            ind = ~isnan(d.rho(:,ir,is));
            fid = fopen([d.site{is},'_edited_',d.responses{ir},'_',date,'.txt'],'w');
            A = [d.T(ind) d.rho(ind,ir,is) d.pha(ind,ir,is) d.rhoerr(ind,ir,is) d.phaerr(ind,ir,is)];
            fprintf(fid,'%e %e %e %e %e\r\n',A');
            fclose(fid);
        end
        
    elseif save_menu == 5 %WSINV
        
    end
    
else %QUIT-----------------------------------------------------------------
    irun = 0;
    close(1)
      
end

end

end % END FUNCTION edit_data
%%
function [ax] = get_ax
% sub function to identify axis handles for each data component
% only intended for use with edit_data.m

    h = gcf;
    axesObjs = get(h,'Children');
    get(axesObjs,'Children');
    ax = gobjects([1 4]); % only works for matlab 2013a and later
    
    if numel(findobj(axesObjs,'type','Legend')) == 2 && numel(findobj(axesObjs,'type','axes')) == 5 % plot has xx, xy, yx, yy
        
        ax(1) = axesObjs(4); % rho xx and yy
        ax(2) = axesObjs(7); % rho xy and yx
        ax(3) = axesObjs(2); % pha xx and yy
        ax(4) = axesObjs(5); % pha xy and yx
    
    elseif numel(findobj(axesObjs,'type','Legend')) == 1 && numel(findobj(axesObjs,'type','axes')) == 3 % plot has xy, yx only
        % temporary fix - just leave ax(1) and ax(3) as type graphicsplaceholder
        ax(2) = axesObjs(4); % rho xy and yx
        ax(4) = axesObjs(2); % pha xy and yx
    end

end
%%
function [ax] = set_ax_data(d, ax, is, don)
% sub function to set data plot while editing data. this replaces existing
% plot data with the updated d data. could be editing data or undoing edit
% data
% isgraphics checks are needed to correctly plot for figures with full Z or 
% off-diagonal only Z
% only intended for use with edit_data.m

% d is the data structure
% ax is the axes object containing axes handles
% is is the current station
% don is offset needed to correctly index the axis children

    if isgraphics(ax(1)) % rho xx and yy
        set(ax(1).Children(d.nf+3),'XData',d.T,'YData',d.rho(:,4,is));
        set(ax(1).Children(d.nf*2+4),'XData',d.T,'YData',d.rho(:,1,is));
    end
    if isgraphics(ax(2)) % rho xy and yx
        set(ax(2).Children(don+d.nf+3),'XData',d.T,'YData',d.rho(:,3,is));
        set(ax(2).Children(don+d.nf*2+4),'XData',d.T,'YData',d.rho(:,2,is));
    end
    if isgraphics(ax(3)) % pha xx and yy
        set(ax(3).Children(1),'XData',d.T,'YData',d.pha(:,4,is));
        set(ax(3).Children(2),'XData',d.T,'YData',d.pha(:,1,is));
    end
    if isgraphics(ax(4)) % pha xy and yx
        set(ax(4).Children(1+don),'XData',d.T,'YData',d.pha(:,3,is)+180);
        set(ax(4).Children(2+don),'XData',d.T,'YData',d.pha(:,2,is));
    end

end