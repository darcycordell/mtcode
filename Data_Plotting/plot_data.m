function stp = plot_data(d,datatype_flag,dpred,s)
%Function to plot MT data
%
% Usage: plot_data(d,datatype_flag) OR plot_data(d,datatype_flag,dpred)
%
% "d" is a data structure containing all relevant MT data
%
%datatype_flag (OPTIONAL): 
%       1 = apparent resistivity and phase (default)
%       2 = real and imaginary impedances
%       3 = tipper
%
% "dpred" is an OPTIONAL data structure which includes another set of MT
% data. If you want to compare two data sets (e.g. observed data and
% predicted data) then include dpred in the function call.
%
% "s" is in OPTIONAL input. it is the structure output from
% detailed_statistics or detailed_statistics_2D. In the first case,
% statistics are computed from impedances and in the latter, statistics are
% computed from rho and phase as in the Mackie 2D inversion.

if ~exist('datatype_flag','var')
    datatype_flag = 1;
end

if exist('dpred','var') && ~exist('s','var')
    disp('No input stats structure detected. Misfit calculated from impedances by default')
    s = detailed_statistics(d,dpred);
end

if isfield(d,'loc')
    d = set_map_projection(d);
    [L] = load_geoboundary_file_list;
end

%Choose site range
site_range={[num2str(1),':',num2str(1),':',num2str(d.ns)]};
answer=char(inputdlg({'Enter site range (ex. : 1:1:13)'},'Plot sites',1,site_range));
%The sites to plot (stp = site to plot)
stp=str2num(answer);  %#ok<*ST2NM>

%Check if user input is valid
if isempty(stp)
stp = 1;
disp('Invalid Entry: Plotting Station 001 as Default');
end

if max(stp) > d.ns || min(stp) < 1
stp = 1;
disp('Invalid Entry: Plotting Station 001 as Default')
end


%Loop over sites to plot
for is = stp
    close all;
    
    if datatype_flag == 1 %Plot apparent resistivity and phase
        set_figure_size(1);
        plot_rho_pha(d,is)
        if exist('dpred','var')
            plot_rho_pha(dpred,is)
        end
        figure_name = 'rho_pha_';
    elseif datatype_flag == 2 %Plot impedances
        set_figure_size(1);
        plot_impedance(d,is)
        if exist('dpred','var')
            plot_impedance(dpred,is)
        end
        figure_name = 'impedance_';
    elseif datatype_flag == 3 %Plot tipper
        set_figure_size(1);
        plot_tipper(d,is)
        if exist('dpred','var')    
            plot_tipper(dpred,is)
        end
        figure_name = 'tipper_';
    end

    %If no data has been plotted then don't plot stations or rms (e.g. if
    %person chose "plot tipper" when there is no tipper).
    if ~isempty(findobj(gcf,'type','axes'))

        if ~strcmp('ZXX',d.responses)
    %     if isempty(find(ismember(d.responses,'ZXX'),1)) %Check if diagonals exist and set subplots accordingly
            subplot(2,2,2)
        else
            subplot(2,3,3)
        end

        if isfield(d,'loc')
            %Plot station locations in map view
            plot_topo(d,1);    
            plot_geoboundaries_geoshow(L)
            geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
            geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
            axis(d.lim)
            xlabel('Longitude'); ylabel('Latitude')
            box on
        else % assumed 2-D inversion
            plot(d.y,d.z,'kv','markerfacecolor','k')
            hold on
            plot(d.y(is),d.z(is),'rv','markerfacecolor','r')
            set(gca,'ydir','reverse')
            xlabel('Distance (km)')
            ylabel('Depth (km)')
            set(gca,'DataAspectRatio',[1 1 1])
            offset = 0.2 * diff(get(gca,'ylim'));
            text(d.y,d.z - offset,strrep(d.site,'_','\_'),'rotation',90)

        end

        %If predicted data is supplied then plot rms for each site
        if exist('dpred','var')
            plot_misfit_rms_by_site(d,dpred,is,s) %Plot rms
        end

    end

    if exist('dpred','var') %If predicted data is supplied alter filenames and folders accordingly before saving figure
        print_figure(['compare_responses_',dpred.niter],[figure_name,num2str(is,'%03.0f'),'_',d.site{is}]); %Save figure
    else
        print_figure(['data_response'],[figure_name,num2str(is,'%03.0f'),'_',d.site{is}]); %Save figure
    end

end

