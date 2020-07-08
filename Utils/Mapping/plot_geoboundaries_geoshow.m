function plot_geoboundaries_geoshow(L)
% Function which plots the geoboundaries stored in the "L" data structure
% using the "geoshow" function in the MATLAB mapping toolbox. For use with
% m_map or UTM model space plotting, see "plot_geoboundaries"
%
% Usage: plot_geoboundaries_geoshow(L)
%
% Inputs: 
%   "L" is a data structure containing all the information about the
%   various geological or geographic boundaries.
%
%   See load_geoboundary_file_list.m for more info on "L"
%


if ~isstruct(L)
    return
end

for i = 1:length(L)

    line = L(i).line;
    
    for il=1:length(line)

        if ~isempty(line{il})
            
            if strcmp(L(i).type,'line')
                geoshow(line{il}(2,:),line{il}(1,:),'DisplayType','line','linestyle',L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
            elseif strcmp(L(i).type,'marker')
                geoshow(line{il}(2,:),line{il}(1,:),'DisplayType','multipoint','Marker',L(i).marker,'markeredgecolor',L(i).color,'markerfacecolor',L(i).color,'MarkerSize',L(i).markersize); hold on
            elseif strcmp(L(i).type,'fill')
                if ispolycw(line{il}(1,:),line{il}(2,:)) % check that the fill is defined clockwise. CCW plots as a hole
                    fill_lat = line{il}(2,:);
                    fill_lon = line{il}(1,:);
                else
                    [fill_lon, fill_lat] = poly2cw(line{il}(1,:), line{il}(2,:));
                end
                geoshow(fill_lat,fill_lon,'DisplayType','polygon','facecolor',L(i).facecolor,'linestyle',L(i).linestyle,'linewidth',L(i).linewidth,'edgecolor',L(i).color,'facealpha',L(i).facealpha); hold on
%                 geoshow(line{il}(2,:),line{il}(1,:),'DisplayType','polygon','facecolor',L(i).facecolor,'linestyle',L(i).linestyle,'linewidth',L(i).linewidth,'edgecolor',L(i).color,'facealpha',L(i).facealpha); hold on
            else
                warning(['Unable to plot geo boundaries from ',L(i).file])
            end
            
        end

    end
end


end% END main