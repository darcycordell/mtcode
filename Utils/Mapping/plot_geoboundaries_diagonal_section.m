function plot_geoboundaries_diagonal_section(d,ew_locations,ns_locations)
% Function to plot geoboundaries on a diagonal section.
% Only tested for well tracks. Reads in text file with longitude,latitude,
% elevation, then projects well track onto the diagonal slice.
%
% Usage: plot_geoboundaries_diagonal_section(d,ew_locations,ns_locations)
%
% "d" is the data structure (necessary for conversion to lat-long)
% ew_locations are the EW locations of the diagonal trace
% ns_locations are the NS locations of the diagonal trace

%%

if ~strcmp(d.name,'None')
[L] = load_geoboundary_file_list;

if ~isstruct(L)
    return
end

[x1, ind] = min(ew_locations);
y1 = ns_locations(ind);

[x2, ind] = max(ew_locations);
y2 = ns_locations(ind);
s= (y2-y1)/(x2-x1);

% make these more finely sampled so projection does not depend on cell size
ns_locations = linspace(ns_locations(1),ns_locations(end),1000);
ew_locations = linspace(ew_locations(1),ew_locations(end),1000);

for count = 1:length(L)
    
    line = L(count).line;
    
    for i=1:length(line)

        line_geo=line{i};

        [line_utm(1,:),line_utm(2,:)] = geo2utm(line_geo(1,:),line_geo(2,:),d.origin(2),d.origin(1));

        line_utm(1,:) = line_utm(1,:)-500000;

        for j = 1:size(line_geo,2) % find indices of projection
            [~,I(j)] = min(sqrt( ((line_utm(1,j)-ew_locations).^2) + ((line_utm(2,j)-ns_locations).^2) ));              
        end

        if s<0 % if northwest to southeast slice
            dx = abs(ew_locations(I) - min(ew_locations)); % find distance to "zero" point on diagonal slice
            dy = abs(ns_locations(I) - max(ns_locations));
        else % if southwest to northeast slice
            dx = abs(ew_locations(I) - min(ew_locations)); % not thoroughly tested
            dy = abs(ns_locations(I) - min(ns_locations));
        end

        dist = sqrt( (dx.^2) + (dy.^2) ); % distance to "zero" point
        if length(line_geo(:,1))>2
            figure(2)
%             plot(dist./1000,line_geo(3,:),L(count).linespec,'LineWidth',L(count).linewidth);
%             plot(dist./1000,line_geo(3,:),L(i).linestyle,'color',L(i).color,'LineWidth',L(i).linewidth,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on

            if strcmp(L(count).type,'line')
                plot(dist./1000,line_geo(3,:),L(count).linestyle,'LineWidth',L(count).linewidth,'color',L(count).color,'marker',L(count).marker,'markersize',L(count).markersize,'markerfacecolor',L(count).color); hold on
            elseif strcmp(L(count).type,'marker')
                plot(dist./1000,line_geo(3,:),'linestyle',L(count).linestyle,'LineWidth',L(count).linewidth,'color',L(count).color,'marker',L(count).marker,'markersize',L(count).markersize,'markerfacecolor',L(count).color); hold on
            elseif strcmp(L(count).type,'fill')
                patch(dist./1000,line_geo(3,:),L(count).facecolor,'linestyle',L(count).linestyle,'linewidth',L(count).linewidth,'edgecolor',L(count).color,'facealpha',L(count).facealpha); hold on
            else
                warning(['Unable to plot geo boundaries from ',L(count).file])
            end
            
            
            
            
            
            

        end
        clear line_utm I        
    end

end

end

end