function plot_geoboundaries(L,varargin)
%
% Apr. 2020 version allows lines, markers, and filled polygons
%
% Function which reads in a text file containing a list of text files and a
% line specification. The file containing the list of text files must be on
% the MATLAB path and must be in the following format:
%
% Valid entries:
% 
% for Markers:
% filename, type
% filename, type, color+marker, markersize
% 
% for Lines:
% filename, type
% filename, type, linestyle+color, linewidth
% filename, type, linestyle+color+marker, linewidth, markersize
% 
% for Fills:
% filename, type
% filename, type, linestyle+color, linewidth, facecolor, facealpha
%
% if only filename and type are specified, default values are used for the
% other parameters
%
% filename = name of file containing longitude and latitude
% type = marker, line, or fill
% color, marker, and linestyle must use the default MATLAB symbols
% linewidth and markersize must be positive values
% facealpha is the transparancy of the golled polygon, a number between 0
% and 1.
%
% Examples:
% 
% for Markers:
% points.txt, marker
% points.txt, marker, kd, 12
% 
% for Lines:
% lines.txt, line
% lines.txt, line, --k, 2
% lines.txt, line, -.ko, 2, 8
% 
% for Fills:
% polygons.txt, fill
% polygons.txt, fill, -k, 2, b, 0.5
%
% Each file (points.txt, lines.txt, etc.) must also be on the MATLAB path.
% These files contain longitudes and latitudes specifying the location of
% points to plot. The line specification tells whether the locations should
% be plotted as markers,lines, or filled polygons.
%
% Each text file (points.txt, lines.txt, etc.) must have the following
% format: first column = longitude, the second column =
% latitude. Optional third column = elevation in m.b.s.l (negative numbers above sea level).
% If third column does not exist, then elevation is 0.
% If you have multiple outlines (e.g. multiple lakes), you can
% separate them using an X in the first column.
%
% Example:
%
% 	-170.654762689	-36.308221971
% 	-170.655393971	-36.305359311
% 	-170.655173541	-36.303431968
% 	-170.651991609	-36.303275250
% X
% 	-170.566124141	-36.025318497
% 	-170.562805249	-36.027113793
% 	-170.560039018	-36.028777844
%
%
% Usage: plot_geoboundaries(origin,z) OR plot_geoboundaries
%
% If you are plotting the geoboundaries in map coordinates (e.g. latitude and longitude), then no input
% is needed. The geoboundary is plotting using m_map.
% 
% However, if you are plotting the geoboundaries in model
% coordinates (e.g. x and y), then it is necessary to provide additional
% information to convert the lat-long coordinates to model coordinates. You
% need an origin (e.g. d.origin from data structure) and a z coordinate
% (e.g. d.z from the data structure) to geo-reference the lat-longs in
% model space
%

if ~isstruct(L)
    return
end

if isempty(get_projection)
    disp('Map projection not yet initialized. Defaults being set')
    m_proj('mercator','long',[-180 180],'lat',[-90 90]);
end

for i = 1:length(L)

    line = L(i).line;
    
    for il=1:length(line)

        if ~isempty(line{il})
            if isempty(varargin)
                
                if strcmp(L(i).type,'line')
                    m_plot(line{il}(1,:),line{il}(2,:),L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                
                elseif strcmp(L(i).type,'marker')
                    m_plot(line{il}(1,:),line{il}(2,:),'linestyle',L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    
                elseif strcmp(L(i).type,'fill')
                    m_patch(line{il}(1,:),line{il}(2,:),L(i).facecolor,'linestyle',L(i).linestyle,'linewidth',L(i).linewidth,'edgecolor',L(i).color,'facealpha',L(i).facealpha); hold on
                    
                else
                    warning(['Unable to plot geo boundaries from ',L(i).file])
                end
                

            else % varargin only to enable plotting in model coordinates
                origin = varargin{1};
                zorigin = varargin{2};
                [x,y] = geo2utm(line{il}(1,:),line{il}(2,:),origin(2),origin(1));
                x = x-500000;
                
                % plot3(x/1000,y/1000,min(zorigin)*ones(1,length(x))/1000,L(i).linestyle,'color',L(i).color,'LineWidth',L(i).linewidth,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on

                if length(line{il}(:,1)) == 3
                    %The z-coordinate is negative for plotting purposes
                    %because all inversion codes use z-positive downwards
%                     plot3(x/1000,y/1000,line{il}(3,:),L(i).linespec,'LineWidth',L(i).linewidth,'MarkerFaceColor',L(i).linespec(end)); hold on 
%                     plot3(x/1000,y/1000,line{il}(3,:),L(i).linestyle,'color',L(i).color,'LineWidth',L(i).linewidth,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    if strcmp(L(i).type,'line')
                        plot3(x/1000,y/1000,line{il}(3,:),L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    elseif strcmp(L(i).type,'marker')
                        plot3(x/1000,y/1000,line{il}(3,:),'linestyle',L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    elseif strcmp(L(i).type,'fill')
                        patch(x/1000,y/1000,line{il}(3,:),L(i).facecolor,'linestyle',L(i).linestyle,'linewidth',L(i).linewidth,'edgecolor',L(i).color,'facealpha',L(i).facealpha); hold on
                    else
                        warning(['Unable to plot geo boundaries from ',L(i).file])
                    end
                else
                    %The z-coordinate is negative for plotting purposes
                    %because all inversion codes use z-positive downwards
%                     plot3(x/1000,y/1000,min(zorigin)*ones(1,length(x))/1000,L(i).linespec,'LineWidth',L(i).linewidth,'MarkerFaceColor',L(i).linespec(end)); hold on
                    if strcmp(L(i).type,'line')
                        plot3(x/1000,y/1000,min(zorigin)*ones(1,length(x))/1000,L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    elseif strcmp(L(i).type,'marker')
                        plot3(x/1000,y/1000,min(zorigin)*ones(1,length(x))/1000,'linestyle',L(i).linestyle,'LineWidth',L(i).linewidth,'color',L(i).color,'marker',L(i).marker,'markersize',L(i).markersize,'markerfacecolor',L(i).color); hold on
                    elseif strcmp(L(i).type,'fill')
                        patch(x/1000,y/1000,min(zorigin)*ones(1,length(x))/1000,L(i).facecolor,'linestyle',L(i).linestyle,'linewidth',L(i).linewidth,'edgecolor',L(i).color,'facealpha',L(i).facealpha); hold on
                    else
                        warning(['Unable to plot geo boundaries from ',L(i).file])
                    end
                end
            end
        end

    end
end

% if isempty(varargin)
%     m_grid('box','fancy','xlabeldir','end','tickdir','in');
% end


end% END main