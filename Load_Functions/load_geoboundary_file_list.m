function [L] = load_geoboundary_file_list(varargin)
%
% This function loads a text file containing geoboundary text file names.
% The text file name is defined in user_defaults.
%
% Outputs: L is a data structure containing all the information about the
% various geological or geographic boundaries. Each row of L is a different
% "type" of geographic or geological datawith an associated marker and
% line specification. 
% For example, city locations could be one "type" plotted as black circles 'ok'
% and these locations are in the first row of L.
% Rivers could be another "type" in the second row specified by blue lines
% '-b' in the second column. "L" contains all this information stored as 
% self-descriptive structure handles.
%
% L.line contains the geographic information itself (e.g. lat, long, elevation)
%
% For plotting lines, L.line is a cell containing latitude/longitude pairs of
% points. Each cell represents a new line. For example, if your dataset
% contains 3 different rivers, then you would have a 3x1 cell where each
% cell contains the lat/long points of each river. This ensures that when
% plotting the rivers as lines, MATLAB does not connect the end points of
% different rivers together.
%
% The text file should have the geoboundary file names in a list, with the
% line specifiers given in a second column. The code stops reading the file
% once it hits a #. So you can add extra lines after a #
% For example:
%
%   Lakes.txt, fill, -k, 1, b, 0.5
%   Roads.txt, line, '-k', 3
%   Cities.txt, marker, 'ro', 4
%
%   #
%   Everything after the '#' is not read and can be used for comments
%
% In the above example Lakes.txt is plotted as a filled polygon with a black
% outline with a linewidth of 1. The polygon is filled blue with a 0.5 transparency
% The Roads.txt file is plotted as a black line with a width of 3
% Cities are plotted as red open circles with size of 4.
%
% Each text file (e.g. Lakes.txt, Roads.txt, etc.) in the list must have the following
% format: first column = longitude, the second column =
% latitude. Optional third column = elevation. If third column does not exist, 
% then elevation is defaulted to 0.
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

u = user_defaults;

%If no geofile is specified in user_defaults, then the
if strcmp(u.geofile,'none') && nargin==0
    L = 0;
    return
end

if nargin == 1
    u.geofile = varargin{1};
end

fid=fopen(u.geofile);

if fid == -1
    error('Error: Geoboundary file does not exist on path!')
end
% object = struct([]);
L = struct([]);
i = 1;
%First read in the list of text files contained within the geofile
while 1
    line_str = fgetl(fid);

    if line_str == -1
        break
    end
    
    %Anything after '#' is not read
    if strcmp(line_str,'#') || strcmp(line_str,'')
        break
    end
    
    C = strsplit(line_str,',');
    
    if strcmp(strtrim(C{2}),'line') %%%%%%%%%%%% LINE
        
        if numel(C) == 2 % only filename and line specified
            L(i).file = strtrim(C{1}); %File name for ith type 
            L(i).type = strtrim(C{2});
            L(i).linestyle = '-'; % default
            L(i).color = 'k';
            L(i).linewidth = 0.5; % default                   
            L(i).marker = 'none';
            L(i).markersize = 0.5; % default
        elseif numel(C) == 4 % a linestyle and linewidth are specified
            L(i).file = strtrim(C{1}); %File name for ith type
            L(i).type = strtrim(C{2});
            tmp = strtrim(C{3}); % third argument
            [L(i).linestyle, L(i).color] = strtok(tmp,{'y','m','c','r','g','b','w','k'}); % no marker type, so this works
            L(i).linewidth = str2double(C{4});
            L(i).marker = 'none';
            L(i).markersize = 0.5; % default                                             
        elseif numel(C) == 5 % a linestyle, linewidth, marker, and markersize are specified
            L(i).file = strtrim(C{1}); %File name for ith type
            L(i).type = strtrim(C{2});
            tmp = strtrim(C{3}); % third argument
            if contains(tmp,{'o','+','*','.','x','s','d','^','v','>','<','p','h'})
%             [L(i).linestyle, ~] = strtok(tmp,{'y','m','c','r','g','b','w','k'});
                ts = regexp(tmp,{'y','m','c','r','g','b','w','k'},'match');
                L(i).color = char(string(ts(~cellfun('isempty',ts))));
                st = strsplit(tmp,{'y','m','c','r','g','b','w','k'});
                L(i).linestyle = st{1};
                L(i).marker = st{2};
                L(i).linewidth = str2double(C{4});
                L(i).markersize = str2double(C{5});
            else
                error(['Unrecognized marker symbol on line ',num2str(i),' in ',u.geofile])
            end
            
%             L(i).linestyle = tmp(1);
%             L(i).color = tmp(2);
%             if numel(tmp) == 2 %  marker                   
%                 L(i).marker = 'none';
%             elseif numel(tmp) == 3
%                 L(i).marker = tmp(3);    
%             end
%             L(i).linewidth = str2double(C{4});
%             if numel(C) == 5
%                 L(i).markersize = str2double(C{5});
%             else
%                 L(i).markersize = 1;
%             end           
        else
            error(['Format of line ',num2str(i),' in ',u.geofile,' is incorrect'])
        end
        
    elseif strcmp(strtrim(C{2}),'marker') %%%%%%%%% MARKER
        
        if numel(C) == 2 % only filename and marker specified
            L(i).file = strtrim(C{1}); %File name for ith type 
            L(i).type = strtrim(C{2});
            L(i).linestyle = 'none'; % default
            L(i).color = 'k';
            L(i).linewidth = 0.5; % default                   
            L(i).marker = 'o';
            L(i).markersize = 8; % default
        elseif numel(C) == 4
            L(i).file = strtrim(C{1}); %File name for ith type 
            L(i).type = strtrim(C{2});
            L(i).linestyle = 'none'; % default
            tmp = strtrim(C{3}); % third argument
%             [~, L(i).color] = strtok(tmp,{'y','m','c','r','g','b','w','k'});
            L(i).color = tmp(1);
            L(i).linewidth = 0.5; % default  
            L(i).marker = tmp(2);
            L(i).markersize = str2double(C{4}); % default           
        else
            error(['Format of line ',num2str(i),' in ',u.geofile,' is incorrect'])           
        end
                             
    elseif strcmp(strtrim(C{2}),'fill') %%%%%%%%%%%%%%%%% FILL
         
        if numel(C) == 2 % only filename and fill specified
            L(i).file = strtrim(C{1}); %File name for ith type 
            L(i).type = strtrim(C{2});
            L(i).linestyle = '-'; % default
            L(i).color = 'k';
            L(i).linewidth = 0.5; % default
            L(i).facecolor = 'b'; % default
            L(i).facealpha = 1; % default
        elseif numel(C) == 6
            L(i).file = strtrim(C{1}); %File name for ith type 
            L(i).type = strtrim(C{2});
            tmp = strtrim(C{3}); % third argument           
            [L(i).linestyle, L(i).color] = strtok(tmp,{'y','m','c','r','g','b','w','k'});
            L(i).marker = 'none';
            L(i).markersize = 0.5; % default
            L(i).linewidth = str2double(C{4});
            L(i).facecolor = strtrim(C{5});
            L(i).facealpha = strtrim(C{6});        
        else
            error(['Format of line ',num2str(i),' in ',u.geofile,' is incorrect'])
        end
        
        
    else
        error(['Format of line ',num2str(i),' in ',u.geofile,' is incorrect'])
            
    end
      
        

            
    i = i+1;

end

fclose(fid);

%Initialize L structure
% L = struct([]);
% L(length(file)).line = {};
% % L(length(file)).linespec = [];
% % L(length(file)).linewidth = [];
% 
% L(length(file)).marker = {};
% L(length(file)).color = [];
% L(length(file)).linestyle = [];
% L(length(file)).linewidth = [];
% L(length(file)).markersize = [];

%Now loop through the files and load the geographic coordinates
for i = 1:length(L)

    fid = fopen(L(i).file);

    if fid == -1
        error(['Specified geoboundary file ',file{i},' does not exist on path!'])
    end

    line_str=fgetl(fid);
    boundary_count=1;
    line_count=1; clear line
    %Loop through file to separate out different sections of a given type
    %into different cells. For example, if your file contains geographic
    %coordinates for multiple rivers, then separate the different rivers by
    %an 'X' in your file and each river will be put in a different cell.
    while 1
        if line_str == -1
            break
        end
        if ~isempty(line_str)
            if line_str(1)~='#'
                if line_str(1) ~= 'X' && line_str(1) ~= '>'
                    line{boundary_count}(:,line_count)=sscanf(line_str,'\t%f%f');
                    line_count=line_count+1;
                else
                    boundary_count=boundary_count+1;
                    line_count=1;
                end
            end
        end

        line_str=fgetl(fid);
    end
    fclose(fid);
    
    L(i).line = line;
%     L(i).linespec = linespec{i};
%     L(i).linewidth = linewidth(i);
    
%     L(i).marker = marker{i};
%     L(i).color = color{i};
%     L(i).linestyle = linestyle{i};
%     L(i).linewidth = linewidth{i};
%     L(i).markersize = markersize{i};

end