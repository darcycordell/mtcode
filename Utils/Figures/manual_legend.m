function h = manual_legend(varargin)
%
% Function to produce a manual legend that does not depend on the current
% plot.
%
% Usage: manual_legend(varargin)
%
% Inputs: varargin is a list of linetype and label pairs as strings
%
% Example:
%   manual_legend('Data #1','or','Data #2','-k','Data #3','*m')
%
%   This will produce a legend which shows red circles, black lines, and
%   magenta stars with their corresponding labels.
%

numel = length(varargin)/2;

if round(numel)-numel ~= 0
    error('Error on Input of Manual Legend: Must have linetype and legend label pairs as strings')
end

h = zeros(numel,1); legend_string = {};
for i = 1:numel
    legend_string{i} = varargin{i*2-1};

    if ischar(varargin{i*2})
        h(i) = plot(NaN,NaN,varargin{i*2}); hold on
    elseif isnumeric(varargin{i*2}) && length(varargin{i*2})==3
        h(i) = plot(NaN,NaN,'.','Color',varargin{i*2}); hold on
    else
        error('You must input either a MATLAB-recognized string specifying the color/linetype (e.g. -k) or a 3x1 vector of colors (e.g. [0 0.5 0.2])')
    end
end

legend(h, legend_string,'AutoUpdate','off');