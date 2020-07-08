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

    h(i) = plot(NaN,NaN,varargin{i*2});
end

legend(h, legend_string,'AutoUpdate','off');