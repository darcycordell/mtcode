function [fig] = set_figure_size(n)
%
% A function to set a standard figure size based on the computer's
% screensize.
%
% Usage: [fig] = set_figure_size(n)
%
% Input:
%       n: the figure number
%
% Output:
%       fig: the figure handle
%


if exist('n','var')~=1
    n = 1;
end

screensize=get(groot,'Screensize');
fig=figure(n);
clf
set(fig,'Position',[0.1*screensize(3) 0.1*screensize(4) 0.8*screensize(3) 0.8*screensize(4)])
