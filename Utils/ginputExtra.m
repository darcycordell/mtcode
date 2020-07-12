function [x y] = ginputExtra(n,pointflag,booText)
% INPUT
% n: Number of points to plot
% pointflag = 1: plots lines between points to form polygon (default)
% pointflag = 0: plots points only (no lines connecting points)
% booText: Boolean (default false) command to display point number in
% the plot.


% Author: Lasse Nørfeldt (Norfeldt) 
% Date: 2012-04-09

if nargin ==3
    bText = booText;
else
    bText = false;
end

if nargin == 1
    pointflag = 1;
end

if pointflag == 1
    pointorline = 'r-';
else
    pointorline = 'o';
end

    
H = gca;
set(H, 'YLimMode', 'manual'); set(H, 'XLimMode', 'manual');
set(H, 'YLim', get(H,'YLim')); set(H, 'XLim', get(H,'XLim'));

numPoints = n; xg = []; yg = [];
for i=1:numPoints
    [xi yi click] = ginput(1);
    if click==3
        break
    end
    xg = [xg xi]; yg = [yg yi];
    if i == 1
        hold on;
        plot(H, xg(i),yg(i),'o');
        if bText text(xg(i),yg(i),num2str(i),'FontSize',14); end
    else
        plot(xg([i-1:i]),yg([i-1:i]),pointorline);
        if bText text(xg(i),yg(i),num2str(i),'FontSize',14); end
    end 
end
hold off;

x = xg; y = yg;