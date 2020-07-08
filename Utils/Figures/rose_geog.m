function rose_geog(E,nbins,rmax,col)
%
% Function to plot a CORRECT rose diagram (with the proper polygons filled
% in). The original rose diagram function often did not fill the polygons
% correctly.
%
% Improved by MJU.
%
% Based on 
% http://www.mathworks.com/matlabcentral/fileexchange/33664-earthrose/content/earth_rose.m
% Cameron Sparr cameronsparr@gmail.com
%
% Usage: rose_geog(E,nbins,rmax,col)
%
% Inputs:
%       E: the angle data to bin in a circular histogram (in degrees)
%       nbins: the number of bins to bin the data
%       rmax: The maximum tick to plot on the radial scale (e.g. the radius
%              of the rose diagram)
%       col: The color to fill the polygons
%

% edited BL 2020 - prepare azimuths separately for phase tensor and
% induction arrows in invoking functions - focus on ploting in this
% function
E = E*pi/180;  % Convert to radians

if nargin > 1  
    % Following added by MJU to fill sectors
    polar([pi/2, pi/2],[0., rmax]) % Plot vector to control radial scale
    set(findobj(gca,'Type','line'),'Color','k');
    hold on
    polar([-pi/2, -pi/2],[0., rmax])
    set(findobj(gca,'Type','line'),'Color','k');
    [tout, rout] = rose(E,nbins);
    polar(tout, rout);
    [xout, yout] = pol2cart(tout, rout);
    set(gca, 'nextplot', 'add');
%     patch(xout, yout, col); % this doesn't work, need to plot each bin separately
    for i = 1:nbins
        patch(xout(4*(i-1)+1:4*i), yout(4*(i-1)+1:4*i), col);
    end
else
    rose(E);
end

hHiddenText = findall(gca,'type','text');
Angles = 0 : 30 : 330;
hObjToDelete = zeros( length(Angles)-4, 1 );
k = 0;
for ang = Angles
   hObj = findall(hHiddenText,'string',num2str(ang));
   switch ang
   case 0
      set(hObj,'string','E','HorizontalAlignment','Left');
   case 90
      set(hObj,'string','N','VerticalAlignment','Bottom');
   case 180
      set(hObj,'string','W','HorizontalAlignment','Right');
   case 270
      set(hObj,'string','S','VerticalAlignment','Top');
   otherwise
      k = k + 1;
      hObjToDelete(k) = hObj;
   end
end
delete( hObjToDelete(hObjToDelete~=0) );

end