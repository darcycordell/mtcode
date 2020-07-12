function P = plot_rectangular_fault_plane(strike, dip, width, length,location)
% Function which plots a rectangular fault plane in 3-D space
%
% Inputs
%
% strike: strike of fault in degrees (right-hand rule)
% dip: dip of the fault in degrees (right-hand rule)
% width: width of the fault specified in axis units (e.g. if you are
% plotting in lat-long, then this should be in lat-long degrees)
% length: length of the fault specified in axis units
% location: The location of the fault in 3-D space specified as a vector [Px, Py, Pz]
%   The location specifies the corner of the prism opposite to the quadrant
%   of the strike. E.g. for strike = 45 (NE quadrant), the location specifies 
%   the SW corner of the rectangle; for strike = 75 (SE quadrant), the
%   location specifies the NW corner of the rectangle.
%
% Note the "width" is defined as the dipping side of the rectangular plane
% while the "length" is the side parallel to the surface. Width can be
% larger than length.
%
% Outputs
%
% P is a 4x3 vector specifying the corners of the rectangle A,B,C,D

A = location; % Pivot point at the surface

H = width*sind(strike+90);
V = width*cosd(strike+90);
Z = -width*sind(dip);

B = [length*sind(strike) length*cosd(strike) 0]+A;
C = [length*cosd(90-strike)+H length*sind(90-strike)+V Z]+A;
D = [width*sind(strike+90) width*cosd(strike+90) -width*sind(dip)]+A;

P = [A; B; C; D];

figure(1); p = patch(P(:,1),P(:,2),P(:,3),'red');
xlabel('EW Distance (km)');ylabel('NS Distance (km)');zlabel('Elevation (km b.s.l.)'); %axis equal
axis equal


end







