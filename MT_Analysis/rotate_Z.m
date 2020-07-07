function Z_rot = rotate_Z(Z,rot)
%
% Function which rotates the 2 x 2 impedance tensor for a particular
% frequency and location
%
% Usage: Z_rot = rotate_Z(Z,rot)
%
% Inputs: "Z" is a 1 x 4 vector of complex impedance: [Zxx Zxy Zyx Zyy]
%       "rot" is an angle (in degrees) to rotate the data clockwise
%       (negative angles is counter-clockwise)
%
% Outputs: "Z_rot" is the rotated 1 x 4 vector of complex impedance [Zxx Zxy Zyx Zyy]

c = cosd(rot);      s = sind(rot);

R = [ c, -s ; s, c]; %Rotation matrix 

Z_rot = R'*[Z(1) Z(2); Z(3) Z(4)]*R;

Z_rot = [Z_rot(1,1) Z_rot(1,2) Z_rot(2,1) Z_rot(2,2)];