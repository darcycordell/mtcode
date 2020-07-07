function tip_rot = rotate_tip(tip,rot)
%
% Function which rotates the tipper transfer function for a particular
% tipper value
%
% Usage: tip_rot = rotate_tip(tip,rot)
%
% Inputs: "tip" is a 1 x 2 vector of complex [Tzx Tzy]
%       "rot" is an angle (in degrees) to rotate the data clockwise
%       (negative angles is counter-clockwise)
%
% Outputs: "tip_rot" is the rotated 1 x 2 vector of complex [Tzx Tzy]

c = cosd(rot);      s = sind(rot);

R = [ c, -s ; s, c]; %Rotation matrix 

tip_rot = (R*tip.').';