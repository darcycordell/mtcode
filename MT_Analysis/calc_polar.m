function [polar] = calc_polar(d)
%
% Function which calculates the polar diagram (i.e, peanut) for MT data and
% puts the necessary x and y coordinates for plotting into a structure
% called "polar"
%
% Usage: [polar] = calc_polar(d)
%
% Inputs: "d" is standard MT data structure
% Outputs: "polar" contains the x and y coordinates for the peanut diagrams
% for all frequencies, stations, and impedance tensor components in the
% input data
%

rad = 180./pi;              
nrot = 72;

%Loop through all angles and for each angle compute x and y
for irot =1:nrot+1
    rot_ang = (irot-1)*360./nrot;

    % rotate impedances 
    [dr]=rotate_d(d,rot_ang);
    c = cos(rot_ang/rad); s = sin(rot_ang/rad);
    
    polar.x(:,:,:,irot) = c*abs(dr.Z);
    polar.y(:,:,:,irot) = s*abs(dr.Z);
    
end
