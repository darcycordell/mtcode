function [dim] = calc_dim_parameters(Z,rot_ang)
%
% Function which calculate dimensionless parameters. See Simpson and Bahr
% (2004) for details or Martyn's 699 notes.
%
% Usage: [dim] = calc_dim_parameters(Z,rot_ang)
%
% Inputs:
%       Z is a (nf x 4) matrix of impedances for a single MT station
%       rot_ang is the angle to which the impedance data are rotated to
%
% Outputs:
%       dim is a structure containing the necessary dimensionless
%       parameters: s1, s2, d1, d2, Swift skew, Bahr skew
%       It also includes the Swift angle
%

disp('Calculating Dimensionless Parameters')
rad = 180./pi;   
dim.s1 = Z(:,1)+Z(:,4);
dim.d1 = Z(:,1)-Z(:,4);
dim.s2 = Z(:,2)+Z(:,3);
dim.d2 = Z(:,2)-Z(:,3);
dim.swift_skew = abs(dim.s1)./abs(dim.d2);
dim.a =2*real(dim.s2.*conj(dim.d1))./(abs(dim.d1).*abs(dim.d1)-abs(dim.s2).*abs(dim.s2));
dim.alpha = rot_ang+rad*0.25*atan(dim.a);

c1 = real(dim.d1).*imag(dim.s2)-real(dim.s2).*imag(dim.d1);
c2 = real(dim.s1).*imag(dim.d2)-real(dim.d2).*imag(dim.s1);
dim.bahr_skew = sqrt(abs(c1-c2))./abs(dim.d2);

