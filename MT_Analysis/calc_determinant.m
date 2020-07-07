function [detZ, detZerr] = calc_determinant(Z,Zerr)
%Function to compute the determinant of a given impedance tensor
%
% Usage: [detZ,detZerr] = calc_determinant(Z,Zerr)
%
% Inputs: Z is a nf x 4 x ns complex vector with each column being [xx,xy,yx,yy]
%       Zerr is an nf x 4 x ns complex vector of the same size as Z
%
% Outputs: detZ is a nf x 1 x ns complex vector of the determinant of Z
%       detZerr is an nf x 1 x ns complex vector of the errors.
%               Note: the error propagation is likely incorrect and should
%               be corrected in the future

sze = size(Z);

if length(sze)<2
    sze(2) = 1;
end

if length(sze)<3
    sze(3) = 1;
end

detZ = zeros(sze(1),1,sze(3));
detZerr = zeros(size(detZ));

for i = 1:sze(3) %Loop over stations
    
    %Determinant = ad - bc
    a = Z(:,1,i); b = Z(:,2,i); c = Z(:,3,i); d = Z(:,4,i);
    da = Zerr(:,1,i); db = Zerr(:,2,i); dc = Zerr(:,3,i); dd = Zerr(:,4,i);
    
    detZ(:,1,i) = sqrt(a.*d-b.*c);
    
    %Not sure how to handle the determinant errors. This is done using
    %basic error propagation rules but I do not account for complex numbers
    %so this is wrong.
    detZerr(:,1,i) = sqrt((abs(a.*d).*sqrt((da./abs(a)).^2+(dd./abs(d)).^2)).^2+(abs(b.*c).*sqrt((db./abs(b)).^2+(dc./abs(c)).^2)).^2);
    
    
end