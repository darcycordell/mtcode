function [ssqZ, ssqZerr] = calc_ssq(Z,Zerr)
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

ssqZ = zeros(sze(1),1,sze(3));
ssqZerr = zeros(size(ssqZ));

for i = 1:sze(3) %Loop over stations
    
    %Determinant = ad - bc
    a = Z(:,1,i); b = Z(:,2,i); c = Z(:,3,i); d = Z(:,4,i);
    da = Zerr(:,1,i); db = Zerr(:,2,i); dc = Zerr(:,3,i); dd = Zerr(:,4,i);
    
    ssqZ(:,1,i) = sqrt(0.5*(a.^2+b.^2+c.^2+d.^2));
    
    %Not sure how to handle the determinant errors.
    %One option is to simply take the average of the relative errors of
    %all the tensor components. However, this often results in
    %disproportionately large errors from the diagonal components
    %ssqZerr(:,1,i) = 0.25*abs(ssqZ).*(da./abs(a)+db./abs(b)+dc./abs(c)+dd./abs(d));

    %A third option is to take the average of the off-diagonal components
    %only
    ssqZerr(:,1,i) = 0.5*abs(ssqZ).*(db./abs(b)+dc./abs(c));
    
end