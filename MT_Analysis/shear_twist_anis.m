function Z_dist=shear_twist_anis(Z, shear,twist,anis)
%
%Z_prime=shear_twist(Z,shear,twist)
%
%    shear_twist(Z,shear,twist) distorts an impedance
%    tensor Z by applying shear and twist according
%    to Groom & Bailey (1989) to the tensor.
%    The impedance tensor has to be a 2x2xN array
%    while shear and twist need to be in degrees.
Zmat = zeros(2,2,size(Z,1));
for i = 1:size(Z,1)
    Zmat(:,:,i) = [Z(i,1) Z(i,2); Z(i,3) Z(i,4)];
end

% Calculating e and t
e = tand(shear);
t = tand(twist);

% Calculating S and T
S = 1/sqrt(1+e^2)*([1 e; e 1]);
T = 1/sqrt(1+t^2)*([1 -t; t 1]);
A = 1/sqrt(1+anis^2)*([1+anis 0; 0 1-anis]);

% Calculating distorted Z
Z_prime = zeros(2,2,size(Zmat,3));
for j = 1:size(Zmat,3)
    Z_prime(:,:,j) = T*S*A*Zmat(:,:,j);
end

Z_dist(:,1) = squeeze(Z_prime(1,1,:));
Z_dist(:,2) = squeeze(Z_prime(1,2,:));
Z_dist(:,3) = squeeze(Z_prime(2,1,:));
Z_dist(:,4) = squeeze(Z_prime(2,2,:));