% ----------------------------------------------------------
%
% GEN_JFile.m
%
% just for strike test -- create a simple dat file w/ a 
% ridiculous impedance tensor of well defined strike (ws)
%
% well - by time not so ridicolous any more:
% 
% 1) choice of basic 2-D tensor
% 2) random sinusoidal variation of this tensor (real and 
%    imaginary parts seperately) as a function of period
% 3) multiplication by a distortion matrix 
% 4) rotation into new coordinate system
% 5) add (multiply) random relative error of chosen floor
% 6) rescale in order to reflect ~stable resistivities
% 7) output J-format file 
%
% - all paramter hardwired.
%
% ----------------------------------------------------------

% ----------------------------------------------------------
% ONCE AND FOR ALL:
%
%   To prevent any further confusion w/ rotations:
%
%   -> Two definitions of rotation matrices:
%
%   (1) R_f = (c -s; s c) 
%       -> rotates fiels/vectors/object cw. for a_f > 0
%
%   (2) R_c = (c  s;-s c)
%       -> rotates coordinates system cw. for a_c > 0
%
%   General:
%       Rotation of coordinate system by an angle a= a_c is
%       equivalent with rotation of the fields by a= -a_f.
%
%   ==> Transfer functions IV,Z rotate like objects! <==
%
%   Rotation of coordinate system by a= a_c (-> a_f= -a_c):
%               -----------------
%   (1) use R=R_f & a=a_f =-a_c
%       -> IV'=R*IV     Z'=R*Z*RT
%
%   (2) use R=R_c & a=a_c =-a_f
%       -> IV'=R*IV     Z'=R*Z*RT
%
%   Note: the classical rotation matrix is defined as R=R_f
%-----------------------------------------------------------
%
% Here:
%   - strike: a_s. 
%     To get from strike direction into 'measured' system: 
%     + fields/TFs are to be rotated by a_f= a_s
%     + coordinate system is to be rotated by a_c= -a_s
%
%   - used: R=R_f and a=a_s
%
%     (program 'strike' uses R=R_c & a=a_c=-a_s
%                                 -> subroutine: estim_imp)
% ----------------------------------------------------------

% definitions

name = 'tst000'
                   
% reference resistivity to rescale z in order to keep rho ~ constant
rimp = 1000;

% 2-D tensor - Zxy & Zyx
z2d  = [.6+i*.3, -.2-i*.1]

% 'regional' strike
beta = 10

% distortion matrix
%dist = [1  .2; -.1  .8]
dist = [1 .3;.3  1]

% periods
per  = logspace(0,4,32);
nper = size(per,2);

% mu_0 * omega
mu0o = 4*pi*10^-7 * 2*pi*per.^-1;

% error floor for random noise (relative error)
ef = .05

c    = cos(beta*pi/180);
s    = sin(beta*pi/180);
r    = [c, -s; s, c]
% -> this is a regular 2-D rotation matrix

% 'impedance' tensors

% make a nice random variation of z2d w/ period
rndcos = rand(4,1)*ones(1,nper).*cos(rand(4,1)*2*pi*[1:nper]/nper) + 1;

for iper=1:nper
    % random sinusoidal variation
    z2 = [0, rndcos(1,iper)*real(z2d(1))+i*rndcos(2,iper)*imag(z2d(1));...
             rndcos(3,iper)*real(z2d(2))+i*rndcos(4,iper)*imag(z2d(2)), 0];
    
    % distortion
    z  = dist*z2;
    % rotation
    zr = r*z*r';
    
    % add relative error to 'measured' tensor
    zr = (ones(2,2)+(rand(2,2)-1)*ef).*zr;
    
    zr = reshape(zr.',1,4);
    
    zrot(iper,:) = zr;
end

% rho = |z|^2 / mu0o + impose random variation
zrot   = sqrt(rimp)*(mu0o').^.5*ones(1,4).*zrot;

rhorot = mu0o'.^-1 * ones(1,4) .* (abs(zrot).^2);
phirot = atan2(imag(zrot),real(zrot))*180/pi;

% write a *.dat file
fid = fopen([name,'.dat'],'w');

fprintf(fid,'%s\n',['#WRITTEN BY GEN_JFile.m: ',name,' on ',date]);
fprintf(fid,'%s\n','>AZIMUTH   =         0.0');
fprintf(fid,'%s\n','>LATITUDE  =   -120.0000');
fprintf(fid,'%s\n','>LONGITUDE =     50.0000');
fprintf(fid,'%s\n','>ELEVATION =         0.0');
fprintf(fid,'%s\n',[name,'    0']);


RZ  = ['RXY';'RYX';'RXX';'RYY'];
ord = [2,3,1,4];

for iel = 1:4
   fprintf(fid,'%s\n',RZ(iel,:));
   fprintf(fid,'%d\n',nper);
   for iper = 1:nper
       fprintf(fid,'%10.2f %12.4f %10.3f %s\n',per(iper),rhorot(iper,ord(iel)),phirot(iper,ord(iel)),...
                                '      00.00      00.00       00.0   00.0    1.00    1.00');
   end
end

ZZ  = ['ZXX';'ZXY';'ZYX';'ZYY'];
for iel = 1:4
   fprintf(fid,'%s\n',[ZZ(iel,:),' S.I. UNITS']);
   fprintf(fid,'%d\n',nper);
   for iper = 1:nper
       fprintf(fid,'%10.2f %12.8f %12.8f %s\n',per(iper),real(zrot(iper,iel)),imag(zrot(iper,iel)),...
                                 '     0.0       1.00');
   end
end

fclose(fid);