function [dplus_vars] = dplus(Z,Zerr,T,err_flr)
% The D+ algorithm takes data (in impedance) and fits
% a 1D model to the data. The 1D model consists of a finite number of
% infinitely-thin, infinitely conductive sheets with finite conductance.
%
% Usage: [dplus_vars] = dplus(Z,Zerr,T)
%
%   Inputs: 
%       Z = vector of complex impedance (Nx1)
%       Zerr = vector of complex errors (Nx1)
%       T = period (s) (1xN)
%       err_flr = impedance error floor percent (e.g. 5%). *OPTIONAL*
%           If err_flr is not supplied, then the function uses
%           user_defaults "u.dplus_percent"
%
%   Outputs:
%       dplus_vars structure contains periods used for fitting (T),
%       impedance D+ fit (Z_d), impedance errors used for fitting (Zerr),
%       admittance (C_mod), estimated apparent resistivity (rho_d),
%       esimated phase (phi_d), rms misfit (rms), spectral function
%       parameters (a_step, lam_f, a_0, a_f), D+ model (tau vs. z), and
%       depth to infinitely conducting base (b).
%
% Note: If your data include out of quadrant phases (OOQP) these are
% inconsistent with 1-D MT D+ assumptions and D+ cannot fit these data
% points. If all your data are OOQP on one (or both) off-diagonal modes,
% then it will throw a warning to you.
%
% The YX mode usually plots in -180:-90 quadrant. So if you are using D+ to
% fit YX mode you need to rotate your impedances to the 0:90 quadrant
% first. This can be done easily with Z = Z*exp(1i*pi)
%
% For more info on the theory of D+ see:
%
%   Parker, R. L. (1980). The inverse problem of electromagnetic induction: 
%       existence and construction of solutions based on incomplete data. 
%       Journal of Geophysical Research: Solid Earth, 85(80), 4421–4428.
%   Parker, R. L. (1994). Geophysical Inverse Theory. Princeton University Press.
%   Parker, R. L., & Whaler, K. (1981). Numerical methods for establishing 
%       solutions to the inverse problem of electromagnetic induction. 
%       Journal of Geophysical Research: Solid Earth, 86, 9574–9584.


u = user_defaults;

if exist('err_flr','var')
    u.dplus_percent = err_flr;
end

nanidx = find(isnan(real(Z)));

Z(nanidx) = [];
Zerr(nanidx) = [];
T(nanidx) = [];

for ifreq=1:length(T)
    if u.dplus_percent ~= 0
        if abs(real((Zerr(ifreq)))/abs(Z(ifreq))) < u.dplus_percent/100
            Zerr(ifreq)=(u.dplus_percent/100)*abs(Z(ifreq));
            Zerr(ifreq) = Zerr(ifreq) + 1i*Zerr(ifreq);
        end
    end
end

% Zerr=(u.dplus_percent/100)*abs(Z);
% Zerr = Zerr + 1i*Zerr;

%K is the number of log frequencies, N is the number of frequencies. The
%errors are scaled by a factor so that Gaussian noise recovers an rms of
%1.0 when comparing noisy to noise-free data. For large N, this factor
%approaches 1. 
K = log10(max(T))-log10(min(T));
N = length(T);
%(N-K)/N
%Zerr = Zerr*((N-K)/N);

mu0=4*pi*10^-7; f=1./T;
w = 2*pi*f;

C = (1./(1i*w*mu0)).*Z;
%delC = real((1./(1i*w*mu0)).*Zerr);
delC = (1./(w*mu0)).*imag(Zerr)+1i*(1./(w*mu0)).*real(Zerr); %Error propagation for delC: Error for real and imaginary components of Z are swapped

% %Convert input rho and phi into complex admittance:
% modC=sqrt((rho_a.*T)./(2*pi*mu0)); %From Parker 1994 definition of complex admittance
% C=modC.*exp(1i*(phi-90)*(pi/180)); %The negative 90 degree phase correction exists because the 
%                                    %natural position of phase is in the bottom right quadrant 
%                                    %(for a halfspace, phi = -pi/4; see GEOPH 424 Notes).
% 
% %Error propagation of C. Rules vary.
% %delmodC=0.5*sqrt(T./(2*pi*mu0.*squeeze(rho_a(t,tf(t),:)))).*squeeze(rho_err(t,tf(t),:)); %This is calculated from error propagation rules of modulus
% delC=0.5.*rho_err.*abs(C)./rho_a;%<----Equivalent

%Create vector of data with 2N points (real and imaginary in separate rows)
N=length(C);d=[];deld=[];
for k=1:N
%     deld(2*k-1)=delC(k); %error in d
%     deld(2*k)=delC(k); %error in d

    deld(2*k-1)=real(delC(k));
    deld(2*k)=imag(delC(k));
    
    d(2*k-1)=real(C(k))./deld(2*k-1); %real data 
    d(2*k)=imag(C(k))./deld(2*k); %imaginary data
end
d=d'; deld=deld';




%*************************************************************************
% THE INVERSION: The D+ inversion solution takes the equation d=L*a and
% solves for the model parameters a. d is the data (2N x 1), L is the operater
% matrix (2N x M+1) and a is the spectral function weights (M+1 x 1). The
% spectral function weights are the numerators of the partial fraction
% representation of the D+ solution where:
%
% C = a0 + sum(a_n / (lambda_n + i*w)).
%
% The inversion is solving for a's and lambda's.

% *****************Choosing the starting lambda vector:*******************

dlamb=0.01; %the spacing between lambdas
lambda_log=[floor(log10(min(f)))-8:dlamb:ceil(log10(max(f)))+8]; %As lambda -> infty and dlamb -> 0
            %the more accurate the inversion. I have taken max lambda as
            %10^8 times the shortest frequency
            %and the min lambda as 10^-8 times smaller than the lowest frequency
            %with a 0.01 log increment. The log increment is arbitrarily
            %chosen.
            
lambda=[0 10.^lambda_log];  M=length(lambda); w=2*pi./T;

%*************** Build L operator matrix which is 2N x M+1*****************

L=zeros(N*2,M+1); 
for k=1:2:2*N
    L(k,1)=1./deld(k); %the first column of L is 1s and 0s for real and imag parts (representing first spectral component, a0)
end
for j=2:M+1
    for k=1:N
        L(2*k-1,j)=lambda(j-1)./(lambda(j-1)^2+w(k)^2)./deld(2*k-1); %real components of spectral function are in odd rows
        L(2*k,j)=-w(k)./(lambda(j-1)^2+w(k)^2)./deld(2*k); %imaginary components of spectral function are in even rows
    end
end

[a,~] = nnls(L,d); %NNLS Inversion to solve for spectral function components, a

%***************************Remove doublets*******************************
% The solution, a, will tend to be mostly zeros except for a small number
% of non-zero "spikes" in the spectral function. Often, pairs of spikes
% will appear very close together, so it is useful to interpolate and
% average these spikes into one spike
a0=a(1);
a_tmp=a(2:end);
%Remove zeros in a and find corresponding lambda
indl = find(a_tmp~=0);  lam_d=[lambda];  a_d=a_tmp(indl); lam_d=lam_d(indl);

k=1; i=1; a_f=[]; lam_f=[]; 
while i<=length(indl)
    if i==length(indl)
        a_f(k)=a_d(i);
        lam_f(k)=lam_d(i);
        i=i+1;
    else
        if indl(i+1)-indl(i)<=1 %if two indices are adjacent in number then remove doublets and interpolate
            a_f(k)=(((a_d(i+1))+(a_d(i))));
            lam_f(k)=10^((log10(lam_d(i))+log10(lam_d(i+1)))/2);
            k=k+1;
            i=i+2;
        else %otherwise leave the modeled value unchanged
            a_f(k)=a_d(i);
            lam_f(k)=lam_d(i);
            k=k+1;
            i=i+1;
        end
    end
end


%*************************************************************************

% FORWARD CALCULATION
%Do the forward calculation of the spectral components as a sum of partial
%fractions. This gives you your predicted value to compare to the data.
tmp=[]; C_mod=[];

for k=1:N
    if isempty(a_f)
        C_mod(k)=a0;
    else
        for i=1:length(a_f)
            tmp(i)=a_f(i)./(lam_f(i)-1i*w(k));
        end
        C_mod(k)=a0+sum(tmp);
    end
end
%Split up real and imaginary components so the vectors are the same length
for k=1:N
        d_pred(2*k-1)=real(C_mod(k))./(deld(2*k-1));
        d_pred(2*k)=imag(C_mod(k))./(deld(2*k));
end

C_mod=C_mod';

[stats] = statistics([real(C); imag(C)],[real(C_mod); imag(C_mod)],[real(delC); imag(delC)]);
rms = stats.rms;
%misfit = 10*sqrt(sum((((d)-(d_pred'))./deld).^2));


Z_d = 1i*w.*mu0.*C_mod;

% [statsZ] = statistics([real(Z); imag(Z)],[real(Z_d); imag(Z_d)],[real(Zerr); imag(Zerr)]);
% statsZ.rms

%Calculate predicted rho and phase and compare to original data
rho_d=(2*pi*mu0)./(T).*abs(C_mod).^2;
phi_d=atan2(imag(Z_d),real(Z_d))*(180/pi);


% PLOTTING A(LAM) VS LAM STEP FUNCTION
% The spectral function has some useful properties (see Weidelt 2005). This
% plots the step function:
if length(a_f)==0
    a_step=a0;
else
    a_step=a_f(1);
    for na=2:length(a_f)
        a_step(na)=a_step(na-1)+a_f(na);
    end
end

if isempty(a_f)
    
    disp('Warning: Your data is not consistent with D+. Either XY or YX (or both) modes have all phases out of quadrant')
    z = NaN;
    tau = NaN;
    b = NaN;
    
else
    
    %**************************************************************************

    % THE D+ MODEL!
    % The spectral function (a vs. lambda) gives you the modelled data (C_mod)
    % by partial fraction summation. However, the spectral function does not
    % have any physical intuition. To get that, you need to convert from the
    % spectral function (a(lambda)) to the D+ model (tau(z)) where tau is the
    % conductance and z is the depth of the conductivity delta function.

    %This is an adaptation of Rutishauser's algorithm (1957). Parker (1980) and
    %Parker and Whaler (1981) both use Rutishauser's algorithm but I have
    %not been able to find a good resource that details how the algorithm
    %works. Weidelt (2005) published an updated algorithm which is more clearly
    %laid out, so I have used it here. It uses the properties of continued
    %fractions and orthogonal polynomials.
    mu0=4*pi*10^-7;
    s0=sum(a_f); B=[]; alpha=[]; tau=[];dz=[]; z=[]; z(1)=a0;
    numlam=length(lam_f);
    p=zeros(numlam+2,numlam);p_til=zeros(numlam+1,numlam); %p and p_til are the constants of orthogonal polynomials
                                       %like p(0) + p(1)*x + p(2)*x^2 + p(3)*x^3 ...
                                       %where x = i*w
    p(2,:)=ones(1,numlam).*(1/sqrt(s0));
    B(2)=0;
    k=1;

    %solve for the polynomial, p which is an N-1 x N matrix where N is the
    %number of lambda. Each column represents each lambda and each row is the
    %polynomial value associated with (i*w)^n. 
    i = 2; count = 1; dz(1) = 0;
    while 1
        alpha(i)=dot(lam_f.*p(i,:),a_f.*p(i,:));
        p_til(i+1,:)=(lam_f-alpha(i)).*p(i,:)-B(i)*p(i-1,:);
        B(i+1)=sqrt(dot(p_til(i+1,:),a_f.*p_til(i+1,:)));
        p(i+1,:)=p_til(i+1,:)./B(i+1);

        tau(i-1)=p(i,1)^2/mu0; %tau is conductance
        dz(i)=-1/(p(i,1)*p_til(i+1,1)); %dz is distance between spikes
        z(i)=z(i-1)+dz(i); %z is depth of spike

        if z(i)<0 || dz(i) > 10^7 || dz(i) < 0
            break
        else
            count = count+1;
        end

        i = i+1;

    end

    count = 1; ind_del = []; dz(end) = [];
    for i = 2:length(dz)
        if dz(i) < 0.01*dz(i-1)
            dz(i-1) = dz(i-1) + dz(i);

            tau(i) = tau(i-1)+tau(i);

            ind_del(count) = i;
            count = count+1;
        end
    end

    if ~isempty(ind_del)
        dz(ind_del) = [];
        tau(ind_del-1) = [];
    end

    if length(dz) > numlam
        dz(numlam+1:end) = [];
        tau(numlam+1:end) = [];
    end


    z = a0+cumsum(dz);

    b = a0+sum(a_f(lam_f>0)./lam_f(lam_f>0));

end

dplus_vars.T = T;
dplus_vars.Z = Z_d;
dplus_vars.Zerr = Zerr;
dplus_vars.C = C_mod;
dplus_vars.rho = rho_d;
dplus_vars.pha = phi_d;
dplus_vars.rms = rms;
dplus_vars.a_step = a_step;
dplus_vars.lambda = lam_f;
dplus_vars.z = z;
dplus_vars.tau = tau;
dplus_vars.a0 = a0;
dplus_vars.a_f = a_f;
dplus_vars.b = b;
dplus_vars.Zorig = Z;

end
