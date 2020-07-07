function [F,ra,ph]=mt_fwd_occ(pr,nl,w,thick)
%This function computes the MT response of HLE for Occam's inversion
%Ersan Turkoglu, 2004------eturk@phys.ualberta.ca
%f=frequency, pr=[res1...resn,t1...tn-1(meters)], ra=computed apparent resistivity, ph=computed phase.
%# of layers angular frequency permittivity and iwm are global variables
%nl=number of layers
%w=2*pi*f;                                 %Angular frequency
mua=4*pi*10^-7;                           %Magnetic permeability
iwm=1i*mua*w; 
pr=([1./pr(1:nl),thick]);
Q(nl,:)=ones(1,length(w));
for n=nl-1:-1:1 
    Q(n,:)=(Q(n+1,:)+sqrt(pr(n+1)/pr(n))*tanh(sqrt(iwm*pr(n))*pr(n+nl)))./(sqrt(pr(n+1)/pr(n))+Q(n+1,:).*tanh(sqrt(iwm*pr(n))*pr(n+nl)));
end
Z=sqrt(iwm/pr(1)).*Q(1,:);                %Empedance
ra=1./(w*mua).*abs(Z).^2;                 %Apparent Resistivity
ph=atan(imag(Z(1,:))./real(Z(1,:)))*180/pi;        %phase
% [f',ra',ph']

F=[real(Z)';imag(Z)']';


