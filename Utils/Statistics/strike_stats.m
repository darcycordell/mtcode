function [avg_strike, std, avg_r] = strike_stats(strike,J)

%Function which calculates the radial average, standard deviation and
%resultant for a given vector of strike angles.
%
% Usage: [avg_strike,std,avg_r] = strike_stats(strike,J)
%
%
% Inputs:
%       strike: A vector of strike angles in degrees
%       J is a factor to multiply all the strike angles to get rid of any
%       quadrant ambiguities. Since MT has a 90 degree ambiguity of strike, then
%       J=4 for (most) MT applications.
%
% Outputs:
%       avg_strike: The mean radial strike in degrees
%       std: The radial standard deviation in degrees
%       avg_r: The mean resultant in degrees
%
%These calculations of radial averages, resultant and standard deviation
%are taken from Probability and Mathematical Statistics by Birnbaum and
%Lukacs 1972.

%Read in strike vector
strike_hold=strike;

% Because of the 90 degree ambiguity of analysis, we must multiply all the
% strike angles by 4 so that they plot in the proper quadrant. For example,
% multiplication by 4 results in a 90 degree strike being equal to 0 degree
% strike.
strike=strike(find(isnan(strike)==0));
strike=strike*J;
strike=strike*(pi/180); %convert to radians
n=length(strike);

%Take the cosine and sine of the strike angles (the x and y projections)
c=cos(strike);
s=sin(strike);

% average the cos and sin projections
avg_c=sum(c)/n;
avg_s=sum(s)/n;

%The resultant is 0 < r < 1 and is a measure of scatter. As you add the
%vectors, if they all pointed the same way, you would have r=1. If you had
%all the vectors pointing different directions, then r=0.
avg_r=sqrt(avg_c^2+avg_s^2);

%The average radial strike
avg_strike=atan2(avg_s,avg_c)*(180/pi)/J;

%The average radial standard deviation in degrees.
std=(1/J)*sqrt((-2*log(avg_r)))*(180/pi);

end

