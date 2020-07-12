function [x, resnorm,outeriter] = nnls(L,d)
% Function to compute the non-negative least squares algorithm. This is a
% modification of a normal least squares solution with the additional
% constraint that the solution cannot be negative
%
% Usage: [x, resnorm, outeriter] = nnls(L,d)
%
% Inputs:
%       L = operator matrix (N x N)
%       d = data vector (N x 1)
%
% Outputs:
%       x = solution vector (N x 1)
%       resnorm = the norm of the residuals (data misfit)
%       outeriter = number of iterations
%
%

tol = 10*eps*norm(L,1)*length(L);%tolerance for the stopping criterion
n = size(L,2);
% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
wz = nZeros;

% Initialize set false logicals
P = false(n,1);
% Initialize set of true logicals
Z = true(n,1);
%Initialize data vectors of zeros
x = nZeros;

%difference between data and modelled data
diff = d - L*x;
%initialize w
w = L'*diff;

% Set up iteration criterion
outeriter = 0;
iter = 0;
itmax = 3*n; %maximum iterations

% Outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > tol)
   outeriter = outeriter + 1;
   % Reset intermediate solution z
   z = nZeros; 
   % Create wz, a Lagrange multiplier vector of variables in the zero set.
   % wz must have the same size as w to preserve the correct indices, so
   % set multipliers to -Inf for variables outside of the zero set.
   wz(P) = -Inf;
   wz(Z) = w(Z);
   % Find variable with largest Lagrange multiplier
   [~,t] = max(wz);
   % Move variable t from zero set to positive set
   P(t) = true;
   Z(t) = false;
   % Compute intermediate solution using only variables in positive set
   z(P) = L(:,P)\d;
   % inner loop to remove elements from the positive set which no longer belong
   while any(z(P) <= 0)
       iter = iter + 1;
       if iter > itmax
           exitflag = 0;
           resnorm = sum(diff.*diff);
           x = z;
           lambda = w;
           return
       end
       % Find indices where intermediate solution z is approximately negative
       Q = (z <= 0) & P;
       % Choose new x subject to keeping new x nonnegative
       alpha = min(x(Q)./(x(Q) - z(Q)));
       x = x + alpha*(z - x);
       % Reset Z and P given intermediate values of x
       Z = ((abs(x) < tol) & P) | Z;
       P = ~Z;
       z = nZeros;           % Reset z
       z(P) = L(:,P)\d;      % Re-solve for z
   end
   x = z;
   diff = d - L*x;
   w = L'*diff;
end

lambda = w;
resnorm = diff'*diff;
output.iterations = outeriter;
end %END NNLS