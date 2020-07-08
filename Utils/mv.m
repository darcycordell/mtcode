function v = mv(varargin)
% Function which takes any number of matrices of different shapes
% and converts it all into a vector (mv = "make vector").
% I don't think a similar function exists in MATLAB yet (if there is one, 
% we should use the built in version).
%
% Usage: v = mv(varargin)
%
% Inputs: Any number of matrices of any shape
%
% Outputs: A vector of all the matrices' entries in a single column
%
% Example: If Matrix A is N x M and Matrix B is P x Q, then this function
% will output a vector v which has length (N x M x B x Q) with entries of A
% and B.
%

v=[];
for i = 1:length(varargin)
    
    tmp = varargin{i};
    
    v = [v; tmp(:)];
    
end

end


    