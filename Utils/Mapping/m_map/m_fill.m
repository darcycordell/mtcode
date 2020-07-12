function h = m_fill(varargin)

% M_FILL
% 
% Matlab's FILL function for m_map package.
%
% USAGE: M_FILL(LON,LAT,[OPTIONS]) 
%
% Author: Mashrab Kuvatov

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if nargin < 2;
  help m_fill
  return
end

[x,y] = m_ll2xy(varargin{1},varargin{2});
varargin = varargin(:);
s = size(varargin,1);
h = fill(x,y,varargin{3:s});

if nargout == 0
  clear h
end
