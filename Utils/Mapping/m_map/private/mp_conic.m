function [X,Y,vals,labI]=mp_conic(optn,varargin);
% MP_CONIC  Conic projections
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and their inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
% These are  conic projections with two standard parallels, useful
% for showing limited areas at mid-latitudes.
%    Albers equal-area - has an equal-area property
%    Lambert conformal - is conformal
%
%  7/6/99 - fixed tendency to re-define .ulongs if .clong set by user
%  3/4/02 - added error if parallels are equidistant from equator (i.e. not conic projection really)

global MAP_PROJECTION MAP_VAR_LIST

name={'Albers Equal-Area Conic','Lambert Conformal Conic'};

pi180=pi/180;

switch optn,

  case 'name',

     X=name;

  case {'usage','set'}

     X=char({['     ''' varargin{1} ''''],...
              '     <,''lon<gitude>'',[min max]>',...
              '     <,''lat<itude>'',[min max]>',...
              '     <,''clo<ngitude>'',value>',...
              '     <,''par<allels>'',[lat1 lat2]>',...
              '     <,''rec<tbox>'', ( ''on'' | ''off'' )>'});

  case 'get',

     X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' MAP_PROJECTION.routine ')'],...
            [' longitudes: ' num2str(MAP_VAR_LIST.ulongs) ' (centered at ' num2str(MAP_VAR_LIST.clong) ')'],...
            [' latitudes: ' num2str(MAP_VAR_LIST.ulats) ],...
            [' standard parallels: ' num2str(MAP_VAR_LIST.parallels) ],...
            [' Rectangular border: ' MAP_VAR_LIST.rectbox ]);

  case 'initialize',

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    MAP_VAR_LIST.ulongs=[-180 -50];
    MAP_VAR_LIST.ulats=[10 85];
    MAP_VAR_LIST.parallels=NaN;
    MAP_VAR_LIST.clong=NaN;
    MAP_VAR_LIST.rectbox='off';
    k=2;longs_def=0;
    while k<length(varargin),  
      switch varargin{k}(1:3),
         case 'lon', 
           MAP_VAR_LIST.ulongs=varargin{k+1};longs_def=1;
           if MAP_VAR_LIST.ulongs(1)>MAP_VAR_LIST.ulongs(2),
             MAP_VAR_LIST.ulongs=MAP_VAR_LIST.ulongs([2 1]);
           end;
         case 'clo',
           MAP_VAR_LIST.clong=varargin{k+1};
         case 'lat',
           MAP_VAR_LIST.ulats=varargin{k+1};
         case 'par',
           MAP_VAR_LIST.parallels=varargin{k+1};
         case 'rec',
           MAP_VAR_LIST.rectbox=varargin{k+1};
         otherwise
           disp(['Unknown option: ' varargin{k}]);
         end;
       k=k+2;
     end;
    if isnan(MAP_VAR_LIST.clong), MAP_VAR_LIST.clong=mean(MAP_VAR_LIST.ulongs); 
    elseif ~longs_def, MAP_VAR_LIST.ulongs=MAP_VAR_LIST.clong+[-180 180];  end;
    if isnan(MAP_VAR_LIST.parallels), MAP_VAR_LIST.parallels=mean(MAP_VAR_LIST.ulats)*[1 1]; end;

    MAP_VAR_LIST.rlongs=MAP_VAR_LIST.ulongs*pi180;
    MAP_VAR_LIST.rlats=MAP_VAR_LIST.ulats*pi180;
    MAP_VAR_LIST.rparallels=MAP_VAR_LIST.parallels*pi180;

    % These are constants used by the projection formulas

    switch MAP_PROJECTION.name
      case name(1),   
        MAP_VAR_LIST.n=sum(sin(MAP_VAR_LIST.rparallels))/2;
	if MAP_VAR_LIST.n==0, error('Your parallels are equidistant from the equator - use a cylindrical projection!'); end;
        MAP_VAR_LIST.C=cos(MAP_VAR_LIST.rparallels(1)).^2+2*MAP_VAR_LIST.n*sin(MAP_VAR_LIST.rparallels(1));
        MAP_VAR_LIST.rho0=sqrt(MAP_VAR_LIST.C-2*MAP_VAR_LIST.n*sin(mean(MAP_VAR_LIST.rlats)))/MAP_VAR_LIST.n;
      case name(2),
        if diff(MAP_VAR_LIST.parallels)==0,
          MAP_VAR_LIST.n=sin(MAP_VAR_LIST.rparallels(1));
        else
          MAP_VAR_LIST.n=-diff(log(cos(MAP_VAR_LIST.rparallels)))/diff(log(tan(MAP_VAR_LIST.rparallels/2+pi/4)));
        end;   
        MAP_VAR_LIST.F=cos(MAP_VAR_LIST.rparallels(1))/MAP_VAR_LIST.n* ...
                       tan(pi/4+MAP_VAR_LIST.rparallels(1)/2).^MAP_VAR_LIST.n;
        MAP_VAR_LIST.rho0=MAP_VAR_LIST.F/tan(pi/4+mean(MAP_VAR_LIST.rlats)/2).^MAP_VAR_LIST.n;
    end;

    % Get X/Y and (if we are in a box) update the lat/long limits.

    mu_util('xylimits');
    if strcmp(MAP_VAR_LIST.rectbox,'on'),  mu_util('lllimits'); end;


  case 'll2xy',

    long=varargin{1};
    lat=varargin{2};
    vals=zeros(size(long));
    
    % Clip out-of-range values (lat/long box)
    
    if ~strcmp(MAP_VAR_LIST.rectbox,'on') & ~strcmp(varargin{4},'off'),
        vals=vals | long<=MAP_VAR_LIST.longs(1)+eps*10 | long>=MAP_VAR_LIST.longs(2)-eps*10 | ...
	              lat<=MAP_VAR_LIST.lats(1)+eps*10 |   lat>=MAP_VAR_LIST.lats(2)-eps*10;
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(1),long<MAP_VAR_LIST.longs(1),lat);
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(2),long>MAP_VAR_LIST.longs(2),lat);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(1),lat<MAP_VAR_LIST.lats(1),long);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(2),lat>MAP_VAR_LIST.lats(2),long);
    end;

    switch MAP_PROJECTION.name
      case name(1),   
        rho=sqrt(MAP_VAR_LIST.C-2*MAP_VAR_LIST.n*sin(lat*pi180))/MAP_VAR_LIST.n;
      case name(2),
        lat(lat==-90)=-89.999; % Prevents /0 problems in next line
        rho=MAP_VAR_LIST.F ./ tan(pi/4+lat*pi180/2).^MAP_VAR_LIST.n;
    end;
    theta=MAP_VAR_LIST.n*(long-MAP_VAR_LIST.clong)*pi180;

    X=real(rho.*sin(theta));    
    Y=real(MAP_VAR_LIST.rho0-rho.*cos(theta));

    % Clip out-of-range values (rectangular box)

    if strcmp(MAP_VAR_LIST.rectbox,'on') & ~strcmp(varargin{4},'off'),
        vals= vals | X<=MAP_VAR_LIST.xlims(1)+eps*10 | X>=MAP_VAR_LIST.xlims(2)-eps*10 | ...
                     Y<=MAP_VAR_LIST.ylims(1)+eps*10 | Y>=MAP_VAR_LIST.ylims(2)-eps*10;
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(1),X<MAP_VAR_LIST.xlims(1),Y);
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(2),X>MAP_VAR_LIST.xlims(2),Y);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(1),Y<MAP_VAR_LIST.ylims(1),X);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(2),Y>MAP_VAR_LIST.ylims(2),X);
    end;


  case 'xy2ll',

    pi180=pi/180; 
    
    switch MAP_PROJECTION.name
      case name(1),   
        rho=sqrt(varargin{1}.^2+(MAP_VAR_LIST.rho0-varargin{2}).^2);
        theta=atan(varargin{1}./(MAP_VAR_LIST.rho0-varargin{2}));
        Y=asin((MAP_VAR_LIST.C-(rho*MAP_VAR_LIST.n).^2)/(2*MAP_VAR_LIST.n))/pi180;

      case name(2),
        rho=sign(MAP_VAR_LIST.n)*sqrt(varargin{1}.^2+(MAP_VAR_LIST.rho0-varargin{2}).^2);
        theta=atan(varargin{1}./(MAP_VAR_LIST.rho0-varargin{2}));
        Y=(2*atan((MAP_VAR_LIST.F./rho).^(1/MAP_VAR_LIST.n))-pi/2)/pi180;
    end;
    % Clip out-of-range values (lat/long box)
    X=MAP_VAR_LIST.clong+(theta/MAP_VAR_LIST.n)/pi180;

  case 'xgrid',
   
    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},3,varargin{2});

  case 'ygrid',
   
    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},31,varargin{2});

  case 'box',

    [X,Y]=mu_util('box',31);
        
end;



