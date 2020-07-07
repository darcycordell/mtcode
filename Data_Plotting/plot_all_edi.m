function plot_all_edi(varargin)

% Function to plot and save ALL edi files in current directory in Portrait direction
% diagonal components use circle symbol
% off-diagonal square
% Error bars are plotted 
% Adapted by Darcy Cordell - June 2017
% Adapted from plot_edi_2.m written by Greg Nieuwenhuis - June 2013
% 

%  USAGES:

% plot_all_edi; - function reads in all edi files and plots the data

% plot_all_edi(period_limits, rhoa_limits, tipper_limits, phase limits, zrot, trot); - function
%               reads all edi files, plot limits and rotations are
%               specified
%
% Example: plot_all_edi([0.01 1000],[1 1000],[0 1],[0 90],0,0)

% Note that the function assumes the tipper and resistivity is rotated the
% same as the impedance (no rotation is done here)
%-------------------------------------------------------------------------


warning off;
mkdir plot_all_edi

%Get List of All EDI files in Current Directory -------------------------
d=dir;
fff=length(d);
count=1;
for iio=1:fff
    [pathstr, name, ext] = fileparts(d(iio).name);
    if strcmp(ext,'.edi')
        edst{count}=d(iio).name;
        count=count+1;
    else
        disp([d(iio).name,' is not an edi file']);
    end
    
end
nsta=length(edst); %nsta is the number of edi files in the directory, fff is the maximum number of characters in the edi filename

%------------------------------------------------------------------------


%Loop over all stations
for iyu = 1:nsta
    
%Scale factor for SI units
%cf = (4*pi*1e-4)*(4*pi*1e-4);
cf =1;

if nargin== 0 % Just read all EDIs with pre-defined plot limits. No rotations.
    clf
    [z1,zvar1,hz1,hzvar1,rhoa1,rhoaerr1,phs1,phierr1,f1,rot1,coords1] = read_edi_imp(edst{iyu});
    
    per_lim=[.03 100000]; %limits for plotting period  
    rhoa_lim=[.1 100000];
    tip_lim=[-1 1];
    phs_lim=[-5,95];
    freq_lim=sort(1./per_lim);
elseif nargin == 6 % Plot and rotation parameters. Will rotate impedance / tipper

    per_lim=varargin{1};
    rhoa_lim=varargin{2};
    tip_lim=varargin{3};
    phs_lim=varargin{4};
    zrot =varargin{5};
    trot= varargin{6};
    
    freq_lim = sort(1./per_lim);
    
    clf
    [z1,zvar1,hz1,hzvar1,r1,rerr1,phs1,p1,f1,rot1,coords1] = read_edi_imp(edst{iyu});
     
    [z1,zvar1] = rot_z_dr(z1,zvar1,length(f1),rot1(1)); % Rotate impedance
    [rhoa1,rhoaerr1,phs1,phierr1]=calc_MT(z1/796,zvar1/796,f1); % Calculate app res / phase
 
    [hz1,hzvar1] = rot_tip(hz1,hzvar1,trot(1));         % Rotate tipper

    
else
    error('Either must have 0 input arguments or 6 input arguments. See code for details')
end

period1=1./f1; 

if length(rot1)> 1;     rot1=rot1(1);  end


[path1, station_name1, extension1] = fileparts(edst{iyu});



%--------------------------------------------------------------------------
subplot(13,1,[1:7])

% go through and make sure errobars never become negative (doesn't work on log-log plot)
for iper=1:length(f1)
    if (rhoa1(1,1,iper)-rhoaerr1(1,1,iper)) <=0; rhoaLerr1(1,1,iper)=rhoa1(1,1,iper)-.00000001; else rhoaLerr1(1,1,iper)=rhoaerr1(1,1,iper)/2;  end
    if (rhoa1(2,2,iper)-rhoaerr1(2,2,iper)) <=0; rhoaLerr1(2,2,iper)=rhoa1(2,2,iper)-.00000001; else rhoaLerr1(2,2,iper)=rhoaerr1(2,2,iper)/2;  end
    if (rhoa1(2,1,iper)-rhoaerr1(2,1,iper)) <=0; rhoaLerr1(2,1,iper)=rhoa1(2,1,iper)-.00000001; else rhoaLerr1(2,1,iper)=rhoaerr1(2,1,iper)/2;  end
    if (rhoa1(1,2,iper)-rhoaerr1(1,2,iper)) <=0; rhoaLerr1(1,2,iper)=rhoa1(1,2,iper)-.00000001; else rhoaLerr1(1,2,iper)=rhoaerr1(1,2,iper)/2;  end
end


%errorbarloglog(f1,squeeze(rhoa1(1,1,:)),squeeze(rhoaLerr1(1,1,:)),squeeze(rhoaerr1(1,1,:)./2),'o',[1 .6 .8],0.5,[1 .6 .8]); hold on;
%errorbarloglog(f1,squeeze(rhoa1(2,2,:)),squeeze(rhoaLerr1(2,2,:)),squeeze(rhoaerr1(2,2,:)./2),'o',[0 .6 1],0.5,[0 .6 1]);
%errorbarloglog(f1,squeeze(rhoa1(2,1,:)),squeeze(rhoaLerr1(2,1,:)),squeeze(rhoaerr1(2,1,:)./2),'s',[0 0 1],0.5,[0 0 1]); hold on;
%errorbarloglog(f1,squeeze(rhoa1(1,2,:)),squeeze(rhoaLerr1(1,2,:)),squeeze(rhoaerr1(1,2,:)./2),'s',[1 0 0],0.5,[1 0 0]);

%loglog(f1, squeeze(cf*rhoa1(1,1,:)),'ro'); hold on;
%loglog(f1, squeeze(cf*rhoa1(2,2,:)),'bo');
loglog(f1, squeeze(cf*rhoa1(1,2,:)),'rs'); hold on;
loglog(f1, squeeze(cf*rhoa1(2,1,:)),'bs');

set(gca,'xdir','reverse')
set(gcf,'Unit','centimeters','position',[1 0.2 12 20])
set(gca,'xscale','log','yscale','log')
set(gca,'xticklabel','')
axis([freq_lim rhoa_lim])
ax=axis;
ylabel('Apparent Resistivity (Ohm-m)')


%tstr = [station_name1,'(Red-Blue)',' and ',station_name2, '(Green-black) \theta =',num2str(rot1,'%.0f'),' ']
tstr = [station_name1,': xy = red; yx = blue : ', ' \theta =',num2str(rot1,'%.0f'),'\circ'];
    

title(strrep(tstr,'_','\_'));   % Force MATLAB to plot underscores

set(gca,'FontSize',8)
grid on
grid minor



 
%--------------------------------------------------------------------------
subplot(13,1,[8:9])
% errorbar(f1,phs1(1,1,:),phierr1(1,1,:),'o','MarkerSize',3,'color',[1 .6 .8]);hold on;
if mean(phs1(1,2,:)) < -100;    phs1(1,2,:)=phs1(1,2,:)+180;   end
if mean(phs1(2,1,:)) < -100;    phs1(2,1,:)=phs1(2,1,:)+180;   end

%errorbar(f1,phs1(1,2,:),phierr1(1,2,:),'r','MarkerSize',1,'Marker','none'); hold on;
%errorbar(f1,phs1(2,1,:),phierr1(2,1,:),'b','MarkerSize',1,'Marker','none'); 

semilogx(f1,squeeze(phs1(1,2,:)),'rs'); hold on;
semilogx(f1,squeeze(phs1(2,1,:)),'bs');

set(gca,'xdir','reverse')
set(gca,'xscale','log','FontSize',8)
set(gca,'xticklabel','')
axis([freq_lim phs_lim])
ylabel('Phase (degrees)')
grid on
grid minor
%--------------------------------------------------------------------------

subplot(13,1,[10:11])
%errorbar(f1,real(hz1(1,:)),real(hzvar1(1,:)),'r','MarkerSize',1,'Marker','none'); hold on;
%errorbar(f1,imag(hz1(1,:)),imag(hzvar1(1,:)),'b','MarkerSize',1,'Marker','none'); 

semilogx(f1,squeeze(real(hz1(1,:))),'rs'); hold on
semilogx(f1,squeeze(imag(hz1(1,:))),'bs');

set(gca,'xdir','reverse')
set(gca,'xscale','log','FontSize',8)
set(gca,'xticklabel','')
axis([freq_lim tip_lim])
% xlabel('Period (s)')
ylabel('Tzx')
grid on
grid minor
%--------------------------------------------------------------------------
subplot(13,1,[12:13])
%errorbar(f1,real(hz1(2,:)),real(hzvar1(2,:)),'r','MarkerSize',1,'Marker','none'); hold on;
%errorbar(f1,imag(hz1(2,:)),imag(hzvar1(2,:)),'b','MarkerSize',1,'Marker','none');

semilogx(f1,squeeze(real(hz1(2,:))),'rs'); hold on
semilogx(f1,squeeze(imag(hz1(2,:))),'bs');

set(gca,'xscale','log','FontSize',8)
axis([freq_lim tip_lim])
xlabel('Frequency(Hz)')
set(gca,'xdir','reverse')
ylabel('Tzy')
grid on
grid minor

set(gcf,'PaperUnits','centimeters','PaperPosition',[1 0.2 12 20],'PaperPositionMode', 'manual')
saveas(gcf,['plot_all_edi\',station_name1,'.png'])

end


end