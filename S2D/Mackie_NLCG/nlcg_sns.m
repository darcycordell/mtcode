function nlcg_sns(outroot,header,outpath)
%=============================================================================
%
% -> basically a copy of NLCG_MDL!
%
%  WARNING: TOPOGRAPHY NOT IMPLEMENTED!!!!
% 
% INPUTS:
%  header is a title for the model plot
%  looks for files: *.mdl, *.rsp, and *.sns
%
%  If a file scales.m exists in working directory, then you can set values
%  for following:
%
%  scale   = [minY maxY minZ maxZ] (kilometers);
%  rholims = [minlog10rho maxlog10rho]) (Ohm-m)
%  phalims
%  VE (vertical exaggeration)
%
%  Otherwise these are set automatically and you can alter the values and
%  save the file for future use. Note that missing or bad data in this
%  file may cause the function to return to defaults or quit with an error.
%  
%  If a file pal.m exists in working directory, it will be executed to
%  load a custom colormap "cmp". Otherwise a modified jet is used.
%

choice=menu('vertical scale','linear','log','go back main');
if choice==1
    zaxis='linear';
elseif choice==2
    zaxis='log';
else 
    % close all
    return
end

version=6;
%header='NLCG Sensitivities';

if exist('scales.m','file')
  scales;
  ymin=scale(1);       ymax=scale(2);
  zmin=scale(3);       zmax=scale(4);
  rhomin=rholims(1);   rhomax=rholims(2);
  phamin=phalims(1);   phamax=phalims(2);
  % VE is also set in file scales.m!!
  set_params=0;
else
  set_params=1;
end

% load the sensitivity file
name=[outroot,'.sns']
if ~exist(name,'file')
    error(['File: ',name,' does not exist'])
end
[ny,nz,dy,dz,sens]=load_sns(name,version);
dy=dy/1000; dz=dz/1000; y0=0;  y0=y0/1000;

% load the site indices
name=[outroot,'.rsp'];
if exist(name,'file')
    [indsit,sitrms]=load_sites(name);
else
    error(['File ',name,' does not exist.'])
end

nsites=max(size(indsit));

% compute y's and z's at block edges
y=zeros(1,ny+1);
for k=1:ny  y(k+1)=y(k)+dy(k);   end
z=zeros(1,nz+1);
for k=1:nz  z(k+1)=z(k)+dz(k);   end

% compute y's at sites
ysite=y(indsit)+dy(indsit)'/2;

% offset site(1) to y0
y=y-ysite(1)+y0;
ysite=ysite-ysite(1)+y0;
%comment out below lines to plot stations in the middle of each cell
%ET,Jan,2008
if exist([outroot,'.stn'])
    stnfile=parfile; lstn=length(stnfile);
    stnfile(lstn-2:lstn)='stn';
    Aa=textread(stnfile,'%s'); Aa(1:2)=[]; Aa=char(Aa); nstn=size(Aa); nstn=nstn(1);
    for oi=1:nstn
        cf=Aa(oi,:);
        cg=find(cf==';');
        dd(oi,:)=str2num(cf(cg+1:length(cf)));
    end
    dd=cumsum(diff(dd));dd=[0.0001;dd];
    ysite;
    %plot(ysite,'o');hold on;
    ysite=dd';
end
%comment out above lines to plot stations in the middle of each cell
%plot(ysite,'*')


% defaults
if set_params
     ymax=max(ysite);        ymin=min(ysite); dely=ymax-ymin;
     ymin=ymin-1*dely;       ymax=ymax+1*dely;
%     zmin=min(z);            zmax=max(z);
     zmin=min(z);            zmax=10;
     % vertical exageration always starts at 1:1
     VE=1;
end
sensmin=-4; sensmax=1;

% adjust colorbar toggle
adjclb = 0;
colbw = 0;
%choice=menu('plot type','flat','interpolated');
% pmv5 : hardwired
choice=1;
if choice==1  typ='flat';   else  typ='interp';   end

% extend model so that pcolor will plot all elements
sens(nz+1,:)=sens(nz,:);   sens(:,ny+1)=sens(:,ny);

% plot the model
done=0;
h1 = figure;
while ~done
    figure(h1);
    clf
    scale=[ymin ymax zmin zmax];
    senslims=[sensmin sensmax];
  
    % specify color map
    % colmap=menu('color map','John Color','Matlab Jet Color');
    colmap=2;
    if colbw==0
        if colmap==1;              
            color=[  0    0    0    0  154  255  255  255  255;...
                     0  191  255  255  205  255  220  165    0;...
                   255  255  255    0   50    0    0    0    0]';
            cmap=color/255;
        else
            cmap=jet(20);
            cmap(20,:) =[];  cmap(19,:)=[];  cmap(13,:)=[]; 
            cmap(7,:)  =[];  cmap(1,:) =[];
            % below: Volcan's favourite
            % cmap=jet(24);
            % cmap(23,:)=[];    cmap(21,:)=[];    cmap(19,:)=[];    cmap(17,:)=[];
            % cmap(15,:)=[];    cmap(13,:)=[];    cmap(11,:)=[];    cmap(9,:)=[];
            %cmap(7,:)=[];     cmap(5,:)=[];     cmap(3,:)=[];
            cmap=jet(24);
            cmap(23,:)=[];    cmap(21,:)=[];    cmap(19,:)=[];    cmap(17,:)=[];
            cmap(15,:)=[];    cmap(13,:)=[];    cmap(11,:)=[];    cmap(9,:)=[];
            cmap(7,:)=[];     cmap(5,:)=[];     cmap(3,:)=[];
        end
    else
        cmap=gray(13);
    end 

    pcolor(y,z,sens);
  
    hold on
    caxis(senslims);
    colormap(cmap);
    if strcmp(zaxis,'log')
        plot(ysite,1.0000001*z(2)*ones(nsites,1),'kv','markerfacecolor','k','markersize',5)
    else
        plot(ysite,z(1)*ones(nsites,1),'kv','markerfacecolor','k','markersize',5)
    end
    title(header,'Interpreter','none');
    axis(scale);
    xlabel('km');ylabel('km');
    %set(gca,'YTick',[0:50:350])
    hold off
  
    if strcmp(zaxis,'log'); axis square; else; set(gca,'DataAspectRatio',[VE 1 1]); end

    if strcmp(typ,'flat'); shading flat; else; shading interp;  end

    if strcmp(zaxis,'log')
        set(gca,'YScale','log','YDir','reverse')
    else
        set(gca,'YScale','linear','YDir','reverse')
    end

    hcb = colorbar('ver');
    pltpos = get(gca,'Position');
  
    if adjclb==1
        pltpos = get(gca,'Position');
        pratio = get(gca,'PlotBoxAspectRatio');
        pratio = pratio(2)/pratio(1)*1.25;
        clbpos = get(hcb,'Position');
        clbpos(2) = clbpos(2)+(1-pratio)/2*clbpos(4);
        clbpos(4) = pratio*clbpos(4);
        if pratio<1
            set(hcb,'Position',[clbpos(1),clbpos(2),clbpos(3),clbpos(4)]);
            % needed for Matlab7: reset the model plot to old position
            set(gca,'Position',pltpos);            
        end
    end
    if strcmp(zaxis,'linear')
        choice=menu('Change:','sensitivity limits','y limits','z limits',...
                    'vertical exageration','adjust colorbar','color/bw','nothing');
    else
        choice=menu('Change:','sensitivity limits','y limits',...
                    'z limits','nothing');
    end
    if choice==1
        set_params=1;
        sensmin=input('Enter min log10(sensitivity): ');
        while 1
            sensmax=input('Enter max log10(sensitivity) > min: ');
            if sensmax>sensmin break; end
        end
    elseif choice==2
        set_params=1;
        ymin=input('Enter left edge (km): ');
        while 1
            ymax=input('Enter right edge (> left edge) (km): ');
            if ymax>ymin break; end
        end
    elseif choice==3
        set_params=1;
        while 1
            zmin= input('Enter top edge (>=0) (km): ');
            if zmin>=0 break; end
        end
        while 1
            zmax= input('Enter bottom edge (> top edge) (km): ');
            if zmax>zmin break; end
        end
    elseif choice==4
        if strcmp(zaxis,'linear')
            kk=menu('Vertical Exageration:','1:1','x2','x3','x10');
            if kk==1;       VE=1;
            elseif kk==2;   VE=2*VE;
            elseif kk==3;   VE=3*VE;
            elseif kk==4;   VE=10*VE;
            end    
        else
            done=1;
        end 
    elseif choice==5
        % adjust colorbar toggle
        adjclb=rem(adjclb+1,2);
    elseif choice==6
        % change: color<->black and white
        colbw=rem(colbw+1,2);
    else
        done=1;
    end
end

scolbw={'color','bw'};
%choice=menu('Make Postscript File?','yes','no');
if strcmp(zaxis,'log') 
    sensstr='sens_LZ_'; 
else
    sensstr='sens_Z_';
end
choice =1;
if choice==1
    eval(['print -dpsc ',outpath,'/',sensstr,char(scolbw(colbw+1)),'.ps']);
end

zp=menu('Make sens(z)plot?','yes','no');
if zp==1
    choice=menu('All sites?','yes','no');
    while 1
        if choice == 1;
            ismin=1; ismax=nsites;
        else
            ns=input('Enter station number  (0 to end): ');
            ismin=ns; ismax=ns;
            if ns==0 break; end
        end
        for is=ismin:ismax
            h2 = figure;
            sensz = sens(:,indsit(is))';
            plot(z,sensz,'o-'); axis([0,100,0,3]); title([' Site: ',num2str(is)]);
            xlabel('depth [km]');
            ylabel('log10(sensitivity)');
            
            % Output   z and log10 (rho(z)) to file
            str=[outpath,'/','sensz_',num2str(is),'.dat'];
            fid1=fopen(str,'w+');
            fprintf(fid1,'%9.3f %9.3f \n',[z;sensz]);
            fclose(fid1);
            
            str=[outpath,'/','sensz_',num2str(is),'.ps'];
            if choice==2 print(h2,'-dpsc',str); end
            
            pause
        end
        if choice==1 break; end
    end
end

% close all
end