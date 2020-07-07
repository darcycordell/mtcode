function cik=nlcg_mdl_iso(outroot,header,outpath)
%=============================================================================
%
%  WARNING: TOPOGRAPHY NOT IMPLEMENTED!!!!
% 
% INPUTS:
%  header is a title for the model plot
%  looks for files: *.mdl and *.rsp
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

% I will never use log scale, ET, Jan 2008
%choice=menu('vertical scale','linear','log','go back main');

disp('IN : nlcg_mdl_iso')

global seed elev


elev=elev/1000;
choice=1;
if choice==1
    zaxis='linear';
elseif choice==2
    zaxis='log';
else 
    % close all
    return
end

version=6;
%header='NLCG Model';

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
% load the model file
name=[outroot,'.mdl'];

if exist(name,'file')
    [ny,nz,dy,dz,model]=load_mdl(name,version);
else
    error(['File: ',name,' does not exist'])
end


% load the another model file?
% I will never compare here, ET, Jan 2008
%lmod2=menu('compare w/ second model?','yes','no');
lmod2=2;
if lmod2==1
    choice=menu('same first model?','same','other');
    if choice==2
        [mod1fil,mod1path] = uigetfile('*.mdl;*.mod','pick first model');
        mod1file=[mod1path,'/',mod1fil];
        disp(['reading first model file: ',mod1file]);
        [ny,nz,dy,dz,model]=load_mdl(mod1file,version);
    end
    [mod2fil,mod2path] = uigetfile('*.mdl;*.mod','pick second model');
    mod2file=[mod2path,'/',mod2fil];
    disp(['reading second model file: ',mod2file]);
    [ny2,nz2,dy2,dz2,model2]=load_mdl(mod2file,version);
    % relative change
    model=model-model2;
    header=['log[m_1] - log[m_2]'];
    rhomin=-1; rhomax=+1;
end
    
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
z=z-max(elev); 
% compute y's at sites
ysite=y(indsit)+dy(indsit)'/2;

% offset site(1) to y0
y=y-ysite(1)+y0;
ysite=ysite-ysite(1)+y0;

%comment out below lines to plot stations in the middle of each cell
%ET,Jan,2008
if exist([seed,'.stn'])
    stnfile=[outpath,seed,'.stn']; 
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
if set_params==0
     if ymin==0 & ymax ==0
        ymax=max(ysite);        ymin=min(ysite); dely=ymax-ymin;
     end
end
% defaults
if set_params
     ymax=max(ysite);        ymin=min(ysite); dely=ymax-ymin;
     ymin=ymin-0.1*dely;       ymax=ymax+0.1*dely;
%     zmin=min(z);            zmax=max(z);
     zmin=min(z);            zmax=10;
     rhomin=0;               rhomax=3;
     phamin=0;               phamax=90;   phalims=[phamin,phamax];
     % vertical exageration always starts at 1:1
     VE=1;
     % different limits for display of relative differences
     if lmod2==1 
         rhomin=-1; rhomax=+1;
     end
end

% adjust colorbar toggle
adjclb = 0;
colbw = 0;
%choice=menu('plot type','flat','interpolated');
% pmv5 : hardwired
choice=1;
if choice==1  typ='flat';   else  typ='interp';   end

% extend model so that pcolor will plot all elements
model(nz+1,:)=model(nz,:);   model(:,ny+1)=model(:,ny);






% plot the model
done=0;
cik=0;
%h1 = figure('PaperPositionMode','auto');
    %good for high res studies
    max_p=180; %maximum profile lenght possible in km
    max_d=180; % maximum depth of the profile possible in km

if header(1)=='c'
    % Compute conductance Added pmv16
    for iy=1:ny
        conductance(iy)=0.0;
        for iz = 1:nz
            if z(iz) <= zfix
                conductance(iy)=conductance(iy)+ 1000.*dz(iz)/(10^model(iz,iy));
            end
        end
    end
    conductance(ny+1)= conductance(ny);
    % Output   z and log10 (rho(z)) to file
    str=[outpath,'/',seed,'_conductance_',num2str(zfix),'_km.dat'];
    fid1=fopen(str,'w+');
    fprintf(fid1,'%9.3f %9.3f \n',[y;conductance]);
    fclose(fid1);
    scale=[ymin ymax zmin zmax];    got=get(0,'screensize');
    figure('position',[170 110 got(3)*0.80 got(4)*0.80],'PaperPositionMode','auto');
    subplot('position',[0.1 0.2 0.8*abs(scale(2)-scale(1))/max_p abs(scale(4)-scale(3))*0.8/max_d*got(3)/got(4)]);
    semilogy(y,conductance)
    axis([ymin,ymax,1,100000]); 
    xlabel ('Distance(km)')
    ylabel ('Conductance (S)')
    title(header)
    %eval(['print -djpeg ',outpath,'/',seed,'_conductance_',num2str(zlim),'_km.jpeg']);
    %end of conductance.
    cik=0;
    return
end

while ~done | ~cik
    %figure(h1);
    %clf
    scale=[ymin ymax zmin zmax];
    %rholims=[rhomin-.15 rhomax+0.15];
    rholims=[rhomin rhomax];
    % specify color map

    if colbw==0
      cmap=jet(64);
    else
      cmap=gray(13);
    end
    
    got=get(0,'screensize');
    figure('position',[170 110 got(3)*0.80 got(4)*0.80],'PaperPositionMode','auto');
    clf

    font_n='arial';
    font_s=12;
    save scale.mat scale
    set(gca,'fontname',font_n)
    subplot('position',[0.1 0.2 0.8*abs(scale(2)-scale(1))/max_p abs(scale(4)-scale(3))*0.8/max_d*got(3)/got(4)]);
    model(find(model==10))=NaN; %Replace air values with NaN, pcolor plots NaN as white
    pcolor(y,z,model);
    %save 'model_etc.mat' y z model ysite indsit
    set(gca,'fontname',font_n,'FontSize',font_s,'linewidth',1,'tickdir','in');
    hold on
    caxis(rholims);
    colormap(flipud(cmap));
    
    if strcmp(zaxis,'log')
        plot(ysite,1.0000001*z(2)*ones(nsites,1),'kv','markerfacecolor','k','markersize',5)
    else
        %change elev-0 to elev-xx when you have topography
        plot(ysite,-elev-0.1,'kv','markersize',7,'markerfacecolor','k')
        %plot(ysite,-elev*ones(nsites,1),'kv','markerfacecolor','k','markersize',5)
    end
    set(gca,'Fontname',font_n,'fontsize',font_s);
    title(header,'Interpreter','none');
    %title('$\rho$ [$\Omega\cdot$m]','Interpreter','latex','FontSize',5);
    axis(scale);
    xlabel('Distance (km)');ylabel('Depth (km)');
    %set(gca,'XTick',[-200:50:800])
    hold off
  
    if strcmp(zaxis,'log'); axis square; else; set(gca,'DataAspectRatio',[VE 1 1]); end

    if strcmp(typ,'flat'); shading flat; else; shading interp;  end

    if strcmp(zaxis,'log')
        set(gca,'YScale','log','YDir','reverse')
    else
        set(gca,'YScale','linear','YDir','reverse')
    end

%    hcb = colorbar('ver')
    mancol=1;
    if mancol==1
    if lmod2~=1  
        hcb = colorbar('hor','position',[(0.1+0.8*abs(scale(2)-scale(1))/max_p/2-0.125) 0.25 0.25 0.025],... % Change here to put color bar in good place
            'XTick',[0 0.5,1,1.5,2,2.5,3],...
            'FontSize',12,'fontname',font_n,'fontsize',font_s,...
            'XTickLabel',{'1','3','10','30','100','300','1000'});
        set(get(hcb, 'Title'), 'String', 'Resistivity (\Omegam)','fontsize',12,'fontname',font_n);
      %text(100,100,'Resistivity (\Omegam');
    elseif lmod2==1
        hcb = colorbar('ver',...
          'YTick',[-1.0,-.5,0,.5,1.0],...
          'YTickLabel',{'-1.0','-0.5','0.0','0.5','1.0'},...
          'FontSize',12);
    end
    end
    print -depsc -painters 'model.eps'
    print -djpeg100 -painters 'model.jpg'
    
    
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
        choice=menu('Change:','resistivity limits','y limits','z limits',...
                    'vertical exageration','adjust colorbar','color/bw','export model','nothing','quit');
    else
        choice=menu('Change:','resistivity limits','y limits',...
                    'z limits','nothing');
    end
    if choice==1
        set_params=1;
        rhomin=input('Enter min log10(resistivity (Ohm-m)): ');
        while 1
            rhomax=input('Enter max log10(resistivity) > min: ');
            if rhomax>rhomin break; end
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
            zmin= input('Enter top edge (-ve if using elev) (km): ');
            break;
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
        
    elseif choice ==7    
        % Write text file of model in geographic co-ords  %MJU 2016-04-09
        export_nlcg_mdl(model,z,y,indsit);
        
    elseif choice==8
        done=1;
        cik=0;
        break
    else
        done=1;
        cik=1;
        return
    end
end

scolbw={'color','bw'};
%choice=menu('Make Postscript File?','yes','no');
if strcmp(zaxis,'log') 
    modstr='model_LZ_'; 
else
    modstr='model_Z_';
end
choice =1;
if lmod2==1 modstr=['d',modstr]; end
if choice==1
    eval(['print -depsc ',modstr,char(scolbw(colbw+1)),'.eps']);
    eval(['print -djpeg ',modstr,char(scolbw(colbw+1)),'.jpg']);
end
if cik==0
    yorz=menu('Make rho(z) or rho(y) plot?','rho(y)','rho(z)','no');
else
    yorz=3;
    return
end
if yorz==1
    while 1
        depth=input('which depth? (negative to end): ');
        if depth<0 break; end
        mind=1e9; 
        for iz=1:nz 
            if abs(depth-z(iz))<mind mz=iz; mind=abs(depth-z(iz)); end
        end
        depth=z(mz)
        rhoy = model(mz,:)';
        h1=figure
        plot(y,10.^rhoy,'o-');
        save depth depth y rhoy
        set(gca,'XLim',[-5,y(indsit(nsites))+5]); 
        xlabel('distance[km]');
        ylabel('resistivity[\Omegam]');
        str=[outpath,'/','rhoy_',num2str(z(mz)),'.ps'];
        print(h1,'-dpsc',str);
        
        % Output   z and log10 (rho(y)) to file            % MJU 2016-03-18
        str=[outpath,'/','rhoy_',num2str(z(mz)),'.dat'];   % MJU 2016-03-18
        fid1=fopen(str,'w+');                              % MJU 2016-03-18
                                                           % MJU 2016-03-18
        for i=1:length(y)                                  % MJU 2016-03-18
        fprintf(fid1,'%9.3f %9.3f \n',[y(i);(10^rhoy(i))]);% MJU 2016-03-18  
        end                                                % MJU 2016-03-18
        fclose(fid1);                                      % MJU 2016-03-18
        
        
        
        
    end
elseif yorz==2
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
            rhoz = model(:,indsit(is))';
            semilogx(10.^rhoz,z,'o-'); axis([1,1000,0,240]); 
            title(['Site: ',num2str(is)]);
            set(gca,'YDir','Reverse');
            xlabel('resistivity[\Omegam]');
            ylabel('depth [km]');
        
            % Output   z and log10 (rho(z)) to file
            str=[outpath,'/','rhoz_',num2str(is),'.dat'];
            fid1=fopen(str,'w+');
            fprintf(fid1,'%9.3f %9.3f \n',[z;10.^rhoz]);
            fclose(fid1);
        
            str=[outpath,'/','rhoz_',num2str(is),'.ps'];
            if choice==2 print(h2,'-dpsc',str); end
        
            pause
        end
        if choice==1 break; end
    end
end



% close all

end