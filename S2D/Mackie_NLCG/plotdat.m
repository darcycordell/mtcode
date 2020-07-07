function plotdat(logfreq,sites,dat,err,imode,outpath)
%==========================================================================
% plot Pseudosection part
% dat:  [modes,freq,sites, obs_1 obs_2 cal_1 cal_2]
% err:  [modes,freq,sites, obs_1 obs_2];

disp('IN : plotdat')
global statics lstatic
global err_floor;
global cmode type 
scales;
rhomin=rholims(1);   rhomax=rholims(2);

% set the interpolation type
choice = menu('Fill NaN values?','Yes','No');
if choice == 1; ifill=1; else; ifill=0; end

nfreq=size(dat,2); 
nsites=size(dat,3);

% copy relevant data in new array
tdat=squeeze(dat(imode,:,:,:));
terr=squeeze(err(imode,:,:,:));

%=============================================================================
% Apply static shifts to measured MT impedance data
%=============================================================================
if imode<3
    istatic = menu('Apply static shifts (values as after <dostatics> menu)',...
                   'Yes','No');
    if istatic==1
        for is=1:nsites-1
           tdat(:,is,1)=tdat(:,is,1)-statics(imode,is);
        end
    end
end 

%=============================================================================
% Set shading type
%=============================================================================
while 1
    choice=menu('shading','interpolated','flat','go back to main');
    if choice==1; shade=1; elseif choice==2; shade=2; else; return
end
% modify y and f so that the pcolor boxes 'span' nodes
y=sites;        
nsites=max(size(sites));
freq=logfreq;    nfreq=max(size(freq)); 
dy=diff(y);
yy(1)=y(1)-dy(1)/2;    yy(2:nsites)=y(1:nsites-1)+dy/2; 
yy(nsites+1)=y(nsites)+dy(nsites-1)/2;
y=yy;   clear yy, dy;

nsites=nsites+1;
f=freq;  df=diff(f);
ff(1)=f(1)-df(1)/2;
ff(2:nfreq)=f(1:nfreq-1)+df/2;  ff(nfreq+1)=f(nfreq)+df(nfreq-1)/2;
f=ff;      nfreq=nfreq+1;  clear ff, df;

% set limits 
phamin=0; phamax=90;     hzmin=-0.3; hzmax=0.3;
colbw=0;

% Control axis in pcolor
plot19 = [min(y),max(y),min(f),max(f)];
resid_plot = [-5,5];
done=0;

while ~done
    lims=[-rhomax -rhomin; phamin phamax; hzmin hzmax];
    %======================================================================
    % Plot data and responses
    %======================================================================
    h1 = figure;
    
    if colbw==0
        % define colormap (old!)
        % colmap=menu('color map','John Color','Matlab Jet Color')
        colmap=2;
        if colmap==1;              
            color=[  0    0    0    0  154  255  255  255  255;...
                     0  191  255  255  205  255  220  165    0;...
                   255  255  255    0   50    0    0    0    0]';
            cmap=color/255;
        else
            cmap=jet(64);
            %   cmap(20,:)=[];  cmap(19,:)=[];  cmap(13,:)=[]; 
            %   cmap(7,:)=[];   cmap(1,:)=[];   cmap(1,:)=[];   
        end
    else
        cmap=gray(13);
    end
    colormap(cmap);
        
    % loop over types, e.g.: [obs_rho, obs_pha, cal_rho, cal_pha] 
    type_match = [1 3 2 4];
    for itype=1:4
%        if rem(itype,2)==0 

        ktype=type_match(itype);
        subplot(4,1,itype)
        %subplot(2,1,itype/2) 

        % Fill NaN values
        if ifill==1 & ktype<3
            % Determine the last frequency w/ data
            for is=1:nsites
                if sum(tdat(:,is,ktype))>0
                    no_data=1;
                    nlast(is)=max(find(tdat(:,is,ktype)>0));
                else
                    no_data=0;
                    nlast(is)=nfreq;
                end
                for ifreq=1:nlast(is)
                    if terr(ifreq,is,ktype) >= 1E+5;
                        if ifreq>1 & no_data==0
                            tdat(ifreq,is,ktype) = tdat(ifreq-1,is,ktype);  
                        elseif is>1 | no_data==1
                            tdat(ifreq,is,ktype) = tdat(ifreq,is-1,ktype);  
                        end
                    end
                end    
            end      
        end
    
        if (itype<3 & imode<3) ksign=-1; else ksign=1; end
        
        pcolor(y, f, ksign.*tdat(:,:,ktype));
        axis(plot19);
        ylabel('log(freq)')  % MJU 2016-03-02
        hold on
    
        plot(sites,f(1)*ones(nsites-1,1),'kv','markerfacecolor','k','markersize',5);
        if imode<3 & itype<3   klims=1;  elseif ...
           imode<3 & itype>2   klims=2;  elseif ....
           imode==3            klims=3;  
        end
        caxis(lims(klims,:)); 
        
        if itype==1 title([char(cmode(imode)),' Mode']); end
    
        colorbar('vert');
        if shade==1; shading interp; else; shading flat; end
        hold off
        % end
    end
    xlabel('distance (km)')  % MJU 2016-03-18
    %====================================================
    % Plot residuals
    %====================================================   
    % Note that here, non-shifted data are needed for non shift inversions
    % and shifted data for shift inversions
    if imode<3
        for is=1:nsites-1
            % reverse static
            if istatic==1
                tdat(:,is,1)=tdat(:,is,1)+statics(imode,is);
            end
            % and redo for shift inversions
            if lstatic(imode,is)==1
                tdat(:,is,1)=tdat(:,is,1)-statics(imode,is);
            end
        end 
    end
    h2 = figure;

    % here: type=[rho abs] or [rel img].
    for itype=1:2
        % subplot(211)
        subplot(2,1,itype)
        temp = (tdat(:,:,itype)-tdat(:,:,itype+2))./terr(:,:,itype);
        if colbw==1 temp=abs(temp); end
        
        pcolor(y, f, temp);
        ylabel('log(freq)')  % MJU 2016-03-18
        
        if colbw==0 colormap(cmap); else colormap(flipud(cmap)); end
        
        axis(plot19);
        hold on
        plot(sites,f(2)*ones(nsites-1,1),'kv','markerfacecolor','k','markersize',5);
        title([char(cmode(imode)),' mode residuals: ',...
               upper(char(type(imode,itype)))]);
       
        if colbw==0 caxis(resid_plot); else caxis([0 resid_plot(2)]); end;

        colorbar('vert')
        if shade==1; shading interp; else; shading flat; end
        hold off
    end
 
    xlabel('distance (km)')  %MJU 2016-03-18
    
    % change scales
    if imode<3
        choice=menu('Change:','resistivity limits','phase limits','color/bw','nothing');
        if choice==1
            set_params1=1;
            rhomin=input('Enter min log10(resistivity (Ohm-m)): ');
            while 1
                rhomax=input('Enter max log10(resistivity) > min: ');
                if rhomax>rhomin  break; end
            end
        elseif choice==2
            set_params=1;
            phamin=input('Enter min phase (deg): ');
            while 1
                phamax=input('Enter max phase > min phase: ');
                if phamax>phamin  break; end
            end
        elseif choice==3
            colbw=rem(colbw+1,2);
        else
            done=1;
        end
    else
        choice=menu('Change:','tipper limits','color/bw','nothing');
        if choice==1
            set_params1=1;
            hzmin=input('Enter min tipper: ');
            while 1
                hzmax=input('Enter max tipper > min: ');
                if hzmax>hzmin  break; end
            end
        elseif choice==2
            % change: color<->black and white
            colbw=rem(colbw+1,2);
        else
            done=1;
        end
    end    
end

scolbw={'color','bw'};
% print to file
choice=menu('Make Postscript File?','yes','no');
if choice==1
       disp(['Making PS files for ',char(cmode(imode)),' response']);
       print(h1,'-dpsc',[outpath,'/','nlcg_rsp_',lower(char(cmode(imode))),'_',char(scolbw(colbw+1)),'.ps']);  
       print(h2,'-dpsc',[outpath,'/','nlcg_rsp_',lower(char(cmode(imode))),'_',char(scolbw(colbw+1)),'_err.ps']);
end

% close all
return
end
end