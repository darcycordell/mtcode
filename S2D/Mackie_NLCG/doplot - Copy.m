function doplot(isite,lfreq,sites,dat,err,imode,istatic,outpath)
%=============================================================================
% dat:  [modes,freq,sites, obs_1 obs_2 cal_1 cal_2]
% err:  [modes,freq,sites, obs_1 obs_2];

% subplot commands changed my MJU 2018-03-14

disp('IN : doplot')
imode
global cmode

if istatic==1 % Only load statics file if static shifts are chosen to be applied, DR 19/04/2010
    load statics_1
else % Otherwise load y_dist file (No guarantee that the rms in this file is correct though), DR 19/04/2010
    load y_dist
end
lper=-lfreq;
nper=length(lper);

% note that dat is log10(rho)
tdat=squeeze(dat(imode,1:nper,isite,:));
terr=squeeze(err(imode,1:nper,isite,:));

%tmin=floor(min(lper));
tmin=-2;
tmax=3;
%tmax=ceil(max(lper));
sym={'bo';'ro';'g.'}; symA={'b';'r';'g'};
    


if imode<3 
    h1=subplot(12,1,[1:6]);
    %======================================================================
    % Apply static shifts to measured MT impedance data
    if istatic==1
        tdat(:,1)=tdat(:,1)-statics(imode,isite);
    end
    %=============================================
    %errorbar(lper,tdat(:,1),terr(:,1),terr(:,1),char(sym(imode)),'markersize',3);
    plot(lper,tdat(:,1),char(sym(imode)),'markersize',3);
    lper;%
    tdat(:,1);%
    terr(:,1);%
    hold on
%%%%plot(lper,tdat(:,3),char(symA(imode)),'linewidth',1);
    ylabel('Apparent Resistivity (\Omega m)');
    %axis equal
    axis([tmin-0.2 tmax+0.2 0 3]);
    set(gca,'YTick',[-1 0 1 2 3 4 5 6],'XTick',[-4 -3 -2 -1 0 1 2 3 4 5]);    hold on
    set(gca,'YTicklabel',[0.1 1 10 100 1000 10000 100000 1000000]);
    set(gca,'XTicklabel',[]);

    sites=round(sites*10)/10;
    if istatic==1 % Only calculate static shift information if static shifts are chosen to be applied, DR 19/04/2010
        tits=round(statics*100)/100;
        sitrms=round(sitrms*10)/10;
        title([' Site-',num2str(isite),', TM-statics=',num2str(tits(1,isite)),', TE-statics=',num2str(tits(2,isite)),', rms=',num2str(sitrmss(isite))]);
    else % Otherwise use simplified title without static shift information, DR 19/04/2010
        sitrms=round(sitrms*10)/10;
        title([' Site-',num2str(isite),', rms=',num2str(sitrms(isite))]);
    end
    box on; grid on; %text(tmin,3.5,['rms= ',num2str(sitrms(isite))])
    
    
    h2=subplot(12,1,[7:10]);
    %errorbar(lper,tdat(:,2),terr(:,2),terr(:,2),char(sym(imode)),'markersize',3);
    plot(lper,tdat(:,2),char(sym(imode)),'markersize',3);
    hold on
%%%%plot(lper,tdat(:,4),char(symA(imode)),'linewidth',1);
%   xlabel('Period (sec.)');
    ylabel('Phase (deg)');
    xlabel('Period (s');
    axis([tmin-0.2 tmax+0.2  -10 100]);
    set(gca,'YTick',[0 45 90],'XTick',[-4 -3 -2 -1 0 1 2 3 4 5]);    hold on
    set(gca,'YTicklabel',[0 45 90]);
    %set(gca,'XTicklabel',[]);   
    set(gca,'XTicklabel',[0.0001 0.001 0.01 0.1 1 10 100 1000 10000 100000]);
    box on; grid on
    eval(['print -depsc ','site_',num2str(isite),'_',lower(char(cmode(imode)))]);
else
    h3=subplot(12,1,[11:12]);
    errorbar(lper,tdat(:,1),terr(:,1),terr(:,1),'b.','markersize',5);    hold on
    errorbar(lper,tdat(:,2),terr(:,2),terr(:,2),'r.','markersize',5);
    %legend('real','imaginery')
    errorbar(lper,tdat(:,3),'b','linewidth',1);
    errorbar(lper,tdat(:,4),'r','linewidth',1);
    xlabel('Period (sec');
    ylabel('Tipper');
    axis([tmin-0.2 tmax+0.2  -.6 .6]);
    set(gca,'YTick',[-0.5 0 0.5],'XTick',[-4 -3 -2 -1 0 1 2 3 4 5]);    hold on
    set(gca,'YTicklabel',[-0.5 0 0.5]);
    set(gca,'XTicklabel',[0.0001 0.001 0.01 0.1 1 10 100 1000 10000 100000]);   
    box on; grid on
    eval(['print -depsc ','site_',num2str(isite),'_',lower(char(cmode(imode)))]);
end
end