function doprofiles(lfreq,sites,dat,err,ymin,ymax,imode,outpath)
% plot profiles
% dat:  [modes,freq,sites, obs_1 obs_2 cal_1 cal_2]
% err:  [modes,freq,sites, obs_1 obs_2];
% loop for ever;

global statics cmode
sym={'b.';'r.';'g.'}; symA={'b-';'r-';'g-'};

nsites=max(size(sites));
iper=input('Enter period to plot (0 to exit): ');
nfreq=max(size(lfreq));

period = 10^(-lfreq(iper))
tdat = squeeze(dat(imode,iper,1:nsites,:));
terr = squeeze(err(imode,iper,1:nsites,:));
 
%======================================================================
% Apply static shifts to measured MT impedance data
choice = menu('Apply static shifts (values as after <dostatics> menu)',...
              'Yes','No');
if choice == 1
    for is=1:nsites
        tdat(is,1)=tdat(is,1)-statics(imode,is);
    end
end

if imode<3
    tdat(:,1:2:3) = 10*tdat(:,1:2:3);
    terr(:,1) = 2*tdat(:,1).*terr(:,1);
end

figure('PaperPositionMode','auto');
subplot(2,1,1)

if imode<3
    semilogy(sites,tdat(:,3),char(symA(imode)));
    hold on
    errorbar(sites,tdat(:,1),terr(:,1),char(sym(imode)))
    hold off
    axis([ymin ymax 0.3 1000])
    ylabel([char(cmode(imode)),' apparent resistivity']);
else
    errorbar(sites,tdat(:,1),terr(:,2),char(sym(imode)))
    hold on
    plot(sites,tdat(:,3),char(symA(imode)))
    axis([ymin ymax -.8 .8])
    ylabel([char(cmode(imode)),' real part']);
end
xlabel('distance(km)')
    
subplot(2,1,2)

errorbar(sites,tdat(:,2),terr(:,2),char(sym(imode)))
hold on
plot(sites,tdat(:,4),char(symA(imode)))
hold off
xlabel('distance(km)')
if imode<3
    axis([ymin ymax 0 90])
    ylabel([char(cmode(imode)),' phase (degrees)']);
else
    axis([ymin ymax -.8 .8])
    ylabel([char(cmode(imode)),' imaginary part']);
end

%set(gca,'YGrid','on');
eval(['print -dpsc ',outpath,'/','profile_f',num2str(iper),'_',lower(char(cmode(imode))),'.ps']);
pause

%close all;
end