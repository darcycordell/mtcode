function docurves(lfreq,sites,dat,err,imode,outpath)
% plot sounding curves:     
% -> calls: doplot
% dat:  [modes,freq,sites, obs_1 obs_2 cal_1 cal_2]
% err:  [modes,freq,sites, obs_1 obs_2];

disp('IN : docurves')

global cmode
%==========================================================================
% Apply static shifts to measured MT impedance data
istatic = menu('Apply static shifts (values as after <dostatics> menu)',...
               'Yes','No');
%==========================================================================
% loop for ever;
while 1
    choice=menu('plot option','plot all','plot single site','all at once','go back to main');
    if choice==1
        plotall=1;
    elseif choice==2
        plotall=0;
    elseif choice==3
        plotall=3;
    else
        return
    end
    nsites=length(sites);
    % figure out how many sites need to be plotted
    if plotall==1 | plotall==3;
        kstart=1;
        kend=nsites;
    else
        while 1
            isite=input('Enter site number to plot (0 to exit): ');
            if ~(isite>nsites) break; end
             fprintf(2,'%s\n','Invalid entry: isite > nsites')
        end
        if isite==0 return; end;
        kstart=isite;
        kend=isite;
    end

    nfreq=length(lfreq);
    % plot curves
    for isite=kstart:kend
      
        if (imode == 1)
          figure(20)
          doplot(isite,lfreq,sites,dat,err,1,istatic,outpath);
        end
        
        if (imode == 2)
          figure(20)
          doplot(isite,lfreq,sites,dat,err,2,istatic,outpath);
        end    
        
        if (imode == 3)
          figure(20)
          doplot(isite,lfreq,sites,dat,err,3,istatic,outpath);
        end    
              
        if imode ==4
             figure(20)
            doplot(isite,lfreq,sites,dat,err,1,istatic,outpath);
            doplot(isite,lfreq,sites,dat,err,2,istatic,outpath);
            doplot(isite,lfreq,sites,dat,err,3,istatic,outpath); %MJU commeted
        end
        
        eval(['print -depsc ','site_',num2str(isite),'_',...
             lower(char(cmode(imode)))]);
         if plotall==1 |plotall==0;
             pause
         end
        close(20);
    end 
end
end
