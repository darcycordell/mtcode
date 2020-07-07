function dostatics(dat,err,outpath)

%Plot static shift as a function of station number
global statics lstatic sitrms
global kmode cmode

nsites=size(statics,2);
% any data not been inverted for static shift?
if length(find(lstatic(:,:)==0))>0
    choice=menu('Data are at least partly not inverted for static shift',...
                'calculate statics','do not calculate statics (default)');
    if choice==1
        % note that dat is log10(rho)
        % data minus response
        stat=dat(:,:,:,1)-dat(:,:,:,3);
        % set NaNs to zero
        stat(find(isnan(stat)))=0;
        % indices of actual data points
        N=stat~=0;
        for imode=1:2
            if kmode(imode)==1
                % indices of data that not inverted for static shift
                no_inv=find(lstatic(imode,:)==0);
                % sum over frequencies, but not sites
                tstat=sum(stat(imode,:,:),2);
                tN=sum(N(imode,:,:),2);
                statics(imode,no_inv)=tstat(no_inv)./tN(no_inv);
            end
        end
    elseif choice==2
        % (re)set to zero
        no_inv=find(lstatic==0);
        statics(no_inv)=0;
    end
end

nsites = size(statics,2);
ls={'+','o'}; cs={'b','r'}; ;
figure
ttitle=[];
for imode=1:2
    if kmode(imode)==1
        % indices of data that not inverted for static shift
        no_inv=find(lstatic(imode,:)==0)
        plot(no_inv,statics(imode,no_inv),[char(ls{imode}),char(cs{1})]);   
        ttitle=[ttitle,'  ',char(ls{imode}),' = ',char(cmode{imode}),' mode  '];
        hold on
        inverted=find(lstatic(imode,:)==1);
        plot(inverted,statics(imode,inverted),[char(ls{imode}),char(cs{2})]);   
    end
end
if length(find(lstatic==1))>1
    ttitle=[ttitle,'  red: inverted '];
end     

plot ([0 nsites+1],[0 0],'k:');    
ylabel('static');
xlabel('site number');
title(ttitle);
axis([0 nsites+1 -2. 2.0]);
load y_dist.mat
save statics_1.mat y statics sitrms
clear y;
delete y_dist.mat;
eval(['print -dpsc ','statics.ps']);
end