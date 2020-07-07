function plotrms(dat,err,imode,outpath)

global cmode
global statics lstatic 

nsites=size(dat,3)

% Note that here, as in the residual plot of plotdat,
% non-shifted data are needed for non shift inversions
% and shifted data for shift inversions
tdat=dat;
if imode<3
    for is=1:nsites-1
        if lstatic(imode,is)==1
            tdat(imode,:,is,1)=tdat(imode,:,is,1)-statics(imode,is);
        end
    end 
end

% get rms
if imode==4
    tdat=tdat(:,:,:,1:2)-tdat(:,:,:,3:4);
    terr=err(:,:,:,:);
    ttitle=['rms all modes'];
else
    tdat=tdat(imode,:,:,1:2)-tdat(imode,:,:,3:4);
    terr=err(imode,:,:,:);
    ttitle=['rms ',char(cmode{imode}),' mode'];
end
res=tdat./terr;

% data points
N = ~isnan(res) & ~isinf(res);
% set others to zero so that they fall out of summation
res(find(isnan(res)))=0; res(find(isinf(res)))=0;
% sum over all but sites
N=sum(sum(sum(N,1),2),4);
s=sum(sum(sum(res.^2,1),2),4);
rms=squeeze(sqrt(s./N));

% plot
figure
plot([1:nsites-1],rms(1:nsites-1),'b+');
hold on
plot([0 nsites+1],[1 1],'k:');
    
ylabel('rms');
xlabel('site number');
title(ttitle);

axis([0 nsites 0 max(rms)+1]);
if imode==4 fname='rms_all'; else fname=['rms_',char(cmode{imode})]; end
eval(['print -dpsc ',fname,'.ps']);
end