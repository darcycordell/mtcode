function pmv20_iso

% Plots the results of the ISOTROPIC NLCG6 inversion
% For full version history see separate Word document

% pmv20_iso updates all prior versions. 

% pmv20_aniso plots the results of the anisotropic inversion

% Major changes from pmv19 are (a) functions in separate files
% and a separate script is used for ISOTROPIC and ANISOTROPIC
% inversions.

close all
clear all

global statics lstatic
global err_floor parfile
global cmode type kmode version seed elev


cmode = {'TM';'TE';'HZ';'TMTE'}; 
type = {'rho','pha';'rho','pha';'rel','img'};
%version = '6_9';
version = '6_11';

[parfile,outpath] = uigetfile({'*.par'},'pick inversion parameter file');
parfile=[outpath,'/',parfile];
if exist(parfile,'file')
    disp(['reading parameter file: ',parfile]);
    [seed,kmode,modfile,datfile,tau,err_floor,statshift,indsitmod] = ...
        load_par_iso(parfile);
else
    error([parfile,' does not exist!'])
end
seedp=[outpath,seed,'_',version];
modfile=[outpath,'/',modfile];

% convert error floors to log10 scale
% drho in log10 space: 
% dlog10(rho) = d(ln(rho))/ln(10) = rel.d(rho)/ln(10)
err_floor(1:2,1) = 0.434.*err_floor(1:2,1)./100;
% rund2inv doc tells me (ws) that phase error is in percent of 
% equivalent ln(rho), and thus: 1/(2*100)*180/pi = 0.286 in deg
err_floor(1:2,2) = .286.*err_floor(1:2,2);
% don't change Hz error floor.

% load the response file which will be called dat
if exist([seedp,'.rsp'],'file')
    disp(['reading response file: ',seedp,'.rsp']);
    [dat,per,indsit,sitrms] = load_nlcgout([seedp,'.rsp'],kmode,indsitmod);
else
    seedp(length(seedp))='1';
    if exist([seedp,'.rsp'],'file')
        disp(['reading response file: ',seedp,'.rsp']);
        [dat,per,indsit,sitrms] = load_nlcgout([seedp,'.rsp'],kmode,indsitmod);
    else
        seedp(length(seedp))='0';
        if exist([seedp,'.rsp'],'file')
            disp(['reading response file: ',seedp,'.rsp']);
            [dat,per,indsit,sitrms] = load_nlcgout([seedp,'.rsp'],kmode,indsitmod);
        else
            seedp(length(seedp)-1:length(seedp))='9 ';
            if exist([seedp,'.rsp'],'file')
                disp(['reading response file: ',seedp,'.rsp']);
                [dat,per,indsit,sitrms] = load_nlcgout([seedp,'.rsp'],kmode,indsitmod);
            else
%                 error(['File ',seedp,'.rsp', ' does not exist!'])
            end
%             error(['File ',seedp,'.rsp', ' does not exist!'])
        end
%         error(['File ',seedp,'.rsp', ' does not exist!'])
    end
end
% get freq
freq = 1./per; logfreq = log10(freq);
nfreq = length(freq);
nsites=length(indsit);

%Define elev=0, overwrite if file elevations.txt exists
elev(1:length(indsit))=0;
if exist([outpath 'elevation.txt'],'file')
    elev=load([outpath 'elevation.txt']);
end

% load the model file
if exist(modfile,'file')
    [ny,nz,dy,dz,model] = load_mdl(modfile,version);
else
    error(['File: ',modfile,' does not exist'])
end

% read objfun (rms, roughness...) as a function of iteration
nit=0;objfun=zeros(1,6);
if exist([seedp,'.log'],'file')
    [objfun] = load_log_iso([seedp,'.log'],err_floor);
    nit=size(objfun,1)-1;
else
    warning(['File: ',[seedp,'.log'],' does not exist'])
end

% compute y's and z's at block edges
y=zeros(1,ny+1);
for k=1:ny  y(k+1)=y(k)+dy(k); end
z=zeros(1,nz+1);
for k=1:nz  z(k+1)=z(k)+dz(k); end

% compute y's at sites
ysite=y(indsit)+dy(indsit)'/2;
diff(ysite);

y0=0;
% sites, y: site location wrt 1st site (+y0) in km
sites=ysite-ysite(1)+y0;
sites=sites/1000;

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
    dd=cumsum(diff(dd));dd=[0.001;dd];
    sites;
    %plot(sites,'o');hold on;
    sites=dd';
end
%comment out above lines to plot stations in the middle of each cell
sites;
%plot(sites,'*')
y=sites;

save y_dist.mat y sitrms

% Initialize statics
lstatic = zeros(2,nsites);
statics = zeros(2,nsites); 

% set defaults
ymin=min(y);     	ymax=max(y);             dy=ymax-ymin;
ymin=ymin-.1*dy;	ymax=ymax+.1*dy;
rhomin=0;        	rhomax=3;                
phamin=0;        	phamax=90;

% load errors matrix
err = zeros(3,nsites,nfreq,2);
for imode=1:3
    if kmode(imode)==1
        name=[outpath,'/',char(datfile(imode))];
        if exist(name,'file')
            [err]=load_err(err,name,imode,indsit,indsitmod,per);
        else
            warning(['File: ',name,' does not exist'])
        end
    end
end
% modes,freqs,sites,[meas_1 meas_2]
err=permute(err,[1 3 2 4]);

% load static shifts
lstatic=zeros(2,nsites);
statics=zeros(2,nsites);
if statshift == 'y'
    [lstatic,statics]=load_stat(kmode,cmode,nsites,indsit,indsitmod,outpath);
    save statics_1.mat statics sitrms y                        
end

dat=reshape(dat(:,:,1:12),[nsites,nfreq,4,3]);
dat=permute(dat,[4 2 1 3]);
% modes,freqs,sites,[obs_1 obs_2 meas_1 meas_2],mode

% high (>1E+5) errors -> NaN
dat(err>1E+5)=NaN;
err(err>1E+5)=NaN;

% modes are 'found' when there's at least one data point ~= NaN
no_nan = ~isnan(dat);
kmode = sum(sum(sum(no_nan,4),3),2)>0 ; % e.g. [1 0 1] <= TM and HZ modes inverted
mtmode = find(kmode(1:2)>0);          % e.g. [1] ...

% scale error of rho & pha & tzy
% d(ln(rho)) => d(log10(rho)) = d(ln(rho))/ln(10)
err(mtmode,:,:,1) = 0.43.*err(mtmode,:,:,1);
err(mtmode,:,:,2) = (180/pi).*err(mtmode,:,:,2);
% err(3,:,:,:)   = err(3,:,:,:).*0.707;
% I (ws) guess the error in the input file does not refer to the 
% absolute error of the complex number but its real and imag parts

% log10(rho), pha=-pha
dat(1:2,:,:,[1:2:3]) = log10(dat(1:2,:,:,[1:2:3]));
dat(1:2,:,:,[2:2:4]) = -dat(1:2,:,:,[2:2:4]);

%Block added by TB May23-2008
%Writes file data.mat with data from response file (*.rsp) organized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periods=per';
for j=1:nfreq
    for k=1:nsites
        obsTMrho(j,k)=10.^dat(1,j,k,1);
        obsTMphi(j,k)=dat(1,j,k,2);
        calTMrho(j,k)=10.^dat(1,j,k,3);
        calTMphi(j,k)=dat(1,j,k,4);
        obsTErho(j,k)=10.^dat(2,j,k,1);
        obsTEphi(j,k)=dat(2,j,k,2);
        calTErho(j,k)=10.^dat(2,j,k,3);
        calTEphi(j,k)=dat(2,j,k,4);
        obsHZrel(j,k)=dat(3,j,k,1);
        obsHZimg(j,k)=dat(3,j,k,2);
        calHZrel(j,k)=dat(3,j,k,3);
        calHZimg(j,k)=dat(3,j,k,4);
    end
end
save([outpath 'data.mat'],'obsTMrho','obsTMphi','calTMrho','calTMphi',...
                          'obsTErho','obsTEphi','calTErho','calTEphi',...
                          'obsHZrel','obsHZimg','calHZrel','calHZimg',...
                          'periods');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shift responses to real responses for static shift inversion
for imode=1:2
    for is=1:nsites
        if lstatic(imode,is)==1
            dat(imode,:,is,3)=dat(imode,:,is,3)-statics(imode,is);
        end
    end
end
% extend data matrices for pcolor
err(:,nfreq+1,:,:) = err(:,nfreq,:,:);
err(:,:,1+nsites,:) = err(:,:,nsites,:);
dat(:,nfreq+1,:,:) = dat(:,nfreq,:,:);
dat(:,:,nsites+1,:) = dat(:,:,nsites,:);

choice=menu('Apply error floor to data?','Yes','No');
%=============================================================================
% Apply error floors
%=============================================================================
if choice==1
    for imode=1:3
        
        if kmode(imode)
            for itype=1:2
                terr = err(imode,:,:,itype);
                lindex = find(terr<=err_floor(imode,itype));
                hindex = find(~isnan(terr));
                terr(lindex)=err_floor(imode,itype);
                err(imode,:,:,itype)=terr;
                nlow = length(lindex); nhigh = length(hindex);
                if nhigh>0 percent=100*nlow/nhigh; else percent=0; end
                disp(['Percent of ',upper(char(type(imode,itype))),...
                  ' data (',char(cmode(imode)),' mode) below error floor ',...
                   num2str(err_floor(imode,itype)),': ',num2str(percent)]);
            end
        end
    end 
end

subset_rms=0;
if subset_rms
    % get subset rms
    sit1=9; freq1=1;     
    sit2=18; freq2=nfreq; 

    tdat=dat(:,freq1:freq2,sit1:sit2,1:2)-dat(:,freq1:freq2,sit1:sit2,3:4);
    terr=err(:,freq1:freq2,sit1:sit2,:);
    res=tdat./terr;
    N=length(find(~isnan(res) & ~isinf(res)))
    res(find(isnan(res)))=0; res(find(isinf(res)))=0;
    s=sum(sum(sum(sum(res.^2,4),3),2));
    rms=sqrt(s/N)
end

hmode=[];
for imode=1:3
    if kmode(imode)>0 hmode=[hmode,' ',char(cmode(imode))]; end
end

% loop for ever
cik=0;
while cik==0
    choice=menu('plot option',...
                'pseudosections','sounding curves','profiles',...
                'rms','statics','resistivity model','sensitivity map',...
                'objective function','conductance','quit');
    switch choice
        case {1,2,3,4}
            lindex=find(kmode>0);
            if length(lindex) > 1
                % TM and TE mode? => TMTE (together)
                cmmode=cmode(lindex);
                if choice==2 & find(lindex==1) & find(lindex==2) 
                    cmmode{length(lindex)+1}='TMTE';
                elseif choice==4
                    cmmode{length(lindex)+1}='All';
                end
                imode = menu('choose mode',cmmode);
                if imode==length(lindex)+1
                    imode=4;
                else
                    imode = lindex(imode);
                end
            else
                imode = lindex(1);
            end
            if choice==1
                plotdat(logfreq,sites,dat,err,imode,outpath)
            elseif choice==2
                docurves(logfreq,sites,dat,err,imode,outpath)
            elseif choice==3
                doprofiles(logfreq,sites,dat,err,ymin,ymax,imode,outpath)
            elseif choice==4
                % rms
                plotrms(dat,err,imode,outpath)
            end
        case 5
            % statics
            dostatics(dat,err,outpath)
        case 6
            header=['model: ',seed,' |',hmode,' | tau=',num2str(tau),...
                    ' | rms=',num2str(objfun(nit+1,2)),' | ',num2str(nit),' iters'];
            %header=''
            cik=nlcg_mdl_iso(seedp,header,outpath);
        case 7
            header=['sens: ',seed,' |',hmode,' | tau=',num2str(tau),...
                         ' | rms=',num2str(objfun(nit+1,2)),' | ',num2str(nit),' iters'];
            %header=''
            nlcg_sns(seedp,header,outpath)
        case 8
            plot_objfun(objfun,outpath)
        case 9 %conductance
            header=['conductance: ',seed,' |',hmode,' | tau=',num2str(tau),...
                    ' | rms=',num2str(objfun(nit+1,2)),' | ',num2str(nit),' iters'];
            %header=''
            cik=nlcg_mdl_iso(seedp,header,outpath);
        case 10
            % close all
            return
    end
end
end