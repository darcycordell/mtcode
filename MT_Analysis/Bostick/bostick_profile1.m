function [p1] = bostick_profile1(fig_number,Z,Zvar,T,azimuth,tlim,rlim,plim,dlim,stn_name,lat,long,elev,maplim,ns,profile_dir)

% Simple version of Bostick pseudo section
% Just does Bostick transform (not find maximum rho)

rad = 180./pi;       nfreq = length(T);
mu = 4*pi*1e-7;      p1 = 1;
prof_azim = mean(azimuth);

mkdir bostick

% Process co-ords
x = long; y=lat;
x_mean = mean(x);  y_mean = mean(y); % km
x = cos(y_mean/rad)*111*(x-x_mean);y = 111*(y-y_mean); % km

irun = 1;

while irun ==1 ;

c = cosd(prof_azim);      s = sind(prof_azim);
R = [ c, -s ; s, c];

for is=1:ns
    loc = [x(is),y(is)];
    loc = R*loc';
    x_rot(is) = loc(1);y_rot(is)=loc(2);
end
    
% Add loop to put stations in order on monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');
     
for is =1:ns
    i_sort = index(is);  
    Z1 = Z(:,:,:,i_sort);    Z1var = Zvar(:,:,:,i_sort);   
    [rho_bostick_av(:,is),depth_bostick(:,is),~,~] = bostick1(Z1,T);
end

%  using interpolation and pcolor to plot the average bostick resistivity model

logdepth=log10(depth_bostick);
d_interp=min(min(logdepth)):0.1:max(max(logdepth));
nd = length(d_interp);
  
% Bostick_pseudo=zeros(length(dsequence),ns);  
bostick_pseudo_av= zeros(nd,ns+1);  % ns+1 to add the last column to the map
  
for is=1:ns
    x_gridding = logdepth(:,is); y_gridding = rho_bostick_av(:,is);
    x_gridding(isnan(y_gridding))=[]; % Remove the NaN values from the vector
    y_gridding(isnan(y_gridding))=[]; % Remove the NaN values from the vector
    if length(x_gridding)>5
        bostick_pseudo_av(:,is)=interp1(x_gridding,y_gridding,d_interp,'nearest',NaN);
    else
        bostick_pseudo_av(:,is)=nan(size(d_interp));
    end
end

% bostick_pseudo_av =   abs(bostick_pseudo_av);


figure(fig_number);

% move the stations to the middle of each column
x_sort_m=x_sort;
for j=1:ns-1
    x_sort_m(j)=(x_sort(j)+x_sort(j+1))/2;
end
x_sort_m(ns)=x_sort(ns)+(x_sort(ns)-x_sort(ns-1))/2;

subplot(2,1,1);
% give the last column the same width with the second last column using
% function 2*x_sort(ns)-x_sort(ns-1)
pcolor([x_sort,2*x_sort(ns)-x_sort(ns-1)],-(10.^d_interp)/1000,log10(bostick_pseudo_av)); 
shading 'flat';title(['\rho_{av}  \theta =', num2str(prof_azim),'^o']);
ylabel ('depth (km)');
xlabel ('distance (km)');
axis([min(x_sort_m),max(x_sort_m),dlim(1),dlim(2)])
caxis([log10(rlim(1)),log10(rlim(2))]);
hold on;  
colorbar('vert');colormap(flipud(jet)); 



plot(x_sort_m,zeros(1,ns),'kv'); hold off;


subplot(2,1,2);
% give the last column the same width with the second last column using
% function 2*x_sort(ns)-x_sort(ns-1)
pcolor([x_sort,2*x_sort(ns)-x_sort(ns-1)],-d_interp,log10(bostick_pseudo_av));   
shading 'flat';title(['\rho_{av}  \theta =', num2str(prof_azim),'^o']);
ylabel ('-log_{10} depth (m)');
xlabel ('distance (km)');
axis([min(x_sort_m),max(x_sort_m),-log10(abs(1000*dlim(1))),-log10(abs(1000*dlim(2)))])
caxis([log10(rlim(1)),log10(rlim(2))]);
colorbar('vert');colormap(flipud(jet)); 
hold on;  

plot(x_sort_m,-min(min(logdepth)),'kv'); hold off;


eval('print -djpeg100 bostick\pseudo_bostick_average.jpg');
%==========================================================================  
% Output to dat file
str1=['bostick\','rho_bostick_pseudo_av.dat'];
fid1=fopen(str1,'w+');
for is=1:ns
    fprintf(fid1,'\n%s\n',['Stn_',num2str(is)]);
    fprintf(fid1,'\n%s\n\n',['depth/','rho_bostick_average']);
    for ifreq=1:length(d_interp)
        fprintf(fid1,'%9.3f %9.3f\n',[10.^d_interp(ifreq);bostick_pseudo_av(ifreq,is)]);
    end
end
fclose(fid1);

ichoice = menu(' ','Increase rot angle','Decrease rot angle','Finish');
if ichoice == 1;  prof_azim = prof_azim+5; end
if ichoice == 2;  prof_azim = prof_azim-5; end
if ichoice == 3;  irun = -1; end
  
end 

end