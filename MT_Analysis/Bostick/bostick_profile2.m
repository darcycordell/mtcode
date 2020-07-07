function [p1] = bostick_profile2(Z,Zvar,T,azimuth,tlim,rlim,plim,dlim,stn_name,lat,long,elev,maplim,ns,profile_dir)
% Feb 13 2014  
% Plot Bostick pseudosection of MT data

rad = 180./pi;     p1 = 1;
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
    [rho_bostick_av(:,is),d_av(:,is),~,~] = bostick1(Z1,T);
end

%==========================================================================
% Using interpolation and pcolor to plot the average bostick resistivity model
logdepth=log10(d_av);
d_interp=min(min(logdepth)):0.1:max(max(logdepth));
nd = length(d_interp);

figure;


%==========================================================================
% Calculate the maximum resistivity corresponding to each layer in d_interp
for is=1:ns
    i_sort = index(is);  
    Z1 = Z(:,:,:,i_sort);    Z1var = Zvar(:,:,:,i_sort);  
    [~,~,~,~,~,rho_bostick_max(:,is),~,d_max,~,~] = bostick2(Z1,Z1var,T,d_interp);
end
bostick_pseudo_max=zeros(nd,ns+1);  % ns+1 to add the last column to the map
bostick_pseudo_max(:,1:ns)=rho_bostick_max;


% move the stations to the middle of each column
x_sort_m=x_sort;
for j=1:ns-1
    x_sort_m(j)=(x_sort(j)+x_sort(j+1))/2;
end
x_sort_m(ns)=x_sort(ns)+(x_sort(ns)-x_sort(ns-1))/2;

subplot(2,1,1);
pcolor([x_sort,2*x_sort(ns)-x_sort(ns-1)],-d_max/1000,log10(bostick_pseudo_max)); 
shading 'flat';title(['\rho_{max}  \theta =', num2str(prof_azim),'^o']);
axis([min(x_sort_m),max(x_sort_m),dlim(1),dlim(2)]);
ylabel ('depth (km)');
xlabel ('distance (km)');
colorbar('vert');colormap(flipud(jet)); 
caxis([log10(rlim(1)),log10(rlim(2))]);hold on;  


plot(x_sort_m,-min(min(logdepth)),'kv'); hold off;

subplot(2,1,2);
pcolor([x_sort,2*x_sort(ns)-x_sort(ns-1)],-log10(d_max),log10(bostick_pseudo_max));   
shading 'flat';title(['\rho_{max}  \theta =', num2str(prof_azim),'^o']);
ylabel ('-log10 depth (m)');
xlabel ('distance (km)');
axis([min(x_sort_m),max(x_sort_m),-log10(abs(1000*dlim(1))),-log10(abs(1000*dlim(2)))]);
caxis([log10(rlim(1)),log10(rlim(2))]);
colorbar('vert');colormap(flipud(jet)); 
hold on;  
plot(x_sort_m,-min(min(logdepth)),'kv'); hold off;

eval('print -djpeg100 bostick\pseudo_bostick_max.jpg');
%==========================================================================  
% Output to dat file
str1=['bostick\','rho_bostick_pseudo_max.dat'];
fid1=fopen(str1,'w+');
for is=1:ns
    fprintf(fid1,'\n%s\n',['Stn_',num2str(is)]);
    fprintf(fid1,'\n%s\n\n',['depth/','rho_bostick_max']);
    for ifreq=1:length(d_interp)
        fprintf(fid1,'%9.3f %9.3f\n',[10^d_interp(ifreq);bostick_pseudo_max(ifreq,is)]);
    end
end
fclose(fid1);

ichoice = menu(' ','Increase rot angle','Decrease rot angle','Finish');
if ichoice == 1;  prof_azim = prof_azim+5; end
if ichoice == 2;  prof_azim = prof_azim-5; end
if ichoice == 3;  irun = -1; end

end

end 