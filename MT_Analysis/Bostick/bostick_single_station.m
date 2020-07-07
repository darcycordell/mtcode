function [p1] = bostick_single_station(fig_number,Z,Zvar,T,rot_init,tlim,rlim,plim,stn_name,lat,long,elev,ns,p_shift)

% Plots 1-D Bostick trnasform for a signle station
% Enci - 2014
rad = 180./pi;       nfreq = length(T);     
p1 = 1;           rot_ang = rot_init;

mkdir bostick

%==========================================================================
% Need to number stations along profile (which is orthogonal to strike)
% Without this they are in alphabetical sequence by name.
% Could change convention here as is usually North (in following line x is east)

x = long ;                              y=lat;
x_mean = mean(x);                       y_mean = mean(y); 
x = cos(y_mean/rad)*111*(x-x_mean);     y = 111*(y-y_mean); % km

c = cosd(rot_ang);      s = sind(rot_ang);     R = [ c, -s ; s, c];

for is=1:ns
  loc = [x(is),y(is)];                 loc = R*loc';
  x_rot(is) = loc(1);                  y_rot(is)=loc(2);
end
% Add loop to put stations in order on monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');

%==========================================================================
is=1;         % Start with station #1
irun = 1; % Variable to loop over stations
while irun == 1  ;  % Loop over MT stations    
    i_sort = index(is);
    for ifreq = 1:nfreq  % Extra impedance for current station from Z 
        Z1(:,:,ifreq) = (Z(:,:,ifreq,i_sort));
        Z1var(:,:,ifreq) = (Zvar(:,:,ifreq,i_sort));
    end
    
    [good,rho_bostick_av,d_av,rho_berd_av,pha_berd_av,rho_bostick_max,rho_bostick_min,d_max,...
        max_dir,min_dir] = bostick2(Z1,Z1var,T);
    if good==0
        disp(['The bostick calculation',' for station ',num2str(is),'-',stn_name{i_sort},' did not perform because of possible negative rho number']);
    else
        figure(fig_number);
        loglog(d_av,rho_bostick_av,'ro');hold on;
        loglog(d_max,rho_bostick_max,'bo');
        loglog(d_max,rho_bostick_min,'go');
        title(['Bostick mapping of station ',num2str(is),'-',stn_name{i_sort}]);
        ylabel ('\rho (\Omega m)');
        xlabel ('depth (m)');
        grid on;    grid minor;
        legend('average','max','min');
        axis([1000,1000000,0,10000]);
        hold off;

        % Save the image to jpg file;
        print('-djpeg100',['bostick\','Stn_',num2str((is)),'-',...
            stn_name{i_sort},'_rho_bostick.jpg']);

        figure(173); 
        subplot(3,1,1);
        semilogx(d_av,rho_bostick_av,'ro');hold on;
        semilogx(d_max,rho_bostick_max,'bo');
        semilogx(d_max,rho_bostick_min,'go');
        title(['Bostick mapping of station ',num2str(is),'-',stn_name{i_sort}]);
        ylabel ('\rho (\Omega m)');
        xlabel ('depth (m)');
        grid on;    grid minor;
        legend('average','max','min');
        axis tight;
    %     axis([1000,1000000,-1000,10000]);
        hold off;

        subplot(3,1,2);
        loglog(T,rho_berd_av,'ro');    
        title(['Berdichevsky average resistivity of station ',num2str(is),'-',stn_name{i_sort}]);
        ylabel ('\rho (\Omega m)');
        xlabel ('period (s)');
        axis tight;
        grid on;    grid minor

        subplot(3,1,3);
        semilogx(T,pha_berd_av,'ro');
        title(['Berdichevsky average phase of station ',num2str(is),'-',stn_name{i_sort}]);

        axis tight;
    %     axis([1,100000,-360,360]);
        ylabel ('phase (deg)')
        xlabel ('period (s)');
        grid on;    grid minor

        % Save the image to jpg file;
        print('-djpeg100',['bostick\','Stn_',num2str((is)),'-',...
            stn_name{i_sort},'_rho_bostick & Berdichevsky averages.jpg']);   


        % Output station name, location ; period, depth, rho to file
        str=['bostick\',stn_name{i_sort},'_rho_bostick.dat'];
        fid2=fopen(str,'w+');
        fprintf(fid2,'\n%s\n',['Stn_',num2str((is))]);
        fprintf(fid2,'\n%s\n',['latitude',num2str(y(is)),'longitude',num2str(x(is))]);
        %the lat, long are wrong, need to be corrected
        fprintf(fid2,'\n%s\n\n',['period','depth','rho_bostick']);
        for ifreq=1:nfreq
            fprintf(fid2,'%9.3f %9.3f %9.3f \n',[T(ifreq);d_av(ifreq);...
                rho_bostick_av(ifreq)]);
        end

        fprintf(fid2,'\n%s\n\n',['depth','rho_bostick_max','max_direction',...
            'rho_bostick_min','min_direction']);
        for idepth=1:length(d_max)
        fprintf(fid2,'%9.3f %9.3f %9.3f %9.3f %9.3f \n',[d_max(idepth);...
            rho_bostick_max(idepth);max_dir(idepth);...
            rho_bostick_min(idepth);min_dir(idepth)]);
        end
        fclose(fid2);
    end
  
    figure (50)
    plot(long,lat,'o'); hold on;
    plot(long(i_sort),lat(i_sort),'ro')
    axis equal

    ichoice = menu(' ','Next station','Previous station','Finish');
%     ichoice=1;
    if ichoice == 1;  is = is+1; end
    if ichoice == 2;  is = is-1; end
    if is < 1;  is=1; end
    if is > ns; is=ns; end
    if ichoice == 3; irun = -1; end
  
end  % End of irun loop

end