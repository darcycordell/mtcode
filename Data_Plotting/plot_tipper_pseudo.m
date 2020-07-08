function plot_tipper_pseudo(d)
%
% Function which plots real and imaginary tipper pseudosections for zx and
% zy modes
%
% Usage: plot_tipper_pseudo(d)
%
% "d" is an MT data structure
%
%%
u = user_defaults;

%Pull a profile to plot the pseudo-section. "sidx" is stations indices
%included on profile and "rot_ang" is the direction of the normal to the
%profile direction
[sidx,midpoint,azimuth] = get_pseudo_section_indices(d);

rad = 180./pi;     

% Process cooordinates
if isfield(d,'loc') && ~all(d.origin==0) % here x is E-W, y is N-S 
%     lon_mean = mean(d.loc(sidx,2));  lat_mean = mean(d.loc(sidx,1)); % km
    lon_mean = midpoint(1);  lat_mean = midpoint(2); % km
    x = cos(lat_mean/rad)*111*(d.loc(sidx,2)-lon_mean);y = 111*(d.loc(sidx,1)-lat_mean); % km
else
    % if data from 2D inversion, then there is no reference to geographic
    % coordinates. just use x, y from inversion.. but need to switch x and
    % y for plotting
    x = d.y; y = d.x;    
end
N = length(sidx);

if ~u.rotate_data_to_azimuth
    dr=d;
else
    [dr]=rotate_d(d,azimuth); %Rotate data   
end

output_data_menu = menu('Output tipper data along profile to text files?','Yes','No');

irun = 1; % Variable to loop over azimuth

while irun == 1  % Loop over azimuth
%%
    c = cosd(azimuth);      s = sind(azimuth);     R = [ c, -s ; s, c]; % this is a positive-counterclockwise rotation matrix
    %Rotate station coordinates to get distance along profile
    x_rot = zeros(N,1); y_rot = zeros(size(x_rot));
    for is=1:N
        loc = [x(is),y(is)];    loc = R*loc';
        x_rot(is) = loc(1);     y_rot(is)=loc(2);
    end
    
    % Add loop to put stations in order on monotonically increasing x_rot
    [x_sort,index] = sort(x_rot,'ascend');
    index = sidx(index);  
%%
    %========================================================================
    %Plot real tipper zx
    set_figure_size(1);   
    subplot(2,2,1)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(real(dr.tip(:,1,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
    else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
        title('Real T_{zx}')
    end
    
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(flipud(u.cmap)) % make large positive values red
    caxis(u.tiplim)
    hcb = colorbar;
    hcb.Label.String = 'Real T_{zx}';
    set(gca,'Layer','top')
    
    %========================================================================
    %Plot real tipper zy
    subplot(2,2,2)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(real(dr.tip(:,2,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
    else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
        title('Real T_{zy}')
    end 
    
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    caxis(u.tiplim)
    hcb = colorbar;
    hcb.Label.String = 'Real T_{zy}';
    set(gca,'Layer','top')
        
    %========================================================================
    %Plot imaginary tipper zx
    subplot(2,2,3)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(imag(dr.tip(:,1,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
    else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
        title('Imag T_{zx}')
    end 
    
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    caxis(u.tiplim)
    hcb = colorbar;
    hcb.Label.String = 'Imag T_{zx}';
    set(gca,'Layer','top')
       
    %========================================================================
    %Plot imaginary tipper zy
    subplot(2,2,4)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(imag(dr.tip(:,2,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
    else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
        title('Imag T_{zy}')
    end 
    
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    caxis(u.tiplim)
    hcb = colorbar;
    hcb.Label.String = 'Imag T_{zy}';
    set(gca,'Layer','top')

    annotation('textbox', [0 0.92 1 0.08], ...
    'String', ['Apparent Resistivity and Phase Pseudo-Section. Profile azimuth = ',num2str(azimuth),char(176),' Data rotated to ',num2str(unique(dr.trot)),char(176)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')  
    
    if isfield(d,'niter') && ~isempty(d.niter) % if using S3D, inversion iteration may have been loaded
        folder = ['tipper_pseudo_',d.niter]; % append iteration number to folder name
    else
        folder = 'tipper_pseudo_000';
    end
    if exist(['./',folder],'dir')~=7
        mkdir(folder);
    end
    
    print_figure(folder,['tipper_',d.name,'_rot_',num2str(azimuth)]); %Save figure
%%
    if output_data_menu == 1
        
        % Output real tipper zx to file
        str=[folder,'\pseudo_realtzx_',num2str(azimuth),'.dat'];
        fid1=fopen(str,'w+');
        fprintf(fid1,'               ');
        for is=1:N
            fprintf(fid1,'%9.3f',x_sort(is));   
        end
        fprintf(fid1,'\n');
        fprintf(fid1,'\n');
        for ifreq=1:d.nf
          fprintf(fid1,'%9.3f      ',d.T(ifreq));
            for is=1:N; 
                fprintf(fid1,'%9.3f',real(dr.tip(ifreq,1,index(is)))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        % Output real tipper zy
        str=[folder,'\pseudo_realtzy_',num2str(azimuth),'.dat'];
        fid1=fopen(str,'w+');
        fprintf(fid1,'               ');
        for is=1:N
            fprintf(fid1,'%9.3f',x_sort(is));   
        end
        fprintf(fid1,'\n');
        fprintf(fid1,'\n');
        for ifreq=1:d.nf
          fprintf(fid1,'%9.3f      ',d.T(ifreq));
            for is=1:N; 
                fprintf(fid1,'%9.3f',real(dr.tip(ifreq,2,index(is)))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        %Output imaginary tipper zx
        str=[folder,'\pseudo_imagtzx_',num2str(azimuth),'.dat'];
        fid1=fopen(str,'w+');
        fprintf(fid1,'               ');
        for is=1:N
            fprintf(fid1,'%9.3f',x_sort(is));   
        end
        fprintf(fid1,'\n');
        fprintf(fid1,'\n');
        for ifreq=1:d.nf
          fprintf(fid1,'%9.3f      ',d.T(ifreq));
            for is=1:N; 
                fprintf(fid1,'%9.3f',imag(dr.tip(ifreq,1,index(is)))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        % Output imaginary tipper zy
        str=[folder,'\pseudo_imagtzy_',num2str(azimuth),'.dat'];
        fid1=fopen(str,'w+');
        fprintf(fid1,'               ');
        for is=1:N
            fprintf(fid1,'%9.3f',x_sort(is));   
        end
        fprintf(fid1,'\n');
        fprintf(fid1,'\n');
        for ifreq=1:d.nf
          fprintf(fid1,'%9.3f      ',d.T(ifreq));
            for is=1:N; 
                fprintf(fid1,'%9.3f',imag(dr.tip(ifreq,2,index(is)))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);    
    
    end
%% rotation loop - check user_defaults to enable      
    if u.rotate_data_to_azimuth % if data rotation is enabled, give option to rotate in +/- 5 degree increments
        ichoice = menu(' ','Increase Rotations','Decrease Rotation','Done');

        if ichoice == 1  
            azimuth = azimuth+5; 
            [dr]=rotate_d(dr,5); %Rotate data
        end
        if ichoice == 2  
            azimuth = azimuth-5; 
            [dr]=rotate_d(dr,-5); %Rotate data
        end
        if ichoice == 3;  irun = -1; end
    else
        irun = -1;
    end    
   
end 