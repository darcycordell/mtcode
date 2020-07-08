function plot_rho_pha_pseudo(d)
%
% Function which plots apparent resistivity and phase pseudosections
%
% Usage: plot_rho_pha_pseudo(d)
%
% "d" is an MT data structure
%
% Note: 180 degrees is only added to the YX phase, not YY phase
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

output_data_menu = menu('Output TE and TM mode data to text files?','Yes','No');

irun = 1; % Variable to loop over azimuth

while irun == 1  % Loop over azimuth
%% set up station indices
    c = cosd(azimuth);      s = sind(azimuth);     R = [ c, -s ; s, c]; % this is a positive-counterclockwise rotation matrix
    %Rotate station coordinates to get distance along profile - note
    %stations rotated OPPOSITE of the profile azimuth so that x-coordinate
    %will correspond to distance along the profile. does this only work if
    %profile runs through center of station grid...?
    x_rot = zeros(N,1); y_rot = zeros(size(x_rot));
    for is=1:N
        loc = [x(is),y(is)];    loc = R*loc';
        x_rot(is) = loc(1);     y_rot(is)=loc(2);
    end
    
    % Add loop to put stations in order on monotonically increasing x_rot
    [x_sort,index] = sort(x_rot,'ascend');
    index = sidx(index);
%% plot off-diagonal components
    %========================================================================
    %Plot XY apparent resistivity
    set_figure_size(1);
    subplot(2,2,1)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),log10(squeeze(dr.rho(:,2,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')  %%%%changed
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\rho_{xy} ')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,u.cmap);
    caxis(log10(u.rholim))
    add_rho_colorbar(u);
    set(gca,'Layer','top')
    
    %========================================================================
    %Plot YX apparent resistivity
    subplot(2,2,2)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),log10(squeeze(dr.rho(:,3,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\rho_{yx} ')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,u.cmap);
    caxis(log10(u.rholim))
    add_rho_colorbar(u);
    set(gca,'Layer','top')
      
    %========================================================================
    %Plot XY Phase
    subplot(2,2,3)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(dr.pha(:,2,index)));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\phi_{xy} ')
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,flipud(u.cmap))
    caxis(u.phalim)
    hcb = colorbar;
    hcb.Label.String = 'Phase (degrees)';
    set(gca,'Layer','top')
        
    %========================================================================
    %Plot YX Phase (in 0-90 quadrant)
    subplot(2,2,4)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(dr.pha(:,3,index))+180);
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\phi_{yx} ')
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,flipud(u.cmap))
    caxis(u.phalim)
    hcb = colorbar;
    hcb.Label.String = 'Phase (degrees)';
    set(gca,'Layer','top')

    annotation('textbox', [0 0.92 1 0.08], ...
    'String', ['Apparent Resistivity and Phase pseudo-section. Profile azimuth = ',num2str(azimuth),char(176),' Data rotated to ',num2str(unique(dr.zrot)),char(176)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

    if isfield(d,'niter') && ~isempty(d.niter) % if using S3D, inversion iteration may have been loaded
        folder = ['rho_pha_pseudo_',d.niter]; % append iteration number to folder name
    else
        folder = 'rho_pha_pseudo_000';
    end
    if exist(['./',folder],'dir')~=7
        mkdir(folder);
    end
    
    print_figure(folder,['rho_pha_xy_yx_',d.name,'_rot_',num2str(azimuth)]); %Save figure
%% plot diagonal components  
    %========================================================================
    %Plot XX apparent resistivity
    set_figure_size(2);
    subplot(2,2,1)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),log10(squeeze(dr.rho(:,1,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')  %%%%changed
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\rho_{xx} ')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,u.cmap);
    caxis(log10(u.rholimdiag))
    add_rho_colorbar(u);
    set(gca,'Layer','top')
    
  %========================================================================
    %Plot YY apparent resistivity
    subplot(2,2,2)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),log10(squeeze(dr.rho(:,4,index))));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\rho_{yy} ')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,u.cmap);
    caxis(log10(u.rholimdiag))
    add_rho_colorbar(u);
    set(gca,'Layer','top')
      
    %========================================================================
    %Plot XX Phase
    subplot(2,2,3)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(dr.pha(:,1,index)));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\phi_{xx} ')
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,flipud(u.cmap))
    caxis(u.phalimdiag)
    hcb = colorbar;
    hcb.Label.String = 'Phase (degrees)';
    set(gca,'Layer','top')
        
    %========================================================================
    %Plot YY Phase (NOT moved to 0-90 quadrant, since phase wraps are normal)
    subplot(2,2,4)
    [hp,vp,C] = fix_pcolor(x_sort,log10(d.T),squeeze(dr.pha(:,4,index)));
    pcolor(hp,vp,C)
    hold on; axis ij
    shading 'flat'
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    end
    title('\phi_{yy} ')
    xlabel('Distance (km)')
    ylabel ('Log10 Period (s)')
    axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
    colormap(gca,flipud(u.cmap))
    caxis(u.phalimdiag)
    hcb = colorbar;
    hcb.Label.String = 'Phase (degrees)';
    set(gca,'Layer','top')

    annotation('textbox', [0 0.92 1 0.08], ...
    'String', ['Apparent Resistivity and Phase pseudo-section. Profile azimuth = ',num2str(azimuth),char(176),' Data rotated to ',num2str(unique(dr.zrot)),char(176)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
    
    print_figure(folder,['rho_pha_xx_yy_',d.name,'_rot_',num2str(azimuth)]); %Save figure
    
%%   option to output TE and TM mode data files 
       
    if output_data_menu == 1
    
        % Output TE apparent resistivity to file
        str=[folder,'\pseudo_app_res_TE_mode_',num2str(azimuth),'.dat'];
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
                fprintf(fid1,'%9.3f',dr.rho(ifreq,2,index(is))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        % Output TE phase to file
        str=[folder,'\pseudo_phase_TE_mode_',num2str(azimuth),'.dat'];
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
                fprintf(fid1,'%9.3f',dr.pha(ifreq,2,index(is))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        %Output TM apparent resistivity to file
        str=[folder,'\pseudo_app_res_TM_mode_',num2str(azimuth),'.dat'];
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
                fprintf(fid1,'%9.3f',dr.rho(ifreq,3,index(is))); 
            end
            fprintf(fid1,'\n');
        end
        fclose(fid1);

        % Output TM phase to file
        str=[folder,'\pseudo_phase_TM_mode_',num2str(azimuth),'.dat'];
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
                fprintf(fid1,'%9.3f',dr.pha(ifreq,3,index(is))); 
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
   
end  % End of irun loop