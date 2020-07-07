function plot_polar_distorted(d)
%
% Function which plots polar diagram of MT impedances. Intended for use on synthetic data
% User can apply twist-shear-anis to see effect of galvanic distortion
% EDI file then written out and tested in tensor decomposition
%
% Inputs: d data structure
%
% Outputs: No function outputs, but automatically saves figures and option
% to save as EDI file
%%
close all
u = user_defaults;
rad = 180/pi;   mu = 4*pi*1e-7;   nrot = 72;

%Set up the periods to plot polar diagrams (set in user_defaults or set
%automatically here)
if length(u.tplot)~=3
    disp('Periods to plot distorted polar diagrams is set automatically. See user_defaults.m')
    tlog = rem(log10(d.T)*10,10);
    idx = find(tlog==0);
    mididx = round(length(idx)/2);
    if length(idx)>=3
        tplot = [idx(mididx-1) idx(mididx) idx(mididx+1)];
    else
        tplot = [idx(mididx) idx(mididx) idx(mididx)];
    end
else
    tplot = nearestpoint(u.tplot,d.T);
end
   
d_orig = d;

w = 2*pi./d.T;
i_stn = 1; % Flag to keep looping over stations
is =1; %Station index

while i_stn == 1;

  d.Z = d_orig.Z(:,:,is);   d.Zerr = d_orig.Zerr(:,:,is);  
  d.tip = d_orig.tip(:,:,is);    d.tiperr = d_orig.tiperr(:,:,is); 
%   d.ns = 1;
  d.site{1} = d_orig.site{is}; d.loc = d_orig.loc(is,:);
  
  [d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);
  
  %[polar] = calc_polar(d);
  
  Z_orig = d.Z; % Save undistorted impedances

  irun = 1;   twist = 0;  shear = 0; anis = 0; strike = 0;

  while irun == 1;  % Start loop that allows change in twist, shear, strike etc.
 
    % Distort Z (does all frequencies)
    d.Z = shear_twist_anis(Z_orig,shear,twist,anis);

    % Rotate from strike to measurement co-ordinates (note negative)
    for ifreq = 1:d.nf
        d.Z(ifreq,:) = rotate_Z(d.Z(ifreq,:),-strike);
    end
  
    % This loop rorates through 360 degrees to generate polar diagram of
    % the distorted impedance data 
    rho = nan(d.nf,4,nrot); pha = nan(size(rho)); phaerr = nan(size(rho)); rhoerr = nan(size(rho));
    rho_av_berd = nan(d.nf,nrot); pha_av_berd = nan(d.nf,nrot);
    polar_xy = nan(d.nf,nrot,2); polar_xx = nan(d.nf,nrot,2);
    for irot =1:nrot+1
      rot_ang = (irot-1)*360./nrot;

      % rotate impedances 
      Z_rot = nan(size(d.Z));
      for ifreq = 1:d.nf
        Z_rot(ifreq,:) = rotate_Z(d.Z(ifreq,:),rot_ang);
      end

      c = cosd(rot_ang); s = sind(rot_ang);
      Z_av_berd = nan(d.nf,1)+1i*nan(d.nf,1);
      s1 = nan(d.nf,1); d1 = nan(d.nf,1); s2 = nan(d.nf,1); d2 = nan(d.nf,1);
      swift_skew = nan(d.nf,1); alpha = nan(d.nf,1); bahr_skew = nan(d.nf,1);
      max_polar = nan(d.nf,1);
      for ifreq=1:d.nf
        Z_av_berd(ifreq) =  (Z_rot(ifreq,2)-Z_rot(ifreq,3))/2;
      
        s1(ifreq) = Z_rot(ifreq,1)+Z_rot(ifreq,4);
        d1(ifreq) = Z_rot(ifreq,1)-Z_rot(ifreq,4);
        s2(ifreq) = Z_rot(ifreq,2)+Z_rot(ifreq,3);
        d2(ifreq) = Z_rot(ifreq,2)-Z_rot(ifreq,3);
        swift_skew(ifreq) = abs(s1(ifreq))/abs(d2(ifreq));
        a =2*real(s2(ifreq)*conj(d1(ifreq)))/(abs(d1(ifreq))*abs(d1(ifreq))-abs(s2(ifreq))*abs(s2(ifreq)));
        alpha(ifreq) = rot_ang+rad*0.25*atan(a);
    
        c1 = real(d1(ifreq))*imag(s2(ifreq))-real(s2(ifreq))*imag(d1(ifreq));
        c2 = real(s1(ifreq))*imag(d2(ifreq))-real(d2(ifreq))*imag(s1(ifreq));
        bahr_skew(ifreq)=sqrt(abs(c1-c2))/abs(d2(ifreq));
        
        [rho(:,:,irot), pha(:,:,irot), rhoerr(:,:,irot), phaerr(:,:,irot)] = calc_rho_pha(Z_rot,d.Zerr,d.T);
      
        rho_av_berd(ifreq,irot) = abs(Z_av_berd(ifreq)*Z_av_berd(ifreq))/(w(ifreq)*mu);
        pha_av_berd(ifreq,irot)= 90-rad*atan2(real(Z_av_berd(ifreq)),imag(Z_av_berd(ifreq)));
      
        polar_xy(ifreq,irot,1) = c*abs(Z_rot(ifreq,2));   polar_xy(ifreq,irot,2) = s*abs(Z_rot(ifreq,2));
        polar_xx(ifreq,irot,1) = c*abs(Z_rot(ifreq,1));   polar_xx(ifreq,irot,2) = s*abs(Z_rot(ifreq,1));
     
        % Find maximum amplitude of polar ellipse
        m1 = max(abs(polar_xy(ifreq,:,1)));  m2 = max(abs(polar_xy(ifreq,:,2)));
        max_polar(ifreq) = 1.2*max([m1,m2]);
      
      end
    end

% Now plot results.  Note x is North
    set_figure_size(1);

    rot_ang = 0;
  
    subplot(3,2,1)   % PLOT APPARENT RESISTIVITY
    loglog(d.T,rho(:,2,irot),'r.') ; hold on;
    loglog(d.T,rho(:,3,irot),'b.') ;
    loglog(d.T,rho_av_berd(:,irot),'g-') ;
    
    loglog(d.T(tplot(1)),rho(tplot(1),2,irot),'ro') ;
    loglog(d.T(tplot(2)),rho(tplot(2),2,irot),'ro') ;
    loglog(d.T(tplot(3)),rho(tplot(3),2,irot),'ro') ;

    loglog(d.T(tplot(1)),rho(tplot(1),3,irot),'bo') ;
    loglog(d.T(tplot(2)),rho(tplot(2),3,irot),'bo') ;
    loglog(d.T(tplot(3)),rho(tplot(3),3,irot),'bo') ;
    
    loglog(d.T,rho(:,1,irot),'r:') ;   loglog(d.T,rho(:,4,irot),'b:') ; 
    
    %DC Not sure what this is....
%     for i1 = 1:length(r_anot); loglog([tlim],[r_anot(i1),r_anot(i1)],'k:');  end
%     for i1 = 1:length(t_anot); loglog([t_anot(i1),t_anot(i1)],[rlim],'k:');  end
%       
    hold off
    
    ylabel ('\rho_a (\Omega m)');
    axis([u.Tlim u.rholim])
    
    tstr = [char(d.site(is)),'  \theta = ', num2str(rot_ang),'^o','; Twist =', num2str(twist),'^o','; Shear = ',num2str(shear),'^o',';Anis =',num2str(anis),';Strike =',num2str(strike)];
    title(strrep(tstr,'_','\_'));   % Force MATLAB to plot underscores

    subplot(3,2,3)  % PHASE PLOT
    semilogx(d.T,pha(:,2,irot),'r.') ; hold on;
    %semilogx(d.T,pha(:,1,irot),'r:') ;
    semilogx(d.T,pha(:,3,irot)+180,'b.') ;
    %semilogx(d.T,pha(:,4,irot),'b:') ;
    semilogx(d.T,pha_av_berd(:,irot),'g-') ;
    
    semilogx(d.T(tplot(1)),pha(tplot(1),2,irot),'ro') ; semilogx(d.T(tplot(2)),pha(tplot(2),2,irot),'ro') ; semilogx(d.T(tplot(3)),pha(tplot(3),2,irot),'ro') ;
    semilogx(d.T(tplot(1)),pha(tplot(1),3,irot)+180,'bo') ; semilogx(d.T(tplot(2)),pha(tplot(2),3,irot)+180,'bo') ; semilogx(d.T(tplot(3)),pha(tplot(3),3,irot)+180,'bo') ;
      
    %DC Not sure what this is...
%     for i1 = 1:length(p_anot); loglog([tlim],[p_anot(i1),p_anot(i1)],'k:');  end
%     for i1 = 1:length(t_anot); loglog([t_anot(i1),t_anot(i1)],[plim],'k:');  end
%    
    axis([u.Tlim u.phalim]); hold off;
    ylabel ('Phase (deg)');
 
    
    subplot(3,2,5)   % PLOT SKEW
    semilogx(d.T,swift_skew,'k.-') ; hold on;
    semilogx(d.T,bahr_skew,'g.-') ;
    axis([u.Tlim -0.05 0.5])
    xlabel ('T(s)'); ylabel ('Skew');
    title('Black = Swift; Green = Bahr')
    hold off;
    
    % plot polar diagrams
    frame_plot = [3,6,9];
    for iplot =1:length(tplot)
    
      subplot(3,3,frame_plot(iplot))
      
      % Plot whole ellipse
      plot(polar_xy(tplot(iplot),:,2),polar_xy(tplot(iplot),:,1),'k-') % Off-diagonal
             hold on
      plot(polar_xx(tplot(iplot),:,2),polar_xx(tplot(iplot),:,1),'k:') % Diagonal
      
      % Now plot axis showing co-ordinate frame
      r = max_polar(tplot(iplot));
      xc1 = [0, r*cos(rot_ang*pi/180.)];      yc1 = [0, r*sin(rot_ang*pi/180.)] ;  
      xc2 = [0, r*cos(pi/2+rot_ang*pi/180.)]; yc2 = [0, r*sin(pi/2+rot_ang*pi/180.)];   
      plot(yc1,xc1,'r-'); plot(yc2,xc2,'b-');
      if ~isnan(r)
        axis([-r r -r r])
        title(['T = ',num2str(d.T(tplot(iplot))),' s'])
      else
          title(['No Data at T = ',num2str(d.T(tplot(iplot))),' s'])
      end
           
      axis equal
      
      hold off            
    end
    
    d.zrot(d.zrot==0)=strike;
    
    fname_edi = [char(d.site{is}),'-distort-st=',num2str(strike,'%03.0f'),'tw=',num2str(twist,'%02.0f'),'sh=',num2str(shear,  '%02.0f'),'an=',num2str(anis,  '%04.2f'),'.edi'];
  
    print_figure('polar_distort',[fname_edi(1:end-4),'_polar_distorted']);
  
    %Change Distortion
    prompt={'Twist Angle (e.g. 15)','Shear Angle (e.g. 30)','Anisotropy (e.g. 0.1)','Strike Direction (E of N)'};
    dlg_title='Distort Data';
    def={num2str(twist),num2str(shear),num2str(anis),num2str(strike)};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(dinp)
        irun = -1;
    else
        twist = str2double(dinp{1});
        shear = str2double(dinp{2});
        anis = str2double(dinp{3});
        strike = str2double(dinp{4});
    end
  

  end  % End of main loop (while ...)
% 
% Add Gaussian noise (0.01 = 1%)


    end_menu = 1;
    while end_menu
      ichoice = menu('  ','Next station','Previous station','Add Noise','Save Data to EDI','Exit');
      if ichoice == 1;   is=is+1; end_menu=0; end
      if ichoice == 2;   is=is-1; end_menu = 0; end
      
      if is < 1;  is =1; disp('Station #1, no more stations'); end                  
      if is > d.ns; is =d.ns; disp('Final Station reached, no more stations'); end

      if ichoice == 3
          d = add_noise(d,0);
      end

      if ichoice == 4;
          d.site = {fname_edi(1:end-4)};
          write_data_edi(d,1);
      end
      
      if ichoice == 5;   i_stn = -1; end_menu = 0; end
    end

end

% Output period, rho,pha to file
%str=[char(stn_name),'_rho_pha.dat'];
%fid1=fopen(str,'w+');
%fprintf(fid1,'%9.3f %9.3f %9.3f %9.3f %9.3f \n',[T;rho_xy;rho_yx;pha_xy;pha_yx]);
%fclose(fid1);

end