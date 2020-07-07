function OCCAM1D
% This code solves the OCCAM 1D Solution for the MT problem. This was
% primarily written by Ersan but re-worked and re-written by Darcy with
% comments. Ersan's GUI was also removed to make for an easier and more
% usable code
%
% User is able to load EDI data and choose to invert either the TE, TM or
% determinant impedance
%
% User is also able to load a text file containing columns of frequency,
% apparent resistivity, phase, apparent resistivity error and phase error
%
% Finally, user is able to make synthetic data within the program and
% invert it.
%
% For more info on OCCAM inverison see:
%   Constable, S., Parker, R. L., & Constable, C. (1987). Occam’s inversion: 
%   A practical algorithm for generating smooth models from electromagnetic 
%   sounding data. Geophysics, 52(3), 289–300. 
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

clear all; close all
curdir = pwd;
% Plotting options
reslims = [1 1000]; %resistivity limits 
depthlims = [0 20000]; %depth limits for plotting

get_data = '(1) Get Data';
mesh_design = '(2) Mesh Design';
inv_params = '(3) Inversion Parameters';
inv_run = '(4) Run Inversion';

track = 0; count = 1;
while 1

main_menu = menu('',get_data,mesh_design,inv_params,inv_run,'Plot Results','Save Results','Quit');


if main_menu == 1
    track(count) = main_menu;
    get_data = '(1) Get Data (DONE)';
    % DATA INPUT OPTIONS  ----------------------------------------------------
    data_menu = menu('','Load EDI','Load Text File','Make Synthetic Data (default)');

    if data_menu == 1 % DATA FROM EDI
        [edifile,edipath] = uigetfile({'*.edi'},'Pick EDI file');

        if edifile == 0
            return
        end
            
        cd(edipath)
        DS = load_data_edi(edifile);
        lat = DS.loc(1);
        lon = DS.loc(2);
        elev = DS.loc(3);
        rot = median(DS.zrot);
        f = DS.f; %Frequencies (Hz)
        T = 1./f; %Periods (s)
        nd = length(f);
        cd(curdir)

        mode_menu = menu('Type of 1-D Data to Use:','XY mode','YX mode','Determinant Average (default)');

        if mode_menu == 1 %Fit the TE mode
            datatype = 'XY';
            rhoa = DS.rho(:,2)'; %apparent resistivity data
            Z = DS.Z(:,2); %impedance
            dZ = real(DS.Zerr(:,2)); %impedance error
            pha = DS.pha(:,2)'; %phase data
        elseif mode_menu == 2 %Fit the TM mode
            datatype = 'YX';
            rhoa = DS.rho(:,3)'; %swap modes
            pha = DS.pha(:,3)'; %swap modes
            %Depending on processing and EDI formatting, it is sometimes necessary
            %to add +180 to the TM EDI phases
            if mean(pha)>0 && mean(pha)<90
                %do nothing
            else
                pha = pha+180;
            end
            Z = DS.Z(:,3)*exp(1i*-pi); %Rotates yx impedance by 180 degree to make the TM mode "look" like the 1-D TE mode
            dZ = real(DS.Zerr(:,3)); %Not sure if I should rotate TM errors...
        else %Fit the determinant impedance
            datatype = 'DETERMINANT';
            %Determinant impedance as defined by Rung-Arunwan et al., 2017
            Z = (sqrt((DS.Z(:,1).^2+DS.Z(:,2).^2+DS.Z(:,3).^2+DS.Z(:,4).^2)./2));

            %The error is propagated using standard error propagation rules
            % A^2 --> dA^2
            % sqrt(A) --> sqrt(dA)
            % a+b+c --> sqrt(da^2 + db^2 + dc^2)
            dZ = real(sqrt(0.5*sqrt(DS.Zerr(:,1).^4+DS.Zerr(:,2).^4+DS.Zerr(:,3).^4+DS.Zerr(:,4).^4)));

            rhoa = ((1./(2*pi*DS.f'*4*pi*10^-7)).*abs(Z.').^2); %Apparent resistivity formula
            pha = (atan2(imag(Z),real(Z))*(180/pi)).'; %Phase formula

        end

        dZ_orig = dZ;

        %Set errors based on error floors. If the measured value (e.g. apparent
        %resistivity) is less than err_flr*value, then the error is boosted to that
        %level. This method makes error floors larger for small errors but leaves
        %large errors unchanged.
        if mode_menu == 1 || mode_menu == 2
            floorZ=abs(Z);
        else
            floorZ = sqrt(abs(DS.Z(:,2)).*abs(DS.Z(:,3)));
        end
        
    elseif data_menu == 2 %DATA FROM TEXT FILE WITH COLUMNS: f, rhoa, pha, rhoerr, phaerr
        [datfile,datpath] = uigetfile({'*'},'Pick text file');
        
        if datfile == 0
            return
        end
        
        cd(datpath)
        [~,~,ext] = fileparts(datfile);
        
        if strcmp(ext,'.dat') || strcmp(ext,'.txt') || strcmp(ext,'.csv') || strcmp(ext,'.xls')
            d = table2array(readtable(datfile));
        elseif strcmp(ext,'.resp')       
            fid = fopen(datfile);
            d = fscanf(fid,'%e %e %e %e %e',[5 Inf])';
            fclose(fid);
        else
            error('Error: File type not supported. Must be .dat, .txt, .csv, .xls, or .resp')
        end
        
        cd(curdir)
        
        lat = 0;
        lon = 0;
        elev = 0;
        rot = 0;
        
        datatype = 'TXT';
        f = d(:,1); %Frequencies (Hz)
        T = 1./f; %Periods (s)
        nd = length(f);
        rhoa = d(:,2)'; %apparent resistivity data
        pha = d(:,3)'; %phase data
        rhoerr = d(:,4); %apparent resistivity error
        phaerr = d(:,5); %phase error
        
        [Z,dZ] = calc_Z(rhoa',rhoerr,pha',phaerr,T);
        
        floorZ=abs(Z);

    else % SYNTHETIC DATA
        datatype = 'SYNTHETIC';
        prompt={'Number of Frequencies','Minimum Frequency','Maximum Frequency','Gaussian Error to Add','Depth to Top of Model Layers','Layer Resistivites'};
        dlg_title='Synthetic Data Parameters';
        def={num2str(80),num2str(0.001),num2str(1000),num2str(0.05),'[0 150 250 5000 5500]','[100 1 50 1 10]'};
        %def={num2str(50),num2str(0.001),num2str(1000),num2str(0.1),'[0 2000 10000]','[100 1 1000]'};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);
        
        if isempty(dinp)
            return
        end

        num_freq = str2double(dinp{1});
        min_freq = str2double(dinp{2});
        max_freq = str2double(dinp{3});
        err = str2double(dinp{4}); %Gaussian error to add to data
        model_depth = str2num(dinp{5}); %Depth to top of each model layer
        model_res = str2num(dinp{6}); %Resistivity of each layer

        % FROM MT1D SYNTHETIC
        f = 10.^linspace(log10(min_freq),log10(max_freq),num_freq).';
        nd = length(f);
        T = 1./f;
        lat = 0;
        lon = 0;
        elev = 0;
        rot = 0;

        [fwd]=calc_fwd_1d(model_depth,model_res,f',err);

        Z = fwd.Z.';
        rhoa = fwd.rho;
        pha = fwd.phi;

        dZ = err*abs(Z);
        floorZ = abs(Z);
        

    end

    %Set constants and frequencies
    w=2*pi*f';mua=4*pi*10^-7; iwm=1i*mua*w;  
    ndR=2*nd;
end

if main_menu == 2 && max(track) >= 1
    track(count) = main_menu;
    mesh_design = '(2) Mesh Design (DONE)';
    % MESH DESIGN ------------------------------------------------------------
    %The mesh can be designed manually using skin depth, but the way the code
    %is currently implemented has a mesh generated automatically using Bostick
    %inversion. (Bostick is automated and works well)
    mesh_menu = menu('Mesh Design Using:','Bostick Automated (default)','Manually Set','Load Text File');

    if mesh_menu == 2 %Manually Set Thicknesses
        meshtype = 'Manually Entered';
        prompt={'Number of Layers','First Thickness in meters (default is 1/3 min skin depth)','Maximum Depth in meters (default is 3*(max skin depth))','Geometric Factor'};
        dlg_title='Mesh Parameters';
        %def={num2str(80),num2str(503*sqrt(min(T)*100)/3),num2str(503*sqrt(max(T)*100)*3),num2str(1.2)};
        def={num2str(137),num2str(10),num2str(200000),num2str(1.1)};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        nl = str2double(dinp{1});
        first_layer = str2double(dinp{2});
        maxdepth = str2double(dinp{3});
        geo_factor = str2double(dinp{4});

        thick=[]; thick(1) = first_layer; i =2;
        while sum(thick)<=maxdepth
            thick(i)=thick(i-1)*geo_factor;
            i=i+1;
        end

    elseif mesh_menu == 3
        [modfile,modpath] = uigetfile({'*'},'Pick Text File with Model Thicknesses');

        meshtype = modfile;
        cd(modpath)
        thick = load(modfile)';
        cd(curdir)

    else
        meshtype = 'Bostick';
        %--------------thickness guess using Bostick-----------------
        Dinv=(rhoa./(mua*w)).^0.5;
        roaD=abs(rhoa.*(180./(2*pha)-1));

        if Dinv(1)>Dinv(2)
            D=Dinv(length(Dinv):-1:1);
        else 
            D=Dinv;
        end
        thick = [D(1)/4 D(1)/2 D(1:end-1)]; %Thicknesses (from Bostick)
    end

    nl = length(thick)+1; %Number of layers
    depth=[0 cumsum(thick)]; %Depth of layer tops
end

if main_menu == 3 && max(track) >= 2
    % INVERSION PARAMETERS ---------------------------------------------------

    prompt={'Max Iterations','Error Floor','Desired RMS','Starting Halfspace Resistivity','Set Discontinuity (1 or 0)'};
    dlg_title='Inversion Parameters';
    def={num2str(200),num2str(0.1),num2str(0.5),num2str(100),num2str(0)};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(dinp)
        track(count) = main_menu;
        inv_params = '(3) Inversion Parameters (DONE)';
        itmax = str2double(dinp{1}); %Maximum number of iterations (usually 10 - 30 is ok)
        err_flr = str2double(dinp{2}); %The data error floor (usually 0.01 to 0.1 is ok). This is a true error floor as a direct
                %percentage of the apparent resistivity and phase, independently.
                %There is no scaling or error propagation carried out.
                %The choice of error floor has a significant impact on the
                %inversion fit and needs to be tuned accordingly.

        %taumax = 100; %Maximum conductance (for plotting purposes)
        rmstreshold = str2double(dinp{3}); %RMS misfit threshold (usually 1 is good for synthetic)
        inres = str2double(dinp{4}); %Starting model resistivity in Ohm m.
        discontinuity = str2double(dinp{5});
        
        if discontinuity
            prompt={'Depth of Disontinuity or Discontinuities (km)'};
            dlg_title='Discontinuity Parameters';
            def={num2str([1 5])};
            num_lines=1;
            disc_inp = inputdlg(prompt,dlg_title,num_lines,def);
            
            if isempty(disc_inp)
                discontinuity = 0;
            else
                disc_depth = str2num(disc_inp{1});
                for i = 1:length(disc_depth)
                    disc_dist = abs(depth - disc_depth(i)*1000);
                    [~, disc_idx(i)] = min(disc_dist);
                end
            end
                    
        end
        
    end

end

if main_menu == 4 && max(track) >= 3
    track(count) = main_menu;
    inv_run = '(4) Run Inversion (DONE)';
    % OCCAM INVERSION----------------------------------------------------------

    %Apply error floor-------------------------------
    [dZ] = apply_errorfloor(dZ,err_flr,floorZ);

    rhoerr = (2*rhoa'.*dZ)./abs(Z);
    phaerr = (180/pi)*(dZ)./abs(Z);

    W = diag(1./[dZ; dZ]); %Data error weighting matrix used in inversion
    d=[real(Z);imag(Z)]'; %Data vector used in inversion

    %-------------starting model-------------------
    clear m J F
    m(1,:)=ones(1,nl)*inres;

    %FORWARD CALCULATION ---------------------------------------------
    %"F" is the forward of real and imaginary impedance values (across columns). Each
    %inversion iteration is saved in a new row.
    [F(1,:),ra,ph]=mt_fwd_occ(m(1,:),nl,w,thick); 
    %-----------------------------------------------------------------

    %Plot initial data and starting model--------------------------------------
    set_figure_size(1);
    subplot(2,2,1)
    logerrorbar(T,rhoa',rhoerr,'.k','-k'); hold on
    loglog(T,ra','--k');
    xlabel('Period (s)')
    ylabel('App Res (\Omega m)')
    axis([min(T) max(T) reslims])
    grid on

    subplot(2,2,3)
    errorbar(T,pha',phaerr,'.k'); hold on
    semilogx(T,ph','--k')
    set(gca,'XScale','log')
    xlabel('Period (s)')
    ylabel('Phase (deg)')
    axis([min(T) max(T) 0 90])
    grid on

    subplot(2,2,[2 4])
    stairs(m(1,:),depth/1000,'-k')

    %--------------------------------------------------------------------------

    %Begin inversion
    X2 = []; rms = []; R1= []; D1 = [];
    rms(1) = norm(W*(F(1,:))'-d')/sqrt(2*nd); %Initial rms of starting model
    rgh1=diag([0;ones(nl-1,1)]) + diag(-ones(nl-1,1),-1);         %delta matrix
    
    if discontinuity
        rgh1(disc_idx,:) = 0;
    end
    
    rgh2=rgh1'*rgh1;                                              %delta' x delta  
    mumax = []; mumax(1)=10000;
    iter = 2; 

    %Iterate to solve the inversion problem using Occam algorithm
    for iter=2:itmax

        if rms(iter-1)>2
            rmsdes=rms(iter-1)/1.5;
        else
            rmsdes=rmstreshold;
        end

        %Build the Jacobian matrix --------------------------------------------
        %I am not sure why Ersan chose to do this in logarithmic space but it
        %seems strange to me that he takes logarithms of both apparent
        %resistivity AND phase (since phase is a linear quantity). I modified
        %the code so that it works by inverted real and imaginary impedance
        %data instead of rho and pha.
        parameter = log10(m(iter-1,:));
        apt = (F(iter-1,:));
        for I=1:length(parameter);
            parameter(I)=parameter(I)+0.005;
            [cpt,~,~] = mt_fwd_occ(10.^parameter,nl,w,thick);      
            turev=(cpt-apt)/0.005;
            J(:,I)=turev';
            parameter(I)=parameter(I)-0.005;
        end

        %----------------------------------------------------------------------

        %Inverse algorithm (from Constable et al., 1987)
        son=(W*J)'*W*J;
        b=(W*J)'*W*(d-(F(iter-1,:))+(J*log10(m(iter-1,:))')')';
        %The solution is found using a "golden section search" algorithm to
        %find the minimum
        [mumax(iter), rms(iter), m(iter,:), F(iter,:),x3,f4]=golden_section(son,b,W,rgh2,d,rmsdes,0.0001,1,mumax(1), 0.01,2*nd,w,thick,nl);
        X3(:,iter-1)=x3';  F4(:,iter-1)=f4';
        R1(iter)=((rgh1*log10(m(iter,:))')'*rgh1*log10(m(iter,:))');                 
        DD(iter)=(log10(m(iter,:))-log10(m(iter-1,:)))*(log10(m(iter,:))-log10(m(iter-1,:)))';   
        if DD(iter)<0.01
            break
        end

        if abs(rms(iter)-rms(iter-1)) < 0.01
            break
        end

    end

    X2=rms.^2*nd;
    [~,ra_mod,ph_mod]=mt_fwd_occ(m(end,:),nl,w,thick); 

    %PLOT RESULTS--------------------------------------------------------------
    figure(1)
    subplot(2,2,1)
    h(1) = loglog(T,ra_mod,'-m','LineWidth',1);
    title(['Total Iterations = ',num2str(iter),'. Final Occam RMS = ',num2str(rms(end))])

    figure(1)
    subplot(2,2,3)
    h(2) = semilogx(T,ph_mod,'-m','LineWidth',1);

    figure(1)
    subplot(2,2,[2 4])
    h(3) = stairs(m(iter,:),depth/1000,'-m','LineWidth',1);
    set(gca,'XScale','log')
    axis([reslims depthlims/1000])
    xlabel('Resistivity (\Omega m)')
    ylabel('Depth (km)')
    axis ij
    grid on
    
end

if main_menu == 5 && max(track) >= 4
    track(count) = main_menu;
    %PLOT ALL ITERATIONS--------------------------------------------------
    im = 1;
    while 1

        [~,ra_view,ph_view]=mt_fwd_occ(m(im,:),nl,w,thick); 

        h(1).YData = ra_view;
        subplot(2,2,1)
        title(['Total Iterations = ',num2str(im),'. Final Occam RMS = ',num2str(rms(im))])
        h(2).YData = ph_view;

        figure(1)
        subplot(2,2,[2 4])
        stairs(m(im,:),depth/1000,'-m','LineWidth',1); hold off
        set(gca,'XScale','log')
        axis([reslims depthlims/1000])
        xlabel('Resistivity (\Omega m)')
        ylabel('Depth (km)')
        axis ij
        grid on

        next_menu = menu('','Next Iteration','Quit');

        if next_menu == 1

            im = im+1;
            if im > iter
                break
            end

        else
            break
        end

    end

end

if main_menu == 6 && max(track) >= 4
    track(count) = main_menu;
    %SAVE RESULTS
    prompt={'Filename to Save'};
    dlg_title='Filename to Save';
    def={'filename_occam1d'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    filename = dinp{1};
    
    rms_mod = rms(end);
    model = m(end,:)';
    starting_model = m(1,:)';
    save([filename,'.mat'],'T','ra_mod','ph_mod','rms_mod','model','thick','depth','iter','lat','lon','elev','rot','datatype','meshtype','err_flr','starting_model');
    save([filename,'_data_input.mat'],'rhoerr', 'phaerr', 'rhoa', 'pha', 'Z', 'dZ');
    
    fid = fopen([filename,'_ef_',num2str(err_flr*100),'.mod'],'w');
    A = [thick' model(1:end-1)];
    fprintf(fid,'%12.3f %e\r\n',A');
    fclose(fid);
    
    
    fid = fopen([filename,'_ef_',num2str(err_flr*100),'.resp'],'w');
    write_error = 10^-9*ones(1,nd);
    A = [f ra_mod' ph_mod' write_error' write_error'];
    fprintf(fid,'%e %e %e %e %e\r\n',A');
    fclose(fid);
    
end

if main_menu > 6
    break
end

count = count+1;
    
end
    
    
    





