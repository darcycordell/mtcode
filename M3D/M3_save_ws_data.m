function M3_save_ws_data(hObject, ~, ~)
    
H=guidata(hObject);

if ~isfield(H,'dat') % check if data has been loaded
    
    disp('No data to save! Load data first.');
    warndlg('Enter data parameters first.')
    
else % data has been loaded
    
    if ~isfield(H,'max_per')||~isfield(H,'min_per')||~isfield(H,'per_skip')
        warndlg('Error: enter values for min per, max per, or periods to skip.')
        return
    end
    
    % give variables shorter names
    H.dex = H.d; % give temp name to d struct that will be exported
    %H.dex = rmfield(H.dex,{'zrot','trot','rho','pha','rhoerr','phaerr'});
    ZZZ = permute(H.d.Z,[2 1 3]); % to handle new d struct
    EEE = permute(H.d.Zerr,[2 1 3]);
    ZZZ_HZ = permute(H.d.tip,[2 1 3]);
    EEE_HZ = permute(H.d.tiperr,[2 1 3]);
    site = H.d.site;
    T = H.d.T;
    coord = [H.d.loc(:,2) H.d.loc(:,1) H.d.loc(:,3)];
    ns = length(site);    
    max_per = str2double(get(H.max_per,'string'));
    min_per = str2double(get(H.min_per,'string'));
    per_skip = str2double(get(H.per_skip,'string'));
    err_flr_Zd = .01*str2double(get(H.errflr_Zdiag,'string')); % diagonal Z
    err_flr_Zod = .01*str2double(get(H.errflr_Zodiag,'string')); % off-diagonal Z
    err_flr_tip = str2double(get(H.errflr_tip,'string')); % tipper
    err_flr_type = str2double(get(H.errflr_type,'string')); %error denominator type (value of 1 or 2. Default = 1)
    resp = get(H.resp,'value');
%     H.rot_cor_ang = H.mesh_rot;
      
    if get(H.conj_Z,'value')
        ZZZ = conj(ZZZ);
    end
    if get(H.conj_tip,'value') % GN - added ability to take conjugate of tipper (April 2012)
        ZZZ_HZ=conj(ZZZ_HZ);
    end
    
    if H.mesh_rot ~= 0 % Rotate impedance data to opposite direction of the mesh if needed
        warndlg(['MT data rotated ',num2str(H.mesh_rot),' degrees.'])
%         c = cos(H.rot_cor_ang*pi/180.);
%         s = sin(H.rot_cor_ang*pi/180.);
%         U = [ c*c ,  c*s ,  c*s , s*s  ; -c*s ,  c*c , -s*s , c*s  ; -c*s , -s*s ,  c*c , c*s  ; s*s , -c*s , -c*s , c*c ]; % Equation 5.4 in Simpson and Bahr
%         Uv = [ c  , s ; -s  , c ];% Regular rotation matrix
%         for tpr=1:length(site)
%             ZZZ(:,:,tpr)=U*squeeze(ZZZ(:,:,tpr));
%             EEE(:,:,tpr)=(U.^2*squeeze(EEE(:,:,tpr)).^2).^.5;
%             ZZZ_HZ(:,:,tpr)=Uv*squeeze(ZZZ_HZ(:,:,tpr));
%             EEE_HZ(:,:,tpr)=(Uv.^2*squeeze(EEE_HZ(:,:,tpr)).^2).^.5;
%         end
%         warndlg(['MT impedances rotated ', num2str(H.rot_cor_ang),' degrees NE! Make sure when you export edis from S3DET_v5 you keep this azimuthal information in the edi files'])
%         %winglink assumes mesh rotation is negative in NE and positive
%         %in NW. Here program understands this convention but shows
%         %always positive angles. If your winglink exported model has
%         %-20 rotation, this program will correctly interpret the angle but will show you as positive
    end
    
    EEE=sqrt(EEE)/796; % sqrt is necessary to convert to std dev
    ZZZ=ZZZ./796; % conversion from field units... necessary for field data
    
    % Impose error floor of the type specified by user:
    % Edits by DC 18/11/14
        if err_flr_type==1 % Denominator is sqrt(abs(Zxy*Zyx)) for ALL Z components
        
            [EEE,EEE_HZ,ZZZ,calc_err,calc_err_edi]=M3_err_flr_type_1(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod,err_flr_tip);
        
        elseif err_flr_type==2 % Denominator is abs(Zxy) for Zxx,Zxy. Denominator is abs(Zyx) for Zyx,Zyy
        
            [EEE,EEE_HZ,ZZZ,calc_err,calc_err_edi]=M3_err_flr_type_2(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod,err_flr_tip);
            
        elseif err_flr_type==3 % The OLD WAY: Original code untouched
            nfim=length(T);
            % DC EDIT Jan 21, 2015 for plotting purposes
            calc_err_edi=[];calc_err=[];
            for ic=1:4
                if ic==1 || ic==4
                    calc_err_edi(ic,:,:)=(EEE(ic,:,:))./sqrt(abs((ZZZ(2,:,:).*(ZZZ(3,:,:)))));
                else
                    calc_err_edi(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(ic,:,:));
                end
            end
            %end DC edits Jan 21, 2015
            for ic=1:4
                for itrpo=1:length(H.d.x)
                    for iu=1:nfim   %number of frequencies in *mat file
                        if ic==1 || ic==4
                            if EEE(ic,iu,itrpo)<err_flr_Zd*sqrt((ZZZ(2,iu,itrpo).*conj((ZZZ(3,iu,itrpo))))) %Conjugate is irrelevant. I think he should be comparing modulus --- DC OCT 28/14
                                EEE(ic,iu,itrpo)=err_flr_Zd*sqrt(abs(ZZZ(ic,iu,itrpo).*ZZZ(ic,iu,itrpo)));   %10 errorfloor
                            end
                        elseif ic==2 || ic==3
                            EEE(ic,iu,itrpo)=err_flr_Zod*sqrt(abs(ZZZ(ic,iu,itrpo).*ZZZ(ic,iu,itrpo)));   %10 errorfloor
                        end
                        if (ic ==1 || ic == 2) && (resp == 2 || resp == 4) %tipper
                            if EEE_HZ(ic,iu,itrpo) < err_flr_tip
                                EEE_HZ(ic,iu,itrpo)=err_flr_tip;   %absolute error floor for tipper error
                            end
                        end
                    end
                end
            end 
    
                %DC EDITS Jan 21, 2015 for plotting purposes
            for ic=1:4
                if ic==1 || ic==4
                    calc_err(ic,:,:)=(EEE(ic,:,:))./sqrt(abs((ZZZ(2,:,:).*(ZZZ(3,:,:)))));
                else
                    calc_err(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(ic,:,:));
                end
            end
            
            %Option to plot error percent to observe error floor being applied
            %plot_err_compare(T,calc_err,calc_err_edi,err_flr,err_flr2,site,H,err_flr_type)
            
       %End DC EDITS Jan 21, 2015 
        
        else %Displays error: err_flr_type must be 1 or 2. Defaults to err_flr_type=1
            fprintf('Error Floor Type must be 1, 2 or 3. \nYour entry was invalid: Default error floor is Type 1\nSee M3DET documentation for more details')
            [EEE,EEE_HZ,ZZZ,calc_err,calc_err_edi]=M3_err_flr_type_1(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod);   
        end
        
        %End Edits by DC 18/11/14
        
    minid=max(flipud(find(min_per>T))) + 1; % adding 1 ensures that all periods are greater than per_min
    maxid=min(find(max_per<T)) - 1; % subtracting 1 ensures that all periods are less than per_max
    if isempty(minid)==1
        minid =1;
    end

    if isempty(maxid)==1
        maxid=length(T);
    end
    frtp = minid:per_skip:maxid; % frequencies to plot    
    set(H.data_text,'string',sprintf('%s\n%s\n%s',['Number of periods = ',num2str(numel(frtp))],['Minimum: ',num2str(T(min(frtp))),' s'],['Maximum: ',num2str(T(max(frtp))),' s']))
                           
    if exist('frtp','var')
            T = T(frtp);         %updating periods
            
            Z = permute(ZZZ(:,frtp,:),[3 1 2]); % re-arrange the impedance array into ns x nr x np
            Z = [real(Z(:,1,:)) imag(Z(:,1,:)) real(Z(:,2,:)) imag(Z(:,2,:)) real(Z(:,3,:)) imag(Z(:,3,:))  real(Z(:,4,:)) imag(Z(:,4,:))];
            Znan = isnan(Z);
            Zzero = (Z==0);
            Z(Znan) = 4.6326e-005; % NaN in impedance set to this value
            Z(Zzero) = -1.8452e-005; % zeros in impedance set to this value
            
            EZ = permute(EEE(:,frtp,:),[3 1 2]); % re-arrange error array  
            EZ = [EZ(:,1,:) EZ(:,1,:) EZ(:,2,:) EZ(:,2,:) EZ(:,3,:) EZ(:,3,:) EZ(:,4,:) EZ(:,4,:)];
            EZ(Znan) = 0.1; % NaN in impedance
            EZ(Zzero) = 0.1; % zeros in impedance
            
            ERZ = ones(size(Z)); % errormap
            ERZ(Znan) = 999; % NaN in impedance
            ERZ(Zzero) = 999; % zeros in impedance
            
            if exist('ZZZ_HZ','var')
                Tip = permute(ZZZ_HZ(:,frtp,:),[3 1 2]);
                Tip = [real(Tip(:,1,:)) imag(Tip(:,1,:)) real(Tip(:,2,:)) imag(Tip(:,2,:))];
                Tipnan = isnan(Tip);
                Tip(Tipnan) = 0;

                ETip = permute(EEE_HZ(:,frtp,:),[3 1 2]);
                ETip = [ETip(:,1,:) ETip(:,1,:) ETip(:,2,:) ETip(:,2,:)];
                ETip(Tipnan) = 10; % find NaN in the tipper- make the errors 10

                ERTip = ones(size(Tip)); % errormap
                ERTip(Tipnan) = 999; % NaN in tipper
            end

            %------------write data---------------------------------------
            def = {'data_name'};
            prompt = {'enter data file name must be same as model name'};
            titles  = 'Data file name';
            datf = char(inputdlg(prompt,titles,1,def)); %model file name
            datf=[datf,'.data'];
            while exist(datf)>0
                newdata_name = questdlg('This file exist! Do you want to overwrite the file?', ...
                    'The data file exist', ...
                    'Yes', 'No','Yes');
                switch newdata_name
                    case 'Yes'
                        delete(datf);
                        break
                    case 'No'
                        def = {'data_name'};
                        prompt = {'enter a different data file name (must be same as model name)'};
                        titles  = 'Data file name';
                        datf = char(inputdlg(prompt,titles,1,def)); %model file name
                        datf=[datf,'.data'];
                end % switch
            end

            fid=fopen(datf,'w'); % Print header and station locations
            fprintf(fid,'%d %d %d\n',length(H.d.x),length(frtp),H.NR);
            fprintf(fid,'%s\n','Station_Location: N-S');
            for kf=1:length(H.d.x)
                fprintf(fid,'%.2f %s',H.d.x(kf),' ');
            end
            fprintf(fid,'\n');
            fprintf(fid,'%s\n','Station_Location: E-W');
            for kf=1:length(H.d.y)
                fprintf(fid,'%.2f %s',H.d.y(kf),' ');
            end
            fprintf(fid,'\n');

            for I=1:length(frtp) % Print Data
                text1='DATA_Period:       ';
                fprintf(fid,'%s %5d\n',text1,T(I));                
                if resp == 2 % Case added by TB and DR on 23/01/09, edited by GN May 2012
                    fprintf(fid,'%1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d\n',[Z(:,:,I)'; Tip(:,:,I)']);  %Z is conjugated
                elseif resp == 1 % Full Z tensor only                      
                    fprintf(fid,'%1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d %1.4d\n',Z(:,:,I)');  %Z is conjugated, need transpose since fprintf prints column by column
                elseif resp == 3 % off-diag Z only
                    fprintf(fid,'%1.4d %1.4d %1.4d %1.4d\n',Z(:,3:6,I)');  %Z is conjugated
                elseif resp == 4 % Tipper only
                    fprintf(fid,'%1.4d %1.4d %1.4d %1.4d\n',Tip(:,:,I)');
                end
            end
            
            for I=1:length(frtp) % Print Errors
                text1='ERROR_Period:       ';
                fprintf(fid,'%s %5d\n',text1,T(I));
                if resp == 2                       
                    fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d\n ',[EZ(:,:,I)'; ETip(:,:,I)']);
                elseif resp == 1                        
                    fprintf(fid,'%d %d %d %d %d %d %d %d\n',EZ(:,:,I)');
                elseif resp == 3                        
                    fprintf(fid,'%d %d %d %d\n',EZ(:,3:6,I)');
                elseif resp == 4
                    fprintf(fid,'%d %d %d %d\n',ETip(:,:,I)');
                end
            end
            
            for I=1:length(frtp) % Print Errormap 
                text1='ERMAP_Period:       ';
                fprintf(fid,'%s %5d\n',text1,T(I));
                if resp == 2                        
                    fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d\n ',[ERZ(:,:,I)'; ERTip(:,:,I)']);
                elseif resp == 1                       
                    fprintf(fid,'%d %d %d %d %d %d %d %d\n',ERZ(:,:,I)');
                elseif resp == 3                        
                    fprintf(fid,'%d %d %d %d\n',ERZ(:,3:6,I)');
                elseif resp == 4
                    fprintf(fid,'%d %d %d %d\n',ERTip(:,:,I)');
                end
            end
            
            fclose(fid);
            
            % Commented this out 2016: not always accurate, and can
            % calculate yourself from manual... and now you don't have to load a model
            % before exporting data.
%             N=length(H.x)*length(frtp)*H.NR;
%             M=H.nx*H.ny*H.nz;
%             N1=length(H.x)*H.NR; % number of stations times number of responses - this is responses per period
%             MEM_parallel= ( ( (M.*N1.*2) + (N.*N.*2) ).*8 )./(1.024E+09) ; % 2 is required to swap the data; 8 is required because all data is in double precision
%             MEM_serial=( 1.2.*( (8.*N.*N) + (8.*N.*M) ) ) / (1.024E09); %the 1.2 term varies anywhere from 1.1 to 1.4 (or larger)
%             warndlg(['Serial code needs ',num2str(MEM_serial,'%.2f'),' GB, Parallel needs ',num2str(MEM_parallel,'%.2f'),' GB (+- 10%)'])
    end
end

    guidata(hObject, H);
              
end