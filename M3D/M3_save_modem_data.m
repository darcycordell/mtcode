function M3_save_modem_data(hObject, ~, ~)
    
    H=guidata(hObject);
    
    if ~isfield(H,'dat') % check if data has been loaded
    
        disp('No data to save! Load data first.');
        warndlg('Enter data parameters first.')
    
    else % data has been loaded
    
        if ~isfield(H,'max_per')||~isfield(H,'min_per')||~isfield(H,'per_skip')
            warndlg('Error: enter values for min per, max per, or periods to skip.')
            return
        end
        
        [modem_datname,path]=uiputfile('*','Save ModEM Data as:',[pwd,'\','modem.data']);
        modem_dat=[path,modem_datname];
        
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
        Zxx_erflr = .01*str2double(get(H.errflr_Zdiag,'string')); % diagonal Z
        Zxy_erflr = .01*str2double(get(H.errflr_Zodiag,'string')); % off-diagonal Z
        Zyx_erflr = .01*str2double(get(H.errflr_Zodiag,'string')); % off-diagonal Z
        Zyy_erflr = .01*str2double(get(H.errflr_Zdiag,'string')); % diagonal Z
        tip_erflr = str2double(get(H.errflr_tip,'string')); % tipper
        err_flr_type = str2double(get(H.errflr_type,'string')); %error denominator type (value of 1 or 2. Default = 1)
              
        % this origin is the center of the stations... z=0 is assumed, but not actually written in data file
        % this origin value is not used directly by ModEM, but makes it clear where the data (and model) is located
        origin=[(max(coord(:,2))-min(coord(:,2)))/2+min(coord(:,2)),(max(coord(:,1))-min(coord(:,1)))/2+min(coord(:,1)),0];
        
        x = H.d.x; % these are in m
        y = H.d.y;
        
        if H.mesh_rot ~= 0
            
            warndlg(['MT data rotated ',num2str(H.mesh_rot),' degrees.'])  
            % first rotate the data
%             [Z,Zvar,Tip,Tvar] = fromersan(ZZZ,EEE,ZZZ_HZ,EEE_HZ); %just changing data format....
%             [Z,Zvar]=rot_z_dr(Z,Zvar,length(T),rot);
%             
%             for is=1:length(Tip(1,1,:))
%                 [Tip(:,:,is),Tvar(:,:,is)] = rot_tip(Tip(:,:,is),Tvar(:,:,is),rot);
%             end
%             
%             [ZZZ,EEE,ZZZ_HZ,EEE_HZ] = toersan(Z,Zvar,Tip,Tvar); % and changing data format back again!  
% 
%             % now rotate the station locations:
%             R=[cosd(rot),sind(rot);-sind(rot),cosd(rot)]; % clockwise rotation matrix
%             
%             for is=1:length(Tip(1,1,:))
%                 tmp=R*[x(is);y(is)];
%                 x(is)=tmp(1);
%                 y(is)=tmp(2);
%             end           
        end       
        
        % Old versions of M3 (pre-2018) had some unit conversions. The new
        % MTcode is standardized to only include SI units and thus, unit
        % conversions are no longer necessary and Z_scale is set to 1.
        Z_scale = 1;
        data(1:4,:,:)=ZZZ(1:4,:,:)./Z_scale;
        dataerr(1:4,:,:)=EEE(1:4,:,:)./Z_scale;
        data(5:6,:,:)=ZZZ_HZ(1:2,:,:);
        dataerr(5:6,:,:)=EEE_HZ(1:2,:,:);
        
        dataerr(isnan(dataerr)) = 1e-9 ; % MJU May 24, 2014
        
        if err_flr_type == 1 % sqrt(abs(Zxy*Zyx))
            
            for iper=1:length(T)
                for istn=1:ns               
                    val=sqrt(abs(data(2,iper,istn)*data(3,iper,istn)));
         
                    if dataerr(1,iper,istn) < val*Zxx_erflr
                        dataerr(1,iper,istn)=val*Zxx_erflr;
                    end
                    if dataerr(2,iper,istn) < val*Zxy_erflr
                        dataerr(2,iper,istn)=val*Zxy_erflr;
                    end
                    if dataerr(3,iper,istn) < val*Zyx_erflr
                        dataerr(3,iper,istn)=val*Zyx_erflr;
                    end
                    if dataerr(4,iper,istn) < val*Zyy_erflr
                        dataerr(4,iper,istn)=val*Zyy_erflr;
                    end
                    if dataerr(5,iper,istn) < tip_erflr
                        dataerr(5,iper,istn)=tip_erflr;
                    end
                    if dataerr(6,iper,istn) < tip_erflr
                        dataerr(6,iper,istn)=tip_erflr;
                    end               
                end
            end
            
        elseif err_flr_type == 2 % Zxx error = % of |Zxy|, and so on
            
            for iper=1:length(T)
                for istn=1:ns                
                             
                    if dataerr(1,iper,istn) < abs(data(2,iper,istn))*Zxx_erflr
                        dataerr(1,iper,istn)=abs(data(2,iper,istn))*Zxx_erflr;
                    end
                    if dataerr(2,iper,istn) < abs(data(2,iper,istn))*Zxy_erflr
                        dataerr(2,iper,istn)=abs(data(2,iper,istn))*Zxy_erflr;
                    end
                    if dataerr(3,iper,istn) < abs(data(3,iper,istn))*Zyx_erflr
                        dataerr(3,iper,istn)=abs(data(3,iper,istn))*Zyx_erflr;
                    end
                    if dataerr(4,iper,istn) < abs(data(3,iper,istn))*Zyy_erflr
                        dataerr(4,iper,istn)=abs(data(3,iper,istn))*Zyy_erflr;
                    end
                    if dataerr(5,iper,istn) < tip_erflr
                        dataerr(5,iper,istn)=tip_erflr;
                    end
                    if dataerr(6,iper,istn) < tip_erflr
                        dataerr(6,iper,istn)=tip_erflr;
                    end               
                end
            end
            
        else
            
            warndlg('Choose Error Type 1 or 2')
            
        end
        
        if T(1)<T(end)
            minid=find(min_per<=T,1,'first'); %DC Corrected so that it gets proper indices
            maxid=find(max_per>=T,1,'last'); %DC Corrected so that it gets proper indices
            frtp = minid:per_skip:maxid; % frequencies to plot  
        elseif T(1)>T(end)
            minid=find(min_per>=T,1,'last'); %DC Corrected so that it gets proper indices
            maxid=find(max_per>=T,1,'first'); %DC Corrected so that it gets proper indices
            frtp = minid:-per_skip:maxid; % frequencies to plot  
        end
                  
        set(H.data_text,'string',sprintf('%s\n%s\n%s',['Number of periods = ',num2str(numel(frtp))],['Minimum: ',num2str(T(min(frtp))),' s'],['Maximum: ',num2str(T(max(frtp))),' s']))

        disp('Periods being used:');
        T(frtp)
        
        % NOTE that winglink exports edi files as e^(+iwt), and wsinv3d uses e^(-iwt) (data in the matfile is in e^(+iwt))
        if get(H.conj_Z,'value')            
            data=conj(data);
        end
                    
        siteLoc=[x,y,coord(:,3)]; % DC Edit --> included z coordinates (coord(:,3)) instead of zeros(length(y),1))
        coord_tmp = coord;
        siteLoc_tmp = siteLoc; % column 1 is N-S, column 2 is E-W
            
        %DC Edit: Added lines to move the stations to center of grid cells and
        %project sites onto topography (so that stations aren't underground or
        %in the air)
%         midHx = round(H.nx/2); %mid point in x and y
%         midHy = round(H.ny/2);
%         dx=H.XX(midHx+1)-H.XX(midHx); %station spacing
%         dy = H.YY(midHy+1)-H.YY(midHy);       
        if isfield(H,'m') % check to see if there is a model
            for i = 1:length(coord_tmp)

                %DC: Comment for synthetics when stations are exactly on grid vertices

    %             siteLoc_tmp(i,1)=(H.YY(indx)+dy/2)*1000; %Move the station to the middle of that cell and over-write station location info
    %             siteLoc_tmp(i,2) = (H.XX(indy)+dx/2)*1000;           

                indx1 = nearestpoint(siteLoc(i,1),H.XX*1000,'previous'); %(dependency "nearestpoint")
                indx2 = indx1 + 1;
                indy1 = nearestpoint(siteLoc(i,2),H.YY*1000,'previous'); %Find the cell index (x,y) that the station is located in
                indy2 = indy1 + 1;

                siteLoc_tmp(i,1) = ((H.XX(indx1)+H.XX(indx2))/2)*1000; 
                siteLoc_tmp(i,2) = ((H.YY(indy1)+H.YY(indy2))/2)*1000; %Move the station to the middle of that cell and over-write station location info

                if ~isfield(H,'AAt')
                     [indz]= min(find(H.AA(indx1,indy1,:)<10^16)); %Move station elevation coordinate to the surface of the topo
                     if isfield(H,'top') % currently, this will only be the case if you generate a model and add topo, or if you load a ModEM data file 
                        siteLoc_tmp(i,3) =((H.Z(indz)*1000)); %Note this is m.b.s.l. so station elevations are negative!
                     else
                        siteLoc_tmp(i,3) = H.Z(indz); %If there is no topography added then set elevations to zero.
                     end
                elseif isfield(H,'AAt')
                    [indz]= min(find(H.AAt(indx1,indy1,:)<10^16)); %Move station elevation coordinate to the surface of the topo
                    siteLoc_tmp(i,3) =((H.Z(indz)*1000)); %Note this is m.b.s.l. so station elevations are negative!
                end
            end
        else
            warndlg('Data file is being written without a corresponding mesh!')
            siteLoc_tmp(:,3) = 0; % just set the elevations to zero
        end
        %End DC Edit
        
        data = permute(data,[2 1 3]);
        H.dex.Z = data(frtp,1:4,:);
        H.dex.tip = data(frtp,5:6,:);
        
        dataerr = permute(dataerr,[2 1 3]);
        H.dex.Zerr = dataerr(frtp,1:4,:);
        H.dex.tiperr = dataerr(frtp,5:6,:);
        
        if get(H.resp,'value') == 2 % Full tensor and tipper
            nr = 12;
        elseif get(H.resp,'value') == 1 % Full impedance tensor
            nr=8;
        elseif get(H.resp,'value') == 3 % off diagonals only
            nr=4;
        elseif get(H.resp,'value') == 4 % tip only
            nr=2; % not sure why this is 2?
        end 
        
        nr=nr/2; %instead of 4,8,or 12 responses, we have 1, 2, 4, or 6 responses 
        if nr == 2
            H.dex.responses={'','ZXY' 'ZYX',''};
            H.dex.tip = NaN+1i*NaN;
            H.dex.tiperr = NaN+1i*NaN;
            H.dex.Z(:,[1 4],:) = NaN+1i*NaN;
            H.dex.Zerr(:,[1 4],:) = NaN+1i*NaN;
            H.dex.rho(:,[1 4],:) = NaN;
            H.dex.rhoerr(:,[1 4],:) = NaN;
            H.dex.pha(:,[1 4],:) = NaN;
            H.dex.phaerr(:,[1 4],:) = NaN;
        elseif nr == 1 %DC Edit: Bug Fix to write out tipper only
            H.dex.responses = {'','','','','','TX ' 'TY '};
            H.dex.Z = NaN+1i*NaN;
            H.dex.Zerr = NaN+1i*NaN;
            H.dex.rho = NaN;
            H.dex.rhoerr = NaN;
            H.dex.pha = NaN;
            H.dex.phaerr = NaN;
        elseif nr == 6
            H.dex.responses = {'ZXX' 'ZXY' 'ZYX' 'ZYY' 'TX ' 'TY '};
        else
            H.dex.responses = {'ZXX' 'ZXY' 'ZYX' 'ZYY'};
            H.dex.tip(:,:,:) = NaN+1i*NaN;
            H.dex.tiperr(:,:,:) = NaN+1i*NaN;
        end
                   
        H.dex.origin = [H.cent_lat H.cent_long];
        
        H.dex.zrot = H.mesh_rot;
        H.dex.trot = H.mesh_rot;
        
        H.dex.T = T(frtp);
        H.dex.f = 1./H.dex.T;
        H.dex.nf = numel(T(frtp));
        H.dex.ns = ns;
        
        H.dex.x = siteLoc_tmp(:,1);
        H.dex.y = siteLoc_tmp(:,2);
        H.dex.z = siteLoc_tmp(:,3);
                     
        % only coords in meters are centered? original lat/lon are being exported          
        write_data_modem(modem_dat,H.dex);
                        
    end
    
    guidata(hObject, H);
    
end % end save_modem_data_Callback