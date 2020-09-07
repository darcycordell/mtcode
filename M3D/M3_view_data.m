function M3_view_data(hObject, ~, ~)

H=guidata(hObject);

if ~isfield(H,'dat')  % check if data has been loaded
    disp('No data to view! Load data first.');
    warndlg('Enter data parameters first.')
    
else % data has been loaded
    
    % check for input for min per max per and per skip
    if ~isfield(H,'max_per')||~isfield(H,'min_per')||~isfield(H,'per_skip')
        warndlg('Error: enter values for min per, max per, or periods to skip.')
        return
    end
    
    % give variables shorter names
    ZZZ = permute(H.d.Z,[2 1 3]); % to handle new d struct
    EEE = permute(H.d.Zerr,[2 1 3]);
    ZZZ_HZ = permute(H.d.tip,[2 1 3]);
    EEE_HZ = permute(H.d.tiperr,[2 1 3]);
    site = H.d.site;
    T = H.d.T;
    max_per = str2double(get(H.max_per,'string'));
    min_per = str2double(get(H.min_per,'string'));
    per_skip = str2double(get(H.per_skip,'string'));
    err_flr_Zd = .01*str2double(get(H.errflr_Zdiag,'string')); % diagonal Z
    err_flr_Zod = .01*str2double(get(H.errflr_Zodiag,'string')); % off-diagonal Z
    err_flr_tip = str2double(get(H.errflr_tip,'string')); % tipper
    err_flr_type = str2double(get(H.errflr_type,'string')); %error denominator type (value of 1 or 2. Default = 1)
%     H.rot_cor_ang = H.mesh_rot;
         
    if get(H.conj_Z,'value')
        ZZZ = conj(ZZZ);
    end
    if get(H.conj_tip,'value') % GN - added ability to take conjugate of tipper (April 2012)
        ZZZ_HZ=conj(ZZZ_HZ);
    end
           
    if H.mesh_rot ~= 0 % Rotate impedance data to opposite direction of the mesh if needed
        warndlg(['MT impedance data rotated ',num2str(H.mesh_rot),' degrees.'])
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
    
    % Older versions of M3 (e.g. pre-2018) performed two conversions: 1) changed from field
    % units to SI and 2) changed errors from variances to standard errors (std dev)
    %
    % The new MTcode is always in SI units and stores standard errors rather
    % than variances. So it is no longer necessary to perform any conversions.
    %
    % The old way: EEE = sqrt(EEE)/796;
    %
    
    % Impose error floor of the type specified by user:
    % Edits by DC 18/11/14
        if err_flr_type==1 % Denominator is sqrt(abs(Zxy*Zyx)) for ALL Z components
        
            [EEE,EEE_HZ,ZZZ,calc_err,~,percent_err_replaced]=M3_err_flr_type_1(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod,err_flr_tip);
        
        elseif err_flr_type==2 % Denominator is abs(Zxy) for Zxx,Zxy. Denominator is abs(Zyx) for Zyx,Zyy
        
            [EEE,EEE_HZ,ZZZ,calc_err,~,percent_err_replaced]=M3_err_flr_type_2(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod,err_flr_tip);
            
        elseif err_flr_type==3 % The OLD WAY: Original code untouched
            nfim=length(T);
            % DC EDIT Jan 21, 2015 for plotting purposes
            %calc_err_edi=[];
%             for ic=1:4
%                 if ic==1 || ic==4
%                     calc_err_edi(ic,:,:)=(EEE(ic,:,:))./sqrt(abs((ZZZ(2,:,:).*(ZZZ(3,:,:)))));
%                 else
%                     calc_err_edi(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(ic,:,:));
%                 end
%             end
            %end DC edits Jan 21, 2015
            N_nan = length(find(isnan(ZZZ(:))==1)); %Number of nan data points
            N = length(ZZZ(:))-N_nan; %Total number of non-nan data points
            count = 0; tip_count = 0;
            for ic=1:4
                for itrpo=1:length(H.d.x)
                    for iu=1:nfim   %number of frequencies in *mat file
                        if ic==1 || ic==4
                            if EEE(ic,iu,itrpo)<err_flr_Zd*sqrt((ZZZ(2,iu,itrpo).*conj((ZZZ(3,iu,itrpo))))) %Conjugate is irrelevant. I think he should be comparing modulus --- DC OCT 28/14
                                count = count + 1;
                                EEE(ic,iu,itrpo)=err_flr_Zd*sqrt(abs(ZZZ(ic,iu,itrpo).*ZZZ(ic,iu,itrpo)));   %10 errorfloor
                            end
                        elseif ic==2 || ic==3
                            count = count + 1;
                            EEE(ic,iu,itrpo)=err_flr_Zod*sqrt(abs(ZZZ(ic,iu,itrpo).*ZZZ(ic,iu,itrpo)));   %10 errorfloor
                        end
                        if (ic ==1 || ic == 2) && (H.resp == 2 || H.resp == 4) %tipper
                            if EEE_HZ(ic,iu,itrpo) < err_flr_tip
                                tip_count = tip_count + 1;
                                EEE_HZ(ic,iu,itrpo)=err_flr_tip;   %absolute error floor for tipper error
                            end
                        end
                    end
                end
            end 

            percent_err_replaced = count/N;
            %percent_nan = N_nan/N;
            
            %DC EDITS Jan 21, 2015 for plotting purposes
            calc_err=nan(size(ZZZ));
            for ic=1:4
                if ic==1 || ic==4
                    calc_err(ic,:,:)=(EEE(ic,:,:))./sqrt(abs((ZZZ(2,:,:).*(ZZZ(3,:,:)))));
                else
                    calc_err(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(ic,:,:));
                end
            end

       %End DC EDITS Jan 21, 2015 
        
        else %Displays error: err_flr_type must be 1 or 2. Defaults to err_flr_type=1
            fprintf('Error Floor Type must be 1, 2 or 3. \nYour entry was invalid: Default error floor is Type 1\nSee M3DET documentation for more details')
            [EEE,EEE_HZ,ZZZ,calc_err,~,percent_err_replaced]=M3_err_flr_type_1(EEE,ZZZ,EEE_HZ,T,H,err_flr_Zd,err_flr_Zod);   
        end
        
        %End Edits by DC 18/11/14
           
    minid=max(flipud(find(min_per>T))) + 1; % adding 1 ensures that all periods are greater than per_min
    maxid=min(find(max_per<T)) - 1; % subtracting 1 ensures that all periods are less than per_max
    if isempty(minid)
        minid = 1;
    end
    if isempty(maxid)
        maxid = numel(T);
    end
    frtp = minid:per_skip:maxid; % frequencies to plot
    
    % update the heading with number of periods
    set(H.data_text,'string',sprintf('%s\n%s\n%s',['Number of periods = ',num2str(numel(frtp))],['Minimum: ',num2str(T(min(frtp))),' s'],['Maximum: ',num2str(T(max(frtp))),' s']))
    
    per_range = sum(isnan(ZZZ),3);
    perplot = unique(per_range,'rows');
    if size(perplot,1) ~= 1
        disp('Warning: 1 or sites has an empty impedance tensor element (e.g. Zxx is missing but the other tensor elements are present). This type of dataset has not been thoroughly tested.')
        perplot = max(perplot,[],1);
        %return
    end
    figure
    bar(log10(T),numel(site)-perplot,0.6);
    hold on
    plot(log10(T(frtp)),numel(site)-perplot(frtp),'rs','markerfacecolor','r')
    plot(log10(T),ones(numel(T)).*numel(site),'k--')
    legend('No. of Stations','Period selected','Total No. of Stations')
    title(['Number of stations with data. ',num2str(numel(T)),' Periods, ',num2str(numel(frtp)),' selected for inversion'])
    xlabel('Log10 Period (s)')
    ylabel('Number of Stations with Data')
%     bar(log10(T(frtp)),numel(site)-perplot(frtp),0.7,'r')
    
    %Print data
    
    is=1;
    chk=1;
    while chk==1
        %%
        set_figure_size(101);
        clf
        
        mu=4*pi*10^-7; Z_plot = nan(size(ZZZ)); Zerr = nan(size(ZZZ));
        for ir=1:4
            Z_plot(ir,:,is)=real(sqrt(T(:,1)./(2*pi*mu))'.*(ZZZ(ir,:,is)))+1i*imag(sqrt(T(:,1)./(2*pi*mu))'.*(ZZZ(ir,:,is)));  %MVL June 28 add 'quote'
            Zerr(ir,:,is)=real(sqrt(T(:,1)./(2*pi*mu))'.*(EEE(ir,:,is)))+1i*imag(sqrt(T(:,1)./(2*pi*mu))'.*(EEE(ir,:,is)));
        end
        
        subplot(5,2,[1,3])
        logerrorbar(T,abs(real(Z_plot(1,:,is)))',abs(real(Zerr(1,:,is)))','go', '-g'); hold on;
        loglog(T(frtp),squeeze(abs(real(Z_plot(1,frtp,is)))),'k*', 'MarkerSize',3); hold on;
        logerrorbar(T,abs(real(Z_plot(4,:,is)))',abs(real(Zerr(4,:,is)))','cs', '-c'); hold on;
        loglog(T(frtp),squeeze(abs(real(Z_plot(4,frtp,is)))),'k*', 'MarkerSize',3); hold on; grid on
        
%         if is==1
            ax=axis;
            if ax(3)<min(min(min(abs(real(Z_plot(:,:,is))), abs(imag(Z_plot(:,:,is))))))
                ax(3) = 0.5*min(min(min(abs(real(Z_plot(:,:,is))), abs(imag(Z_plot(:,:,is))))));
            end
            
            if ax(4)>max(max(max(abs(real(Z_plot(:,:,is)))+abs(real(Zerr(:,:,is))),abs(imag(Z_plot(:,:,is)))+abs(imag(Zerr(:,:,is))))))
                ax(4) = 2*max(max(max(abs(real(Z_plot(:,:,is)))+abs(real(Zerr(:,:,is))),abs(imag(Z_plot(:,:,is)))+abs(imag(Zerr(:,:,is))))));
            end
%             ax(3)=ax(3)*10^2;
%             ax(4)=ax(4)*10^2;
%         end
        axis([T(frtp(1))-T(frtp(1))/4 T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 ax(3) ax(4)]);
        %set(gca,'YTick',[0.1 1 10 100 1000 10000]); set(gca,'XTick',[0.01 0.1 1 10 100]);
        %set(gca,'YTickLabel',{'-1';'0';'1';'2';'3';'4'}); set(gca,'XTickLabel',{'-2';'-1';'0';'1';'2'});
        title(['Station: ',site{is},'Real impedance XX & YY']);
        
        subplot(5,2,[2 4])
        logerrorbar(T,abs(real(Z_plot(2,:,is)))',abs(real(Zerr(2,:,is)))','ro', '-r'); hold on;
        loglog(T(frtp),squeeze(abs(real(Z_plot(2,frtp,is)))),'k*', 'MarkerSize',3); hold on;
        logerrorbar(T,abs(real(Z_plot(3,:,is)))',abs(real(Zerr(3,:,is)))','bs', '-b'); hold on;
        loglog(T(frtp),squeeze(abs(real(Z_plot(3,frtp,is)))),'k*', 'MarkerSize',3); hold on; grid on
               
        axis([T(frtp(1))-T(frtp(1))/4 T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 ax(3) ax(4)]);
        %set(gca,'YTick',[0.1 1 10 100 1000 10000]); set(gca,'XTick',[0.01 0.1 1 10 100]);
        %set(gca,'YTickLabel',{'-1';'0';'1';'2';'3';'4'}); set(gca,'XTickLabel',{'-2';'-1';'0';'1';'2'});
        title('Real impedance XY & YX');
        
        subplot(5,2,[5 7])
        logerrorbar(T,abs(imag(Z_plot(1,:,is)))',abs(imag(Zerr(1,:,is)))','go', '-g'); hold on;
        loglog(T(frtp),squeeze(abs(imag(Z_plot(1,frtp,is)))),'k*', 'MarkerSize',3); hold on;
        logerrorbar(T,abs(imag(Z_plot(4,:,is)))',abs(imag(Zerr(4,:,is)))','cs', '-c'); hold on;
        loglog(T(frtp),squeeze(abs(imag(Z_plot(4,frtp,is)))),'k*', 'MarkerSize',3); hold on; grid on
        
        axis([T(frtp(1))-T(frtp(1))/4 T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 ax(3) ax(4)]);
        %set(gca,'YTick',[-180 -135 -90 -45 0 45 90 135 180]); set(gca,'XTick',[0.01 0.1 1 10 100]);
        %set(gca,'YTickLabel',{'-180';'-135';'-90';'-45';'0';'45';'90';'135';'180'}); set(gca,'XTickLabel',{'-2';'-1';'0';'1';'2'});
        title('Imaginary impedance XX & YY');
        
        subplot(5,2,[6 8])
        logerrorbar(T,abs(imag(Z_plot(2,:,is)))',abs(imag(Zerr(2,:,is)))','ro', '-r'); hold on;
        loglog(T(frtp),squeeze(abs(imag(Z_plot(2,frtp,is)))),'k*', 'MarkerSize',3); hold on;
        logerrorbar(T,abs(imag(Z_plot(3,:,is)))',abs(imag(Zerr(3,:,is)))','bs', '-b'); hold on;
        loglog(T(frtp),squeeze(abs(imag(Z_plot(3,frtp,is)))),'k*', 'MarkerSize',3); hold on; grid on
        
        axis([T(frtp(1))-T(frtp(1))/4 T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 ax(3) ax(4)]);
        %set(gca,'YTick',[-180 -135 -90 -45 0 45 90 135 180]); set(gca,'XTick',[0.01 0.1 1 10 100]);
        %set(gca,'YTickLabel',{'-180';'-135';'-90';'-45';'0';'45';'90';'135';'180'}); set(gca,'XTickLabel',{'-2';'-1';'0';'1';'2'});
        title('Imaginary impedance XY & YX');
        
        subplot(5,2,9)
        loglog(T,real(calc_err(2,:,is)),'ro','MarkerSize',3); hold on
        loglog(T,real(calc_err(3,:,is)),'bx','MarkerSize',3); 
        loglog(T,err_flr_Zd*ones(length(T),1),'--k'); grid on
%         axis([T(frtp(1))-T(frtp(1))/4
%         T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 10^-3
%         2*max(real(calc_err(2,:,is)))]); % commented 12/2019 - error if
%         inversion resp. data with NaN error
        xlabel('Period (s)')
        ylabel('Relative Error')
        title(['Impedance Relative Error with ',num2str(percent_err_replaced*100,3),'% of data errors re-written']);
        manual_legend('Error Floor','--k','XY %Error','ro','YX %Error','bx');
        
        
        subplot(5,2,10)
        errorbar(T,real(ZZZ_HZ(1,:,is)),real(EEE_HZ(1,:,is)),'ro', 'MarkerSize',3); hold on;
        set(gca,'xscale','log')
        semilogx(T(frtp),squeeze(real(ZZZ_HZ(1,frtp,is))),'k*','MarkerSize',3)
        errorbar(T,real(ZZZ_HZ(2,:,is)),real(EEE_HZ(2,:,is)),'bs', 'MarkerSize',3); hold on;
        set(gca,'xscale','log')
        semilogx(T(frtp),squeeze(real(ZZZ_HZ(2,frtp,is))),'k*','MarkerSize',3); grid on
        
        axis([T(frtp(1))-T(frtp(1))/4 T(frtp(length(frtp)))+T(frtp(length(frtp)))/4 -1 1]);
        set(gca,'Yscale','linear');
        %set(gca,'YTick',[-180 -135 -90 -45 0 45 90 135 180]); set(gca,'XTick',[0.01 0.1 1 10 100]);
        %set(gca,'YTickLabel',{'-180';'-135';'-90';'-45';'0';'45';'90';'135';'180'}); set(gca,'XTickLabel',{'-2';'-1';'0';'1';'2'});
        title('Observed Tipper Tx (red) & Ty (blue)');
        
        menu_view_data=menu('','Next Station','Previous Station','Finish');
        if menu_view_data==1
            if is>=length(site)
                is=length(site);
                disp('Reached end of stations')
            else
                is=is+1;
            end
        elseif menu_view_data==2
            if is<=1
                is=1;
                disp('Viewing first station')
            else
                is=is-1;
            end
        else
            chk=0;
            close(figure(101))
        end
        
    end      
end

end % end view_data_Callback