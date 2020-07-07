function M3_generate_model(hObject, eventdata, ~)
% Display model in M3DET window using parameters entered on screen.

    H=guidata(hObject);
    H.lay=1;
    if ~isfield(H,'dat')
        warndlg('You have to load a data file before you can generate the model.')
    else
        xspacing = str2double(get(H.xspacing,'string')); % Get the xspacing
        moutcX = str2double(get(H.moutcX,'string')); % Number of mesh out of the core in X
        yspacing = str2double(get(H.yspacing,'string')); % Get the yspacing
        moutcY = str2double(get(H.moutcY,'string')); % Number of mesh out of the core in Y
        firstthick = str2num(get(H.firstthick,'string')); % Get the first thickness
        H.nz = str2double(get(H.nl,'string')); % Get the number of layers
        xcoef = str2double(get(H.Xinc,'string')); % Get the coefficient for off model in X direction
        ycoef = str2double(get(H.Yinc,'string')); % Get the coefficient for off model in Y direction
        zcoef = str2double(get(H.Zinc,'string')); % Get the coefficient for model in Z direction
        res = str2double(get(H.res,'string')); 

        lomax=(max(H.d.y)-min(H.d.y))/2000 ; % long. max... horizontal distance (km) between farthest E-W stations
        H.ny=(ceil(lomax/yspacing)+4)*yspacing; % not sure where the 4 comes from, buffer so that padding doesn't start right at last station
        Y1=linspace(0,H.ny,(ceil(lomax/yspacing)+5))-yspacing; Y1(1:2)=[]; Y2=-Y1; % not sure where the 5 comes from
        lamax=(max(H.d.x)-min(H.d.x))/2000 ;
        H.nx=(ceil(lamax/xspacing)+4)*xspacing; % not sure where the 5 or 6 come from - just make them the same 
        X1=linspace(0,H.nx,(ceil(lamax/xspacing)+5))-xspacing; X1(1:2)=[]; X2=-X1;
       
        H.XX=sort([X2,X1,0]); H.YY=sort([Y2,Y1,0]);
        
        if ismember(unique(H.d.x),H.XX.*1000) & ismember(unique(H.d.y),H.YY.*1000) % synthetic data check
            H.XX = H.XX - uniquetol(diff(H.XX),10^-10)/2;
            H.YY = H.YY - uniquetol(diff(H.YY),10^-10)/2; 
            H.XX(1) = [];
            H.YY(1) = [];
        end        
        
        xpn(1)=H.XX(1)+(H.XX(1)-H.XX(2))*xcoef; % add negative padding
        xpn(2)=xpn(1)+(xpn(1)-H.XX(1))*xcoef;
        for i=3:moutcX
            xpn(i)=xpn(i-1)+(xpn(i-1)-xpn(i-2))*xcoef;
        end
        
        xpp(1)=H.XX(end)+(H.XX(end)-H.XX(end-1))*xcoef; % add positive padding
        xpp(2)=xpp(1)+(xpp(1)-H.XX(end))*xcoef;
        for i=3:moutcX
            xpp(i)=xpp(i-1)+(xpp(i-1)-xpp(i-2))*xcoef;
        end
        
        ypn(1)=H.YY(1)+(H.YY(1)-H.YY(2))*ycoef; % add negative padding
        ypn(2)=ypn(1)+(ypn(1)-H.YY(1))*ycoef;
        for i=3:moutcY
            ypn(i)=ypn(i-1)+(ypn(i-1)-ypn(i-2))*ycoef;
        end
        
        ypp(1)=H.YY(end)+(H.YY(end)-H.YY(end-1))*ycoef; % add positive padding
        ypp(2)=ypp(1)+(ypp(1)-H.YY(end))*ycoef;
        for i=3:moutcY
            ypp(i)=ypp(i-1)+(ypp(i-1)-ypp(i-2))*ycoef;
        end
        
        %thicknesses
        if length(firstthick)==1
            z = cumprod([firstthick ones(1,H.nz-1).*zcoef]);
        else
            z = firstthick;
            H.nz = length(z);
        end
        z=cumsum(z);z=sort([z,0]); % convert from thicknesses to depths
        H.Z=round(z)/1000; H.XX=sort([H.XX,xpn,xpp]); H.YY=sort([H.YY,ypn,ypp])'; %x, y and z defined as vectors

        H.nx=length(H.XX); H.ny=length(H.YY); H.nz=length(H.Z);
        H.AA=ones(H.nx,H.ny,H.nz)*res; %Default model
       
        m.x = H.XX';
        m.y = H.YY;
        m.z = H.Z';
        m.A = H.AA;
        m.A(end,:,:) = [];
        m.A(:,end,:) = [];
        m.A(:,:,end) = [];
        m.origin = [min(H.XX) min(H.YY) min(H.Z)].*1000; % put in m
        m.dx = abs(diff(H.XX))'.*1000;
        m.dy = abs(diff(H.YY)).*1000;
        m.dz = abs(diff(H.Z))'.*1000;
        m.cx = (m.x(1:end-1)+m.x(2:end))/2;
        m.cy = (m.y(1:end-1)+m.y(2:end))/2;
        m.cz = (m.z(1:end-1)+m.z(2:end))/2;
        m.nx = H.nx-1; % m structure counts number of cells, not edges
        m.ny = H.ny-1;
        m.nz = H.nz-1; 
        
        %Create a meshgrid of cell centers
        [m.X,m.Y]=meshgrid(m.cy,m.cx);

        %Determine elevation topography surface
        m.Z = zeros(m.nx,m.ny);
        for i = 1:m.nx
            for j = 1:m.ny

                ind = find(isnan(squeeze(m.A(i,j,:))),1,'last');

                if isempty(ind)
                    ind = 0;
                end

                m.Z(i,j) = m.cz(ind+1);

            end
        end
        
        %Determine padding
        m.npad(1) = (m.nx - length(m.dx(round(m.dx)==round(m.dx( round(m.nx/2)))))) /2;
        m.npad(2) = (m.ny - length(m.dy(round(m.dy)==round(m.dy(round(m.ny/2)))))) /2;

        m.name = ['M3D_',num2str(m.nx),'x',num2str(m.ny),'y',num2str(m.nz),'z'];
        m.niter = ''; %Number of iterations left blank and can be added later
        
        H.top = m.origin(3); % not sure if needed
    
        H.m = m;

        M3_update_model_plot(H)

        if ~isfield(H,'undo')        
            H.undo=0;    
        end
        H.undo=H.undo+1;

%         M3_save_undo(hObject, eventdata, H)

        guidata(hObject, H) 
    end

end % end generate_model_Callback