function M3_load_winglink_mod(hObject, ~, ~)
% Load WinGLink format model file (not heavily tested)
% Dec 2019 - lines added from load_model_modem to make m structure
    
    H=guidata(hObject);
    delete undo*
    
    if ~isfield(H,'dat')
        warndlg('Load data file first')
        return
    end
    
    H.lay=1;

    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1)

    %Load model
    [meshf,~] = uigetfile({'*.out'},'pick winglink *.out file');
    fid=fopen(meshf,'r');
    [NN]=fscanf(fid,'%d %d %d %d'); % first 3 values: no. of EW cells, no. of NS cells, no. of layers
    H.ny = NN(1); H.nx = NN(2); H.nz = NN(3);
    bx=fscanf(fid,'%s',1); % junk string
    
    H.YY=fscanf(fid,'%f',H.ny);
    H.XX=fscanf(fid,'%f',H.nx);
    H.Z=round(fscanf(fid,'%f',H.nz));
           
    slj=[];
    for uu=1:H.nz
       nul=fscanf(fid,'%f',1); %We don't need this
       sli=fscanf(fid,'%f',H.nx*H.ny);
       slj=abs([slj; sli]);
    end
    for ofbe=1:5
       ofb=fgetl(fid);
    end
    coords_struct=textscan(ofb,'%f',2);
    coords=coords_struct{1}; % real world coordinates
    wing_rot = -cell2mat(textscan(fgetl(fid),'%f',1));
    top = cell2mat(textscan(fgetl(fid),'%f',1));
%     wing_rot=-fscanf(fid,'%f',1); %if this is nonzero station coordinate frame will be rotated.
    H.mesh_rot=wing_rot;
%     H.rot_cor_ang=H.mesh_rot;
    res_tmp=unique(slj);

    if length(res_tmp)>9
       warndlg('You cannot have more than 9 resistivities, program will simplify your model');

       resis9=round(logspace(0,3,9));
       prompt={'Enter the resistivity values to index the model - increasing resistivity'};
       name='Maximum 9 values-increasing resistivity order';
       numlines=1;
       defaultanswer={num2str(resis9)};
       answer=inputdlg(prompt,name,numlines,defaultanswer);
       resis9=str2num(char(answer));
       resdif=diff(resis9); resdif=[resdif max(resdif)*100];
       %slj(find(slj<min(resis9)))=resis9(2);

       for k=1:8
           slj(find(slj>resis9(k)+resdif(k)/2 & slj<=resis9(k+1)+resdif(k+1)/2))=resis9(k+1);
       end
       slj(find(slj>resis9(1) & slj<resis9(2)))=resis9(2);
       if length(find(slj<resis9(1)))>0
           warndlg('Model has resistivities lower than your lowest resistivity value. To avoid complications these values are replaced by your second resistivity value.');
           slj(find(slj<resis9(1)))=resis9(2);
       end

    end
    res_tmp=unique(slj);

    H.AA=reshape(slj,H.ny,H.nx,H.nz);
    H.AA = permute(H.AA,[2 1 3]); % get into form of m structure

    %Program will only work for nine distinct resistivity values in model file
    nnr= length(res_tmp); % this in only true for wsinv
    fclose(fid);

    %now reset the model variables so they match the generated model variables
    % H.nx, H.ny, H.nz, H.XX, H.YY, H.Z, H.AA include ALL model edges - so they will have
    % 1 more entry in each direction than the variables in the m structure
       
    m.dx = H.XX;
    m.dy = H.YY;
    m.dz = H.Z;
    m.A = H.AA;
    m.origin = [coords(1) coords(2) top];

    m.nx = H.nx;
    m.ny = H.ny;
    m.nz = H.nz;     
    
    H.XX=((H.XX)./1000)';
    x_array=cumsum(H.XX);
    mid_x_ind=find( abs((max(x_array)/2)-x_array) == min( abs(  (max(x_array)/2)-x_array)  )   ); %find the indice of the middle of the x array
    H.XX=x_array-x_array(mid_x_ind(1)); %put 0 at the center of the mesh
    XX_tmp(1)=-1*H.XX(H.nx);
    XX_tmp(2:H.nx+1)=H.XX;
    H.XX=XX_tmp';
    
    H.YY=(H.YY)./1000;
    y_array=cumsum(H.YY);
    mid_y_ind=find( abs((max(y_array)/2)-y_array) == min( abs(  (max(y_array)/2)-y_array)  )   ); %find the indice of the middle of the y array
    H.YY=y_array-y_array(mid_y_ind(1));%put 0 at the center of the mesh
    YY_tmp(1)=-1*H.YY(H.ny);
    YY_tmp(2:H.ny+1)=H.YY;
    H.YY=YY_tmp';
    
    H.Z=(H.Z./1000)';
    H.Z=cumsum(H.Z);
    Z_tmp(1)=0;
    Z_tmp(2:H.nz+1)=H.Z;
    H.Z=Z_tmp';
    
    m.x = H.XX;
    m.y = H.YY;
    m.z = H.Z;
    
    %Find cell center midpoints
    m.cx = (m.x(1:end-1)+m.x(2:end))/2;
    m.cy = (m.y(1:end-1)+m.y(2:end))/2;
    m.cz = (m.z(1:end-1)+m.z(2:end))/2;

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
    m.npad(1) = (m.nx - length(m.dx(m.dx==m.dx(round(m.nx/2)))))/2;
    m.npad(2) = (m.ny - length(m.dy(m.dy==m.dy(round(m.ny/2)))))/2;
    
    m.name = meshf;
    m.niter = ''; %Number of iterations left blank and can be added later

    H.AA(H.nx+1,:,:)=H.AA(H.nx,:,:);
    H.AA(:,H.ny+1,:)=H.AA(:,H.ny,:);
    H.AA(:,:,H.nz+1)=H.AA(:,:,H.nz);

    H.nx=H.nx+1;
    H.ny=H.ny+1;
    H.nz=H.nz+1;
   
    H.top = m.origin(3); % not sure if needed
    
    H.m = m;
  
    M3_update_model_plot(H)

    guidata(hObject, H); %save H back to gui

end % end load_winglink_mod_Callback