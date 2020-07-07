function M3_make_modem_cov(hObject, ~, ~)
% Write a ModEM .cov file for inversion, editable values
    
    H=guidata(hObject);

    if isfield(H,'m') % must have model resistivities before making cov file

        prompt = {sprintf('Enter .cov file settings\n\nFile name'),'x smoothing (Default: 0.3)','y smoothing (Default: 0.3)','z smoothing (Default: 0.3)',sprintf('number of times smoothing applied\n(Default: 1)')};
        dlg_title = 'Cov File';
        def = {'modem.cov','0.3','0.3','0.3','1'};
        answer = inputdlg(prompt,dlg_title,1,def);
    
        indm = ones(H.nx-1,H.ny-1,H.nz-1);
        
        if isfield(H,'AAt') % if topo was added in current M3 session    
            A = H.AAt(1:H.nx-1,1:H.ny-1,1:H.nz-1);
        else  
            A = H.AA(1:H.nx-1,1:H.ny-1,1:H.nz-1);
        end
        
        if get(H.fix_ocean,'Value') %If "fix ocean" box is checked
            indm(A==str2double(get(H.ocean_res,'string'))) = 9; % fix ocean cell values with 9.
        end

        indm(A==10^17) = 0; % cov file sets air cells to zero smoothing 
        indm(isnan(A)) = 0; % if imported model, air cells are NaN
                                   
        nzz=H.nz-1;
        cfile=answer{1};
        cov = indm; % rot90(fliplr(indm)); % previously the xy swap was done here
        xsmooth = str2double(answer{2})*ones(nzz,1); % these need to be vectors
        ysmooth = str2double(answer{3})*ones(nzz,1);
        zsmooth = str2double(answer{4});
        nsmooth = str2double(answer{5});
        % note cov file format is different than model file... values for
        % each layer are printed as nx x ny
        [status] = write_modem_covariance(cfile,cov,xsmooth,ysmooth,zsmooth,nsmooth);
    else
        warndlg('No matrix of resistivity indices found- save ModEM model file first!','Error')
    end

    guidata(hObject, H);

end % end make_modem_cov_Callback