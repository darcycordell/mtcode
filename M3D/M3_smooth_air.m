function M3_smooth_air(hObject, ~, ~)
    H=guidata(hObject);
    
    if isfield(H,'AAt') % topography was added in current M3 session
        isair = squeeze(sum(sum(H.AAt>1e16,1),2)); % find number of NaN cells per layer
        model = H.AAt;
    elseif sum(sum(isnan(H.AA(:,:,1)))) > 0 % model might be loaded with topography
        isair = squeeze(sum(sum(isnan(H.AA),1),2)); % find number of NaN air cells per layer
        model = H.AA;
    else
        warndlg('No mesh with topography found','Error')
        return
    end
    
    figure
    plot(isair,'ko'); hold on; plot(isair,'k-')
    xlabel('Layer number')
    ylabel('Number of air cells')
    
    smth = questdlg('Smooth the last air layer? This will remove all air cells from the last layer, placing it entirely underground','Smooth last air layer?','Yes','No','No');
    switch smth
        case 'Yes'
            ind = find(isair>0,1,'last'); % this ensures it is the lowest elevation air layer
%             [~,ind] = min(isair(isair>0));
            res = str2double(get(H.res,'string'));
            lay = model(:,:,ind);
            lay(lay>1e16) = res;
            lay(isnan(lay)) = res;
            model(:,:,ind) = lay;
            disp(['Layer ',num2str(ind),' air cells replaced with ',num2str(res),' Ohm-m cells'])
            uiwait(msgbox(['Layer ',num2str(ind),' air cells replaced with ',num2str(res),' Ohm-m cells'],'Last Air Layer Smoothed','modal'))
            H.lay = ind;            
        case 'No' 
            disp('Last air layer not smoothed')
            return
    end   
    
    % update the model
    if isfield(H,'AAt') % topography was added in current M3 session
        H.AAt = model;
    elseif sum(sum(isnan(H.AA(:,:,1)))) < H.nx*H.ny % model might be loaded with topography
        H.AA = model;
    end
    
    M3_update_model_plot(H)
    
    guidata(hObject, H);
end