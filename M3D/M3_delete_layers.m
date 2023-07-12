function M3_delete_layers(hObject, ~, ~)
    H=guidata(hObject);
    
    if isfield(H,'AAt')
        isearth = squeeze(sum(sum(H.AAt<1e16,1),2)); % find earth cells per layer
    elseif isfield(H,'AA')
        isearth = squeeze(sum(sum(H.AA<1e16,1),2)); % find earth cells per layer
    else
        warndlg('No mesh found','Error')
    end
    
    figure
    plot(isearth,'ko'); hold on; plot(isearth,'k-')
    xlabel('Layer number')
    ylabel('Number of non-air cells')        

    list = cellstr(cat(2,repmat('Layer ',[H.nz-1 1]),num2str((1:H.nz-1)'),repmat(': ',[H.nz-1 1]),num2str(H.Z(1:end-1)'))); % H.nz-1 since M3 counts all edges

    [sel,val] = listdlg('PromptString','Select layers to delete (depth in km)','ListString',list,'Name','Delete Layers','ListSize',[300 300]);

    if ~val
        disp('No layers selected. Returning to Main Menu')

    else % layers selected to be deleted
        
        flag = 0;
        for i = 1:length(sel) % check for stations in selected layers
            if sum(H.d.z./1000 >= H.Z(sel(i))) > 0 && sum(H.d.z./1000 <= H.Z(sel(i)+1)) > 0 % if a station is in a selected layer
                flag = 1;
            end
        end
        
        if flag
            answer = questdlg('WARNING: Some stations are located in layers to be deleted. Stations elevations will be moved. Proceed?','WARNING','Yes','No','No');
            if ~strcmp(answer,'Yes')
                return
            end
        end

        uiwait(msgbox([num2str(length(sel)),' layers removed from model'],'Layers Deleted!','modal'))
        disp([num2str(length(sel)),' layers removed from model'])
        if isfield(H,'AAt')
            H.AAt(:,:,sel) = [];
            H.Zsurf = H.Zsurf - length(sel);
            H.Zsurf(H.Zsurf<=0) = 1;
        elseif isfield(H,'AA')
            H.AA(:,:,sel) = [];
        end
        H.Z(sel) = [];
        H.top = H.Z(1)*1000;   % in m               
        H.nz = H.nz-length(sel);
        M3_update_model_plot(H)

    end
        
    guidata(hObject, H);
end