function M3_jump_to_layer(hObject,~,~)

    H=guidata(hObject);
    
    if ~isfield(H,'m')
        warndlg('Generate a model then modify it using these butttons.')
    else
        newXLim = get(H.axes1,'XLim');
        newYLim = get(H.axes1,'YLim');
%         view_model_Callback(hObject, eventdata, H)
        
        if str2double(get(H.lay_num,'String')) < 1 || str2double(get(H.lay_num,'String')) > H.nz -1
            warndlg('Enter a valid layer number.')
            set(H.lay_num,'String','1')
            return
        else
        H.lay = str2double(get(H.lay_num,'String'));
        end
        
        M3_update_model_plot(H)
        set(H.axes1,'XLim',newXLim)
        set(H.axes1,'YLim',newYLim)
        if get(H.gridlines,'value')
            shading flat
        else
            shading faceted
        end
    end

    guidata(hObject, H)

end