function M3_next_layer(hObject,~,~)
    
    H=guidata(hObject);
    
    if ~isfield(H,'m')
        warndlg('Generate a model then modify it using these butttons.')
    else
        newXLim = get(H.axes1,'XLim');
        newYLim = get(H.axes1,'YLim');
%         view_model_Callback(hObject, eventdata, H)
        
        if H.lay<H.nz-1
            H.lay=H.lay+1;
        else
            warndlg('Invalid layer number.')
            H.lay=1;
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