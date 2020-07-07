function M3_view_model(hObject,~,~)
    
    H=guidata(hObject);
       
    if ~isfield(H,'m')
       warndlg('Generate a model then view it.')
    else
       M3_update_model_plot(H)
    end
    
end % end view_model_Callback