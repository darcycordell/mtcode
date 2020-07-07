function M3_min_res(hObject,~)
    H=guidata(gcbo);
    set(H.min_res,'String',get(hObject,'String')) 
    if isfield(H,'XX')
        M3_update_model_plot(H)
    end
    guidata(hObject, H); 
end