function M3_del_col(hObject,eventdata,~)
    
    H=guidata(hObject);
        
    if ~isfield(H,'m')
       warndlg('Generate a model then modify it using these butttons.')
       return
    end
    
    newXLim = get(H.axes1,'XLim');
    newYLim = get(H.axes1,'YLim');
    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1);
       
    while 1
        [za,~,zc]=ginput(1);
        if zc==3
            break
        end
        xbu=abs(H.YY-za); xbu=find(min(xbu)==xbu);
        H.YY(xbu)=[];
        
        if isfield(H,'AAt')
            H.AAt(:,xbu,:)=[];
        else
            H.AA(:,xbu,:)=[];
        end                
        
        H.nx=length(H.XX); H.ny=length(H.YY);
    
        M3_update_model_plot(H)
        set(H.axes1,'XLim',newXLim)
        set(H.axes1,'YLim',newYLim) 
    end
    
    M3_update_model_plot(H)

    if ~isfield(H,'undo')        
        H.undo=0;    
    end
    H.undo=H.undo+1;

%     M3_save_undo(hObject, eventdata, H)

    guidata(hObject, H)
     
end % end del_col_Callback