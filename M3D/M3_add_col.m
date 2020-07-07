function M3_add_col(hObject,eventdata,~)
    
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
        za=round(za*1000)/1000;minza=min(abs(abs(H.YY)-abs(za)));
        if minza<1   %minimum mesh spacing is 1km
            button = questdlg('Columns are too close! Do you want to add anyways?','Are you sure','Yes','No','Yes');
            switch button
                case 'No'
                    return
            end
        end
        [yut, ~]=size(H.YY);
        if yut>1
            H.YY=sort([H.YY;za]);
        else
            H.YY=sort([H.YY,za]);
        end
        lk=find(H.YY==za);
        H.nx=length(H.XX); H.ny=length(H.YY);
        if isfield(H,'AAt')
            ab(:,1:lk,:)=H.AAt(:,1:lk,:); 
            ab(:,lk+1,:)=H.AAt(:,lk,:); 
            ab(:,lk+2:H.ny,:)=H.AAt(:,lk+1:H.ny-1,:); 
            H.AAt=ab;            
        else           
            ab(:,1:lk,:)=H.AA(:,1:lk,:); 
            ab(:,lk+1,:)=H.AA(:,lk,:); 
            ab(:,lk+2:H.ny,:)=H.AA(:,lk+1:H.ny-1,:); 
            H.AA=ab;
        end
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

end % end add_col_Callback