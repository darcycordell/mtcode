function M3_save_ws_mod(hObject, ~, ~)
    % Dec 2019 - double check the file format - model blocks are nx rows by
    % ny columns in file
    H=guidata(hObject);
    
if ~isfield(H,'m') % check to see if model has been generated or loaded        
    warndlg('Make or load a model before you save it!')
    return        
else      
    
    sto=0; sto1=0;
    
    if isfield(H,'AAt') || sum(sum(isnan(H.AA(:,:,1)))) > 0 % model might be loaded with topography
        warndlg('Current version of WSINV3DMT does not support topography!','Model Not Saved')
        return
    end

    AAA=H.AA(1:H.nx-1,1:H.ny-1,1:H.nz-1);
    xx=diff(H.XX).*1000; % this is in m
    yy=diff(H.YY).*1000;
    nxx=H.nx-1;% less 1 because the difference of an array has one less element than the original array
    nyy=H.ny-1;
    nzz=H.nz-1;
    ZZ=diff(H.Z).*1000;
    
    fix_ocean =get(H.fix_ocean,'Value'); % fix ocean resistivity?
    ocean_res =str2double(get(H.ocean_res,'String')); % ocean resistivity

    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1)
    %----------------Model Verification---------------------
    if rem(nxx,2)==1 && rem(nyy,2)==1
        warndlg('model has odd number of columns and rows, model not saved. make them even')
        sto1=1;
    elseif rem(nxx,2)==1
        warndlg('model has odd number of rows, model not saved. make it even')
        sto1=1;
    elseif rem(nyy,2)==1
        warndlg('model has odd number of columns, model not saved. make it even')
        sto1=1;
    end
    for i=1:nxx-1
        stv=find(H.XX(i)<(H.d.x./1000) & (H.d.x./1000)<H.XX(i+1)); % stations in the first vertical mesh cell
        if ~isempty(stv)
            for j=1:H.ny-1
                sty=find(H.YY(j)<(H.d.y(stv)./1000) & (H.d.y(stv)./1000)<H.YY(j+1)); % stations fall in the same cell
                if length(sty)>1
                    hold on;plot(H.d.y(stv(sty))./1000,H.d.x(stv(sty))./1000,'vr','markerfacecolor','r','markersize',12);hold off;
                    sto=1;
                end
            end
        end
    end
    if sto==1
        warndlg('multiple stations in one cell, model not saved')
        return
    end
    ress=unique(AAA);

    if length(ress)>9
        warndlg('You cannot have more than 9 resistivities, model not saved');
        sto1=1;
    end
    %----------------saves initial model file---------------
    if sto==0 && sto1==0
        nnr=length(ress);
        for opi=1:length(ress) % replace model resistivity values with indices
            AAA(find(AAA==ress(opi)))=opi;
        end
        def = {'data_name'};
        prompt = {'enter model file name no extensions'};
        titles  = 'Model file name';
        modf = char(inputdlg(prompt,titles,1,def)); %model file name

        while exist(modf)>0
            newmod_name = questdlg('This file exist! Do you want to overwrite the file?', ...
                'The model file exist', ...
                'Yes', 'No','Yes');
            switch newmod_name
            case 'Yes'
                eval(['delete ',modf]);
                break
            case 'No'
                def = {'data_name'};
                prompt = {'enter a new model file name (no extensions)'};
                titles  = 'Model file name';
                modf = char(inputdlg(prompt,titles,1,def)); %model file name
            end % switch
        end

        %Compare each layer with the next to find depth interfaces (dint)
        aa=1;
        slj=reshape(AAA,nxx*nyy*nzz,1);
        for uu=1:nzz-1
            if slj(aa:uu*nxx*nyy)==slj(uu*nxx*nyy+1:(uu+1)*nxx*nyy)
            else
                dint(uu)=uu+1;
            end
            aa=uu*nxx*nyy+1;
        end
        if exist('dint')==1
            dint=dint(find(dint~=0)); %layer numbers
            dint=[1 dint nzz+1];
        else
            dint=[1 5 nzz+1];
            nnr=1;
            ress=[ress];
        end
        AAA=reshape(slj,nxx,nyy,nzz);
        if 1
            ress(find(ress)<1)=ocean_res;
            if fix_ocean==1
                CMI=AAA; CMI(find(CMI~=1))=0; %everything is zero except ocean CMI=control model index
                %generate control model index ---used to fix some part of the
                %model such as ocean 0 means free 1 means freezed/fixed.
                fid=fopen([modf,'_cmi'],'w');
                %fprintf(fid,'%s\n','# CONTROL MODEL INDEX');
                %this is strange here now--- X and Y will change here to match the
                %program requirements. This is very strange but IS correct.
                fprintf(fid,'%d %d %d\n',nxx,nyy,nzz);
                %For layered model with up to 9 distinct resistivities
                for kj=1:length(dint)-1
                    fprintf(fid,'%d %d\n',dint(kj),dint(kj+1)-1);    %layers
                    for i=nxx:-1:1
                        fprintf(fid,'%d ',CMI(i,:,dint(kj)));
                        if kj==length(dint)-1 && i==1
                        else
                            fprintf(fid,'\n');
                        end
                    end
                end
                fclose(fid);
            end
        end
        %         for pl=1:nz
        %             sld(:,:,pl)=slh(:,:,pl)';
        %         end
        %generate initial model
        if find(ress<1)>0 & fix_ocean==0; warndlg('Your model has resistivity less than 1, you probably have ocean in your model, you may fix ocean during the inversion. Use "Ocean Model" Module');end
        fid=fopen(modf,'w');
        fprintf(fid,'%s\n','# INITIAL MODEL');
        %this is strange here now--- X and Y will change here to match the
        %program requirements. This is very strange but IS correct.
        fprintf(fid,'%d %d %d %d\n',nxx,nyy,nzz,nnr);
        fprintf(fid,'%.0f ',xx);
        fprintf(fid,'\n');
        fprintf(fid,'%.0f ',yy);
        fprintf(fid,'\n');
        fprintf(fid,'%.0f ',ZZ);
        fprintf(fid,'\n');
        fprintf(fid,'%.0f ',ress);
        fprintf(fid,' ');
        fprintf(fid,'\n');
        %For layered model with up to 9 distinct resistivities
        for kj=1:length(dint)-1
            fprintf(fid,'%d %d\n',dint(kj),dint(kj+1)-1);    %layers
            for i=nxx:-1:1
                fprintf(fid,'%d ',AAA(i,:,dint(kj)));
                if kj==length(dint)-1 && i==1
                else
                    fprintf(fid,'\n');
                end
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'%s\n','10. 1. 100. 0.1');
        fprintf(fid,'%s','10. 1. 100. 0.1   ! used in data.3d ');
        fclose(fid);      
    end
end

    guidata(hObject, H);

end % end save_ws_mod_Callback