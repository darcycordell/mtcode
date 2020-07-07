function M3_save_modem_mod(hObject,~,~)

    H=guidata(hObject);
    
    if ~isfield(H,'m') % check to see if model has been generated or loaded        
        warndlg('Make or load a model before you save it!')
        return        
    else       
        sto=0; 
        
        % Convert from depths (including 0 and depth to bottom of mesh) to
        % thicknesses. This removes the lowest layer from the res. matrix
        if isfield(H,'AAt')
            AAA=H.AAt(1:H.nx-1,1:H.ny-1,1:H.nz-1); 
        else
            AAA=H.AA(1:H.nx-1,1:H.ny-1,1:H.nz-1); 
        end
        xx=diff(H.XX).*1000; % this is in m
        yy=diff(H.YY).*1000;
        nxx=H.nx-1;% less 1 because the diff of an array has one less element than the original array
        nyy=H.ny-1;
        nzz=H.nz-1;
        ZZ=diff(H.Z).*1000;
     

        set(H.axes1,'HandleVisibility','ON');
        axes(H.axes1)
        %----------------Model Verification---------------------
        for i=1:nxx-1
           stv=find(H.XX(i)<(H.d.x./1000) & (H.d.x./1000)<H.XX(i+1));%  stations in the first vertical mesh cell
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
               
        % below slightly rewritten in the unlikely case that a res index is
        % mistaken for a res value when writing matrix of res indices!
        
        % Commented out Dec 2019 - only practical for halfspace / simple
        % models, and currently we do not save models with resistivity
        % indices
%         ress=unique(AAA);
%         num_ind = zeros(length(ress),1);
%         temp_ind = [];
%         AAA_ind = AAA; % intitialize matrix of res indices
%         for opi=1:length(ress) % loop through number of unique resistivities
%            temp = find(AAA==ress(opi)); % find indices of all entries matching unqiue res
%            num_ind(opi) = length(temp); % store number of indices found per res value
%            temp_ind = [temp_ind; temp]; % add onto end of vector of all indices
%         end
%         num_ind = [0;num_ind];
%         cnum_ind = cumsum(num_ind);
%         for rpi = 1:length(ress) % loop through number of unique resistivities
%            AAA_ind(temp_ind(cnum_ind(rpi)+1:cnum_ind(rpi+1))) = rpi; % replace res value with index
%         end
% 
%         H.AAA_ind = AAA_ind;

        def = {'data_name'};
        prompt = {'enter ModEM model file name'};
        titles  = 'Model file name';
        modf = char(inputdlg(prompt,titles,1,def)); %model file name

        while exist(modf)>0
            newmod_name = questdlg('This file exists! Do you want to overwrite the file?', ...
                'The model file exists', ...
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
                        
        origin = round([-sum(xx)/2 -sum(yy)/2 min(H.Z)*1000]); % numair = 0 with no topo
%         rotation = 0; % this is assumed to be zero...
        write_model_modem(modf,round(xx),round(yy),round(ZZ),AAA,origin)
    end

    guidata(hObject, H);

end % end save_modem_mod_Callback