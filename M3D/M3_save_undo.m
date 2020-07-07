function save_undomodel_Call(hObject, ~, H)

    %If you haven't make the model you can't save it
    sto=0; sto1=0;

    AAA=H.AA(1:H.nx-1,1:H.ny-1,1:H.nz-1);
    xx=diff(H.XX).*1000;% XX and YY have an extra row/column for plotting using pcolor
    yy=diff(H.YY).*1000;
    nxx=H.nx-1;
    nyy=H.ny-1;
    nzz=H.nz-1;
    ZZ=diff(H.Z).*1000;% skip the zero at the beginning of Z

    if length(H.XX)<5
        warndlg('Make a model before you save it!')
        return
    else
        set(H.axes1,'HandleVisibility','ON');
        axes(H.axes1)
        %-------Model Verification needed for backuping too so not everything will be backed up------
        if rem(nxx,2)==1 && rem(nyy,2)==1
            %   warndlg('model has odd number of columns and rows, model not saved. make them even')
            sto1=1;
        elseif rem(nxx,2)==1
            %    warndlg('model has odd number of columns, model not saved. make it even')
            sto1=1;
        elseif rem(nyy,2)==1
            %     warndlg('model has odd number of rows, model not saved. make it even')
            sto1=1;
        end

        for i=1:nxx-1
            stv=find(xx(i)<H.x & H.x<xx(i+1));%  stations in the first vertical mesh cell
            if isempty(stv)==0
                for j=1:H.ny-1
                    sty=find(yy(j)<H.y(stv) & H.y(stv)<yy(j+1)); % stations fall in the same cell
                    if length(sty)>1
                        hold on;plot(H.x(stv(sty)),H.y(stv(sty)),'vr','markerfacecolor','r','markersize',12);hold off;
                        sto=1;
                    end
                end
            end
        end
        ress=unique(AAA);
        %----------------saves initial model file---------------
        if sto==0 && sto1==0
            nnr=length(ress);
            for opi=1:length(ress)
                AAA(find(AAA==ress(opi)))=opi;
            end
            modf = ['undo',num2str(H.undo)]; %model file name
            aa=1;slj=reshape(AAA,nxx*nyy*nzz,1);
            for uu=1:nzz-1
                if slj(aa:uu*nxx*nyy)==slj(uu*nxx*nyy+1:(uu+1)*nxx*nyy)
                else
                    dint(uu)=uu+1;
                end
                aa=uu*nxx*nyy+1;
            end
            if exist('dint')==1
                dint=dint(find(dint~=0));dint=[1 dint nzz+1];
            else
                dint=[1 5 nzz+1];nnr=length(ress);%ress=[ress];
            end
            AAA=reshape(slj,nxx,nyy,nzz);
            fid=fopen(modf,'w');
            fprintf(fid,'%s\n','# INITIAL MODEL');
            %this is strange here now--- X and Y will change here to match the
            %program requirements. This is very strange but IS correct.
            fprintf(fid,'%d %d %d %d\n',nyy,nxx,nzz,nnr);
            fprintf(fid,'%.0f ',yy);
            fprintf(fid,'\n');
            fprintf(fid,'%.0f ',xx);
            fprintf(fid,'\n');
            fprintf(fid,'%.0f ',ZZ);
            fprintf(fid,'\n');
            fprintf(fid,'%.0f ',ress);
            fprintf(fid,' ');
            fprintf(fid,'\n');
            %For layered model with up to 9 distinct resistivities
            for kj=1:length(dint)-1
                fprintf(fid,'%d %d\n',dint(kj),dint(kj+1)-1);    %layers
                for i=nyy:-1:1
                    fprintf(fid,'%d ',AAA(:,i,dint(kj)));
                    if kj==length(dint)-1 && i==1
                    else
                        fprintf(fid,'\n');
                    end
                end
            end
            fclose(fid);

            guidata(hObject, H);

        end
    end

end % end save undo model