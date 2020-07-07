function [lstatic,statics]=load_stat(kmode,cmode,nsites,indsit,indsitmod,...
                                     outpath)
% what data are inverted for static shift? 
% read statics.tx ==> lstatic
% read tx.dat     ==> statics

disp('IN : load_stat')

global statics lstatic
kmode;
for imode=1:2
    if kmode(imode)==1
        % which mode & station has been inverted for static shift?
        name=[outpath,'statics.',lower(char(cmode(imode)))];
        if exist(name,'file')
            ltstatic=load(name);
            for isite=1:length(ltstatic)
                ksite=find(indsit==indsitmod(imode,isite));
                lstatic(imode,ksite)=ltstatic(isite);
            end
            disp(['Found: ',name]);
        else
            % no statics.tx file found 
            lstatic(imode,1:nsites)=ones(1,nsites);
            disp(['Not found: ',name]);
            disp(['Set lstatic(',imode,',:) to ones']);
        end
        disp(imode)
        if sum(lstatic(imode,:))>0
            % what are the actual static shift values?
            name=[outpath,lower(char(cmode(imode))),'_pred.dat'];
            if exist(name,'file')
                disp(['Found: ',name]);
                fid=fopen(name,'r');
                line=1;
                isite=0;
                while line ~= -1
                    line=fgetl(fid);
                    k=findstr(line,'.');
                    if length(k)==1
                        isite=isite+1;
                        temp = sscanf(line,'%f');
                        ksite=find(indsit==indsitmod(imode,isite));
                        % statics found are dat(pred)/dat(obs) => log10
                        statics(imode,ksite)=log10(temp(2));
                    end
                end
                fclose(fid);
            else
                disp(['Not found: ',name]);
                disp(['Set statics(',imode,',:) to zeros']);
                statics(imode,ksite)=zeros(nsites);
            end
        end
    end
end
end