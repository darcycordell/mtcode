function [err]=load_err(err,name,imode,indsit,indsitmod,per)
%==========================================================================
% reads nlcg input data:
%
% data error matrix from input data
% rho errors:   ln(rho)
% phase errors: radians
%

disp('IN : load_err')

fid=fopen(name,'r');

line=1;
isite=0;
while line ~= -1
    line=fgetl(fid);
    k=findstr(line,'.');
    j=length(k);
    if j>2
        if next_site==1
            next_site=0;
            % isite: index for this site
            isite=isite+1; iper=0;
            % ksite: overall site index
            ksite=find(indsit==indsitmod(imode,isite));
            % was data at this station found in rsp file?
        end
        iper=iper+1;
        % get correct period index
        tfreq=str2num(char(strread(line,'%s',1)));
%         for iper=iper:length(per) % Commented out, DR, 10/12/2010
%             if (tfreq/per(iper)-1 < 0.001) break; end % Commented out, DR, 10/12/2010
%         end % Commented out, DR, 10/12/2010
        tchar=line(findstr(line,')')+1:length(line));
        temp = sscanf(tchar,'%f');
        if (imode==3) temp(2)=temp(1); end
        err(imode,ksite,iper,:) = temp;
    else
        next_site=1;
    end
end

fclose(fid);
end