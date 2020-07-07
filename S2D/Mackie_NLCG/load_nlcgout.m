function [dat,per,indsit,sitrms]=load_nlcgout(name,kmode,indsitmod)
%=============================================================================
% this function reads nlcg output.data into a matrix dat;
% dat matrix has 13 colomuns :
%       obsTMrho  obsTMphs    calTMrho  calTMphs ...
%   ... obsTErho  obsTEphs    calTErho  calTEphs ...
%   ... obsHzreal obHzimag    calHzreal calHzimag     Period
% there will be freq*sites rows

disp('IN : load_nlcgout')

fid=fopen(name,'r');
line=1; isite=0;
while line ~= -1
    line=fgetl(fid);
    if findstr(line,'block')
        isite=isite+1;
        word=strread(line,'%s');
        indsit(isite)=str2num(char(word(4)));
        line=fgetl(fid); line=fgetl(fid);
        word=strread(line,'%s');
        sitrms(isite)=str2num(char(word(7)));
    end
    if findstr(line,'===')
        iper=0;
        % which modes for this site (indsitmod from parameter file)
      %smode=sum(indsitmod==indsit(isite),2)>0;
        % results in indices:
      %dind=[reshape((smode*ones(1,4).*[1:4; 5:8; 9:12])',1,[]) 13];
        % only indices>0:
        %dind=find(dind>0);
        % start of data section of site isite
        while 1
            line=fgetl(fid);
            if ~(length(line)>30) break; end
            iper=iper+1;
            word=strread(line,'%s');
            for iword=1:length(word)
                if findstr(char(word(iword)),'*') word{iword}='NaN'; end
            end
      %dat(isite,iper,dind)=str2num(char(word(dind)));
            dat(isite,iper,:)=str2num(char(word));
        end
    end
end

% period
per=sort(reshape(dat(:,:,13),1,[]));
per=per(find(per>0));
per=[per(diff(per)>0), per(length(per))];

% adjust data indices according to period (not sure if ever needed)

% if you have repeated periods (0.03 0.03 for example) in rsp file then you
% need this step to remove the data corresponding to the first one in the
% two repeated periods; debugged by Enci Wang Mar 16, 2016
for isite=1:size(dat,1)
    tper=squeeze(dat(isite,:,13));
    dat_clean(isite,:,13)=per;
    tdat=squeeze(dat(isite,:,1:12));
    % set data to NaN and rewrite
    dat_clean(isite,:,1:12)=NaN;
    for iper=1:length(per)
        for itper=1:length(tper)
            if tper(itper)==per(iper)
                dat_clean(isite,iper,1:12)=tdat(itper,:);
            end
        end
    end
end
dat=dat_clean;
fclose(fid);
end