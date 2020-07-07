function get_freq
% 3 July 2007, get_freq_w2_E2 This version can handle with different types of the *.fix files used in external inversion
%26March2004 latest modification for realtime freqset modification
%Changes applied to all files separately te,tm,tip

clear all;
modec=menu('choose mode','te','tm','tip')
mode  = {'te';'tm';'tip'};
[datfile,datpath] = uigetfile('*.te;*.tm;*.tip','pick any data file');
disp('data files are supposed to have identical name root');
nroot = [datpath,datfile(1:max(find(datfile=='.'))-1)];
C = zeros(5,100,100,10);
mode_found = [0,0,0];
for imode = modec:modec
    fid_in=fopen([nroot,'.',mode{imode}],'r');
    if (fid_in == -1) 
        continue;
    else
        disp([' input file >',[nroot,'.',mode{imode}],'< opened']); 
        mode_found(imode) = 1;
    end
    nsites(imode)=fscanf(fid_in,'%i',1);
    nper = fscanf(fid_in,'%f',1);
    if (imode<3) 
        static(imode,1) = fscanf(fid_in,'%f',1);
    end
    for isite=1:nsites(imode)
        for iper=1:nper;
            C(imode,isite,iper,1)=fscanf(fid_in,'%f',1);
            C(imode,isite,iper,6)=fscanf(fid_in,'%s',1);    % string will be omitted
            C(imode,isite,iper,2)=fscanf(fid_in,'%f',1);
            C(imode,isite,iper,6)=fscanf(fid_in,'%s',1);    % string will be omitted
            C(imode,isite,iper,3)=fscanf(fid_in,'%f',1);
            C(imode,isite,iper,6)=fscanf(fid_in,'%s',1);    % string will be omitted
            C(imode,isite,iper,4)=fscanf(fid_in,'%f',1);
            if (imode<3)
                C(imode,isite,iper,5)=fscanf(fid_in,'%f',1);
            end
        end
        if isite<nsites(imode)
            nper=fscanf(fid_in,'%f',1);
            if (imode<3)
                static(imode,isite+1)=fscanf(fid_in,'%f',1);
            end
        end
    end
    fclose(fid_in);
end
% find a set of periods for all modes
ifound=0; xtol=0.22;
for imode = modec:modec
    if mode_found(imode)==1
        for isite=1:nsites(imode)
            xpers=squeeze(C(imode,isite,:,1));
            % sparse -> contains no zero elements
            xperind=find(sparse(xpers));
            xper=xpers(xperind);
            for ii=1:length(xper)
                if ifound==0
                    ifound=1;
                    xperiod(ifound)=xper(ii);
                else
                    inew=new_freq(xper(ii),xperiod,xtol,ifound);
                    if inew==1
                        ifound=ifound+1;
                        xperiod(ifound)=xper(ii);
                    end  
                end
            end
        end
    end
end
% output periods
xperiod=sort(xperiod');
cmodes = '';
for imode = modec:modec
    if mode_found(imode)==1
        cmodes = [cmodes char(mode{imode})];
    end
end
save(['freqset_',cmodes,'.txt'],'xperiod','-ascii')
         nexttodo=menu('do you want to edit frequency set','edit freq_set','do not edit','keep going'); %Edit freqset
         if nexttodo==1
             open (['freqset_',cmodes,'.txt'])
             nexttodo=menu('do you want to edit frequency set','edit freq_set','do not edit','keep going'); %Realtime Edit freqset
             fid_inn=fopen(['freqset_',cmodes,'.txt'],'r');
             xperiod=fscanf(fid_inn,'%f',length(xperiod));
             fclose(fid_inn);
         end
nperiods=length(xperiod);
% write data file at new periods
for imode = 1:3
    if mode_found(imode)==1
        % indices of missing data
        nodata = ones(2,nsites(imode),nperiods);
        fid_out = fopen([nroot,'.input.',mode{imode}],'w');
        disp([' output file: >',[nroot,'.input.',mode{imode}],'<']); 
        fprintf(fid_out,'%i\n',nsites(imode));
        for isite=1:nsites(imode)
            if (imode < 3)
                fprintf(fid_out,'%i  %f\n',[nperiods static(imode,isite)]);
            else
                fprintf(fid_out,'%i\n',nperiods);
            end
            idum=1:nperiods;
            if (imode<3)
                mtdata(idum,1) = 100;
                mtdata(idum,2) = -45;
            else
                mtdata(idum,1) = 0;
                mtdata(idum,2) = 0;
            end
            mtdata(idum,3) = 1000000;
            mtdata(idum,4) = 1000000;
            xpers=squeeze(C(imode,isite,:,1));
            xperind=find(sparse(xpers));
            xper=xpers(xperind);
            % imode <3 (te,tm): 1,2 -> rhoa, phi; 3,4 error rhoa,phi;
            % imode =3 (tp):    1,2 -> re/im tp;  3 -> error tp;        
            for nn=1:length(xper)
                [iper,ifound]=which_period(xperiod,xper(nn),nperiods);
                if ifound==1
                    mtdata(iper,1) = C(imode,isite,nn,2);
                    mtdata(iper,2) = C(imode,isite,nn,3);
                    mtdata(iper,3) = C(imode,isite,nn,4);
                    if (imode < 3)
                        mtdata(iper,4) = C(imode,isite,nn,5);
                    end
                end            
            end
            a='     (  '; b=',      '; u=')     '; q='    ';
            % where is data missing
            if imode<3
                nofreq = find(mtdata(:,3) == 1000000);
                nodata(1,isite,nofreq)=0;
                nofreq = find(mtdata(:,4) == 1000000);
                nodata(2,isite,nofreq)=0;
            else
                nofreq = find(mtdata(:,3) == 1000000);
                nodata(1,isite,nofreq)=0;
                nodata(2,isite,nofreq)=0;
            end    
            for iper=1:nperiods
                fprintf(fid_out,'%e',xperiod(iper));
                fprintf(fid_out,'%s',a);
                fprintf(fid_out,'%e',mtdata(iper,1));
                fprintf(fid_out,'%s',b);
                fprintf(fid_out,'%e',mtdata(iper,2));
                fprintf(fid_out,'%s',u);
                if (imode < 3)
                    fprintf(fid_out,'%e',mtdata(iper,3));
                    fprintf(fid_out,'%s',q);
                    fprintf(fid_out,'%e\n',mtdata(iper,4));
                else
                    fprintf(fid_out,'%e\n',mtdata(iper,3));
                end
            end
        end
        fclose(fid_out);
        % write out missing data matrix
        fid_nodata = fopen([nroot,'.index_',mode{imode}],'w');
        if imode<3
            fprintf(fid_nodata,'%s\n','% apparent resistivities');
        else
            fprintf(fid_nodata,'%s\n','% real parts');
        end       
        for iperiod=1:nperiods
            for isite=1:nsites(imode)
                fprintf(fid_nodata,'%i',nodata(1,isite,iperiod));
            end
            fprintf(fid_nodata,'\n',' ');    
        end
        if imode<3
            fprintf(fid_nodata,'%s\n','% phases');
        else
            fprintf(fid_nodata,'%s\n','% imaginary parts');
        end   
        for iperiod=1:nperiods
            for isite=1:nsites(imode)
                fprintf(fid_nodata,'%i',nodata(2,isite,iperiod));
            end
            fprintf(fid_nodata,'\n',' ');    
        end
        fclose(fid_nodata);
    end
end
%change format of the *.fix file
fixit(nroot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
function inew=new_freq(xper,xperiod,xtol,ifound)
inew=1;
for j=1:ifound
    diff=abs(xper-xperiod(j))/xperiod(j);
    if diff<xtol
        inew=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [iper,ifound]=which_period(xperiod,xper,nperiods)
tolerance=0.22;
ifound=0;
iper=0;
for j=1:nperiods
    diff=abs((xper-xperiod(j))/xper);
    if diff<tolerance
        iper=j;
        ifound=1;
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function fixit(nroot)
%Ersan Türkoglu, 19 June 2007
%This program looks for a *.fix file in the folder...
%if there is, reads the new format file and writes it for
%external inversion format.
root=[nroot,'.fix'];lr=length(root);slashes=find(root=='\'); lsla=length(slashes);
root=root(slashes(lsla)+1:length(root));
%root='1_1.fix'
warning off
if exist(root)==2
    fid=fopen(root,'r');
    er=fgetl(fid);
    ex=length(er);
    fclose(fid);
    fopen(root,'r');
    s=1;
    z=find(er==',');
    if ex>30
        while 1
            test(s,:)=fgetl(fid);
            if test(s,1)=='0' | test(s,1)=='1'
            else
                break
            end
            s=s+1;
        end
        test=test(1:s,:);
        [i,j]=find(test=='1');
        outputs=[i,j];
        fclose(fid);
        fid=fopen(root,'w');
        fprintf(fid,'%5d %5d\n',outputs');
        fclose(fid);
    elseif z>0
        root
        if isempty(str2num(root(1)))==1
            eval(['load ' root]);
            rost=root(1:length(root)-4);
            rost=round(eval(rost));
            fid=fopen(root,'w');
            fprintf(fid,'%5.0d %5.0d\n',rost');
            fclose(fid);
            %eval(['save -ascii ' root ' rost'])
        else
            eval(['load ' root]);
            rost=['X' root(1:length(root)-4)];
            rost=round(eval(rost));
            fid=fopen(root,'w');
            fprintf(fid,'%5.0d %5.0d\n',rost');
            fclose(fid);
            %eval(['save -ascii ' root ' rost'])
        end
    else
        'your *.fix file is good'
    end
else
    'your model is free to change, no *.fix file found!'
end
