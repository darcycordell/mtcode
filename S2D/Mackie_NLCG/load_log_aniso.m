function [objfun]=load_log_aniso(name,err_floor)
%==========================================================================
% get used error_floor and objective function
% objf(iter,1:6)
% 1            2                         3       4                5               6
% s1 [=chisq], rms [=sqrt(chisq/ndpts)], s2/tau, s2 [=roughness], s3[=closeness], s [=objfunc]

disp('IN : load_log_aniso')

global version
fid=fopen(name,'r');
while 1
    line=fgetl(fid);
    if line==-1 break; end
    word=strread(line,'%s');
    if strfind(char(word(1)),'iter')
        if version == '6_11'
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            fgetl(fid); % Only for anisotropic inversion, DR 15/04/2010
            fgetl(fid); % Only for anisotropic inversion, DR 15/04/2010
            fgetl(fid); % Only for anisotropic inversion, DR 15/04/2010
        end
        while 1
            line=fgetl(fid);
            if line==-1 
                fclose(fid);
                return;
            end
            [temp]=sscanf(line,'%f');
            objfun(temp(1)+1,1:6)=temp([2:3 7:9 11])'; % Equivalent sensitivity parameters for anisotropic inversion, DR 15/04/2010
        end
    end
end