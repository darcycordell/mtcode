function [seed,kmode,modfile,datfile,tau,err_floor,statshift,indsitmod]=...
              load_par_iso(name)
%==========================================================================
% read inversion parameter file
clear datfile err_floor 
indsitmod=zeros(3,1); err_floor=zeros(3,2);

fid=fopen(name,'r');

line=fgetl(fid); seed=sscanf(line,'%s',1);
line=fgetl(fid); kmode(1)=sscanf(line,'%s',1);
line=fgetl(fid); kmode(2)=sscanf(line,'%s',1);
line=fgetl(fid); kmode(3)=sscanf(line,'%s',1);
kmode=kmode=='y';
line=fgetl(fid); modfile=sscanf(line,'%s',1);

% Next line determines if inversion is isotropic or anisotropic, DR 23/12/2010
line=fgetl(fid);
devmod=sscanf(line,'%d',1);

modind=find(kmode>0);
for imode=1:length(modind)
    line=fgetl(fid); datfile{modind(imode)}=sscanf(line,'%s',1);
end
if kmode(3)==1 line=fgetl(fid); conj=sscanf(line,'%s',1); end
line=fgetl(fid); uni_std_lap=sscanf(line,'%d',1);
line=fgetl(fid); grd_lap_reg=sscanf(line,'%d',1);
line=fgetl(fid); tau=sscanf(line,'%f',1);            %MJU  2016-03-18 Was not reading  if tau less than 1

line=fgetl(fid); tear=sscanf(line,'%s',1);
if (tear=='y')                                       % MJU 2016-04-02
    line=fgetl(fid); tearfile=sscanf(line,'%s',1);   % MJU 2016-04-02
end                                                  % MJU 2016-04-02


line=fgetl(fid); sensout=sscanf(line,'%s',1);
if kmode(1)==1 line=fgetl(fid); err_floor(1,1)=sscanf(line,'%f',1)'; 
               line=fgetl(fid); err_floor(1,2)=sscanf(line,'%f',1)';
end
if kmode(2)==1 line=fgetl(fid); err_floor(2,1)=sscanf(line,'%f',1)'; 
               line=fgetl(fid); err_floor(2,2)=sscanf(line,'%f',1)';
end
if kmode(3)==1 line=fgetl(fid); err_floor(3,1)=sscanf(line,'%f',1); 
                                err_floor(3,2)=err_floor(3,1);
end
line=fgetl(fid); parfix=sscanf(line,'%s',1);
if (parfix=='y') 
    line=fgetl(fid); fixfile=sscanf(line,'%s',1);
    line=fgetl(fid); fixdamp=sscanf(line,'%d',1);
end

if (kmode(1)==1 | kmode(2)==1);
    line=fgetl(fid); statshift=sscanf(line,'%s',1);
    if (statshift=='y') 
        line=fgetl(fid); statvar=sscanf(line,'%d',1);
        line=fgetl(fid); statdamp=sscanf(line,'%d',1);
    end
else
    statshift='n';
end
line=fgetl(fid); maxit=sscanf(line,'%d',1);

imode=0; inodep=1E6;
while 1
    line=fgets(fid);
    if line==-1 break; end
    if strfind(line,'END') break; end
    inode=str2num(line);
    if inode<inodep
        imode=imode+1; 
        isite=1;
    else
        isite=isite+1;
    end
    indsitmod(modind(imode),isite)=inode;
    inodep=inode;
end

fclose(fid);
end
