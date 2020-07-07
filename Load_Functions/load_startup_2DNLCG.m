function [par]=load_startup_2DNLCG(mode,name)
% Function which reads the Mackie NLCG 2D Inversion Parameter file. This
% file must have a very specified format otherwise it will not work.
%
% Profile_Name
% TM mode (y/n)
% TE mode (y/n)
% Tipper mode (y/n)
% Input model file
% Regularization (1 or 2)
% Input model files
%   *.tm
%   *.te
%   *.tip
% Conjugate tipper? (y/n)
% Laplacian regularization (1 or 2)
% Regularization order (1 or 2)
% Tau parameter
% Output sensitivity map? (y/n)
%   If y, then must include *.tear file here
% Save intermediate steps? (y/n)
% Error floors
%   TM rho error percent
%   TM phase error percent
%   TE rho error percent
%   TE phase error percent
%   Tipper error (absolute)
% Fixed parameters (y/n)
%   If y, then must include *.fix file here
%   If y, include tau for fixed cells (e.g. 10^6)
% Static shifts (y/n)
%   If y, set static shift variance (e.g. 10)
%   If y, set static shift damping (e.g. 10000)
% Number of iterations (e.g. 200)
%
% The remainder of the par file specifies the cell columns where stations
% are located.


fid=fopen(name,'r');

disp('Loading 2D NLCG Inversion Parameter File')

par.profile_name = fgetl(fid);

for i = 1:3
    line = fgetl(fid);
    if strcmp(line,'y')==1
        par.mode(i) = 1;
    else
        par.mode(i) = 0;
    end
end

if strcmp(mode,'iso') % 2019 for aniso compatibility
    par.m0_file = fgetl(fid);
else
    par.m0_xx_file = fgetl(fid);
    par.m0_yy_file = fgetl(fid);
    par.m0_zz_file = fgetl(fid);
end
par.regularization = str2double(fgetl(fid));

par.data_files = cell(3,1);
if par.mode(1) == 1
    par.data_files{1} = fgetl(fid);
end

if par.mode(2) == 1
    par.data_files{2} = fgetl(fid);
end

if par.mode(3) == 1
    par.data_files{3} = fgetl(fid);
    par.conjugate = fgetl(fid);
end

par.laplacian = str2double(fgetl(fid));
par.regularization_order = str2double(fgetl(fid));
par.tau = str2double(fgetl(fid));
if strcmp(mode,'aniso')
    par.aniso_tau = str2double(fgetl(fid));
end
par.sensitivity = fgetl(fid);

if strcmp(par.sensitivity,'y') == 1
    par.tear_file = fgetl(fid);
end

par.save_steps = fgetl(fid);

if par.mode(1) == 1
    par.tm_errflr(1) = str2double(fgetl(fid));
    par.tm_errflr(2) = str2double(fgetl(fid));
end

if par.mode(2) == 1
    par.te_errflr(1) = str2double(fgetl(fid));
    par.te_errflr(2) = str2double(fgetl(fid));
end

if par.mode(3) == 1
    par.tip_errflr = str2double(fgetl(fid));
end

par.fixed_parameters = fgetl(fid);

if strcmp(par.fixed_parameters,'y')
    par.fix_file = fgetl(fid);
    par.fix_tau = str2double(fgetl(fid));
end

par.statics = fgetl(fid);

if strcmp(par.statics,'y')==1
    par.statics_variance = str2double(fgetl(fid));
    par.statics_damping = str2double(fgetl(fid));
end

par.niter = str2double(fgetl(fid));

imode=0; inodep=1E6;
while 1
    line=fgets(fid);
    if line==-1; break; end
    if strfind(line,'END'); break; end
    
    inode=str2double(line);
    if inode<inodep
        imode=imode+1; 
        isite=1;
    else
        isite=isite+1;
    end
    indsitmod(imode,isite)=inode;
    inodep=inode;
end



par.allsites = unique(indsitmod(:));
par.allsites(par.allsites==0)=[];

sites = zeros(3,length(par.allsites));
if par.mode(1) == 1
    sites(1,:) = indsitmod(1,:);
    if par.mode(2) == 1
        sites(2,:) = indsitmod(2,:);
        if par.mode(3) == 1
            sites(3,:) = indsitmod(3,:);
        end
    else
        if par.mode(3) == 1
            sites(3,:) = indsitmod(2,:);
        end
    end
else
    if par.mode(2) == 1
        sites(2,:) = indsitmod(1,:);
        if par.mode(3) == 1
            sites(3,:) = indsitmod(2,:);
        end
    else
        if par.mode(3) == 1
            sites(3,:) = par.indsitemod(1,:);
        end
    end
end

par.sites = sites;



end
