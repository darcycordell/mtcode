clear all;close all;

% edit March 8 2011,

S_num=51;   % number of stations
F_num=35;   % number of frequencies, can be read in your .te file
rhoerror = 0.05; % Rhoerror =  delta-rho / rho
ssmag = 0.7;    % Static shift magnitude (log10)
str='te_fwd_pred.dat'; %insert your file name here

%----------------------

%read data

m=S_num*F_num;
n=5;

phaerror = rhoerror*100*0.29;   % Pha error in degrees

% Set up statics. Either generate random numbers
% or read in from file (previous random distribution)

log10stat = [-0.2696,    0.1802,    0.5232,    0.1617,    0.6941,    1.4142,   -0.0434, -0.3309,   -0.6985,    0.1245,   -0.0041,   -0.6197,    0.1893,   -0.0647, 0.6880,   -0.2841,   -0.2237,    0.5692,   -0.5376,   -0.0921,    0.6852, -0.7167,    1.0908,   -0.3146,    0.0047,    0.1785,   -0.8789,    0.0119, -1.0148,    1.0568,    0.0064,   -1.3732,   -0.6248,   -0.6847,   -0.9590, -0.4025,   -0.0694,   -0.4560,   -0.3905,    0.7706,   -0.6067,   -0.2067, 0.2495,   -1.0561,   -0.6359,    0.2392, 0.6880, -0.0694, 0.1785, 0.1893, -0.6985]

% These numbers used in inversions

for is =1:S_num
   %log10stat(is) = ssmag*randn(1);
   stat(is) = 10^(log10stat(is));
end

log10stat
hist(log10stat)
pause


M=zeros(m,n);

fid=fopen(str,'r');
stanum=fscanf(fid,'%f',1);

for j=1:S_num
   fnum=fscanf(fid,'%f',1);
   ratio=fscanf(fid,'%f',1);
   j;
   for i=1:F_num
      M(i+(j-1)*F_num,1)=fscanf(fid,'%f',1);
      s=fscanf(fid,'%s',1);
      M(i+(j-1)*F_num,2)=fscanf(fid,'%f',1);
      s=fscanf(fid,'%s',1);
      M(i+(j-1)*F_num,3)=fscanf(fid,'%f',1);
      s=fscanf(fid,'%s',1);
      M(i+(j-1)*F_num,4)=fscanf(fid,'%f',1);
      M(i+(j-1)*F_num,5)=fscanf(fid,'%f',1);
   end
end

fclose(fid);

%add noise

rhonoise=randn(m,1)*rhoerror;
phasenoise=randn(m,1)*phaerror;

T=M(:,1);
Ra=M(:,2)+M(:,2).*rhonoise;
pha=M(:,3)+phasenoise;
err_Ra= zeros(m,1)+rhoerror;
err_pha=zeros(m,1)+phaerror;

for i=1:9:S_num
   str='Apparent Resistivity of Station #00';
   str1='Phase of Station #00';
   s=num2str(i);
      if(i<10)
         str(35:35)=s;
         str1(20:20)=s;
      else
         str(34:35)=s;
         str1(19:20)=s;
      end
   figure(i);
   
   subplot(2,1,1);
   plot(log10(T(1:F_num)),log10(M((i-1)*F_num+1:i*F_num,2)),'r');
   hold on;
   plot(log10(T(1:F_num)),log10(Ra((i-1)*F_num+1:i*F_num)),'b');
   hold on;
   errorbar(log10(T(1:F_num)),log10(Ra((i-1)*F_num+1:i*F_num)),err_Ra(1:F_num));
   hold off;
   title(str);
   axis([-1 5 0 3])

   
   subplot(2,1,2);
   plot(log10(T(1:F_num)),-M((i-1)*F_num+1:i*F_num,3),'r');
   hold on;
   plot(log10(T(1:F_num)),-pha((i-1)*F_num+1:i*F_num),'b');
   hold on;
   errorbar(log10(T(1:F_num)),-pha((i-1)*F_num+1:i*F_num),err_pha(1:F_num));
   hold off;
   title(str1)
   axis([-1 5 0 90])

end


figure(23)
hist(rhonoise)
title('error level')

M(:,2)=Ra;
M(:,3)=pha;
M(:,4)=1.1*err_Ra;
M(:,5)=1.1*err_pha*(3.14/180);

%write data

str='input2.te';    % <<<<<<   insert your output file name here
fid1=fopen(str,'w+');
fprintf(fid1,'%3.0f',S_num);
fprintf(fid1,'\n')

for j=1:S_num
   fprintf(fid1,'%3.0f',F_num);
   fprintf(fid1,' %5.4f',1);
   fprintf(fid1,'\n');
   for i=1:F_num
      fprintf(fid1,'%12.5e',M(i+(j-1)*F_num,1)); % edited to match outputs, 9.3f
      fprintf(fid1,'      (');
      fprintf(fid1,' %7.3f',stat(j)*M(i+(j-1)*F_num,2));
      fprintf(fid1,'     ,');
      fprintf(fid1,' %7.3f',M(i+(j-1)*F_num,3));
      fprintf(fid1,'    )');
      fprintf(fid1,'   %8.4f',M(i+(j-1)*F_num,4));
      fprintf(fid1,'       %8.4f',M(i+(j-1)*F_num,5));
      fprintf(fid1,'\n');
   end
end

fclose(fid1);
