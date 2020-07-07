clear all;close all;

%read data

S_num=61;
F_num=25;
m=S_num*F_num;
n=5;

% Rhoerror =  delta-rho / rho
rhoerror = 0.05
% Pha error in degrees
phaerror = rhoerror*100*0.29


M=zeros(m,n);
str='pr1_prisms61.tm';
fid=fopen(str,'r');
stanum=fscanf(fid,'%f',1)

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
   axis([-1 5 0 4])

   
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

str='pr1_prisms61.noise.tm';
fid1=fopen(str,'w+');
fprintf(fid1,'%3.0f',S_num);
fprintf(fid1,'\n')

for j=1:S_num
   fprintf(fid1,'%3.0f',F_num);
   fprintf(fid1,' %5.4f',1);
   fprintf(fid1,'\n');
   for i=1:F_num
      fprintf(fid1,'%11.5f',M(i+(j-1)*F_num,1));
      fprintf(fid1,'  (');
      fprintf(fid1,' %10.3f',M(i+(j-1)*F_num,2));
      fprintf(fid1,'  ,');
      fprintf(fid1,' %7.2f',M(i+(j-1)*F_num,3));
      fprintf(fid1,'  )');
      fprintf(fid1,'  %6.4f',M(i+(j-1)*F_num,4));
      fprintf(fid1,'  %6.4f',M(i+(j-1)*F_num,5));
      fprintf(fid1,'\n');
   end
end

fclose(fid1);
