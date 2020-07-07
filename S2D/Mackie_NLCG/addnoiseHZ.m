clear all;close all;

%read data

S_num=61;
F_num=25;
m=S_num*F_num;
n=4;

noise_level=0.02;

M=zeros(m,n);
str='pr1_prisms61.tip';
fid=fopen(str,'r');
line=fgetl(fid)

for j=1:S_num
   line=fgetl(fid)
   for i=1:F_num
      line=fgetl(fid);      
      M(i+(j-1)*F_num,1)=  sscanf(line(1:12),'%f');
      M(i+(j-1)*F_num,2)=  sscanf(line(15:25),'%f');
      M(i+(j-1)*F_num,3)=  sscanf(line(29:39),'%f');
      M(i+(j-1)*F_num,4)=  sscanf(line(44:50),'%f');
   end
end

fclose(fid);

%add noise

for i=1:m
   T(i)=M(i,1);
   real_noise(i) =  randn(1)*noise_level;
   imag_noise(i) =  randn(1)*noise_level;

   realHz (i)=M(i,2)+real_noise(i);
   imagHz (i)=M(i,3)+imag_noise(i);
   error  (i)= 1.1*noise_level;
end

figure(23)
subplot(2,1,1)
hist(real_noise)
title('error level')
subplot(2,1,2)
hist(imag_noise)

for i=1:m
  M(i,2)=realHz(i);
  M(i,3)=imagHz(i);
  M(i,4)=error(i);
end

%write data

str='pr1_prisms61.noise.tip';
fid1=fopen(str,'w+');
fprintf(fid1,'%3.0f',S_num);
fprintf(fid1,'\n')

for j=1:S_num
   fprintf(fid1,'%3.0f',F_num);
   fprintf(fid1,'\n');
   for i=1:F_num
      fprintf(fid1,'%11.5f',M(i+(j-1)*F_num,1));
      fprintf(fid1,'      (');
      fprintf(fid1,' %+11.5E',M(i+(j-1)*F_num,2));
      fprintf(fid1,' ,');
      fprintf(fid1,' %+11.5E',M(i+(j-1)*F_num,3));
      fprintf(fid1,' )');
      fprintf(fid1,'   %8.4f',M(i+(j-1)*F_num,4));
      fprintf(fid1,'\n');
   end
end

fclose(fid1);
