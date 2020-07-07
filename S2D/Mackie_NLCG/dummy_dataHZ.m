clear all;close all;

%read data

S_num=18;
%F_num=5;

inc =1/4

index =0
for iper = -3.3:inc:2.6
   index=index+1
   T(index) = 10^iper   
end
   

T_num = length(T)


real_tip  = 0
imag_tip  = 0
tip_error = 0.01

%write data

str='input.tip';
fid1=fopen(str,'w+');
fprintf(fid1,'%3.0f',S_num);
fprintf(fid1,'\n')

for j=1:S_num
   fprintf(fid1,'%3.0f',T_num);
   %fprintf(fid1,' %5.4f',1);
   fprintf(fid1,' %5.4f');
   fprintf(fid1,'\n');
   for i=1:T_num
      fprintf(fid1,'%10.5f',T(i));
      fprintf(fid1,'      (');
      fprintf(fid1,' %7.2f',real_tip);
      fprintf(fid1,'     ,');
      fprintf(fid1,' %7.2f',imag_tip);
      fprintf(fid1,'    )');
      fprintf(fid1,'   %8.4f',tip_error);
      fprintf(fid1,'\n');
   end
end

fclose(fid1);
