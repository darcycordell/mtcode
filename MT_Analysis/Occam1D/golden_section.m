function [xmin, fmmin, m, F,x3,f4] = golden_section(son,b,W,rgh2,d,rmsdes,ax, bx, cx, tol,nd,w,thick,nl)
% Function which computes golden section algorithm to find minimum
%
% Written by Ersan Turkoglu, circa early 2000s
%
% Updated to work for MT impedances by Darcy Cordell, 2019
b=b';

ndG=nd;

C = (3-sqrt(5))/2;  R = 1-C;
x0 = ax;  x3 = cx;
if (abs(cx-bx) > abs(bx-ax)),
  x1 = bx;  x2 = bx + C*(cx-bx);
else
  x2 = bx;  x1 = bx - C*(bx-ax);
end
A1=x1*rgh2+son;
m1=b/A1;
[F1,~,~]=mt_fwd_occ(10.^m1,nl,w,thick);
f1 = norm(W*F1'-W*d')/sqrt(ndG);
A2=x2*rgh2+son;
m2=b/A2;
[F2,~,~]=mt_fwd_occ(10.^m2,nl,w,thick);
f2 = norm(W*F2'-W*d')/sqrt(ndG);
k = 1;
while abs(x3-x0) > tol*(abs(x1)+abs(x2)),
if f2 < f1,
    x0 = x1;
    x1 = x2;
    x2 = R*x1 + C*x3;   
    f1 = f2;
    A2=x2*rgh2+son;
    m2=b/A2;
    [F2,~,~]=mt_fwd_occ(10.^m2,nl,w,thick);
    f2 = norm(W*F2'-W*d')/sqrt(ndG);
else
    x3 = x2;
    x2 = x1;
    x1 = R*x2 + C*x0;   
    f2 = f1;
    A1=x1*rgh2+son;
    m1=b/A1;
    [F1,~,~]=mt_fwd_occ(10.^m1,nl,w,thick);
    f1 = norm(W*F1'-W*d')/sqrt(ndG);
end
%xx1(k)=x1; ff1(k)=f1; xx2(k)=x2; ff2(k)=f2;  
k = k+1;
end
if f1 < f2,
  xmin = x1;
  fmmin = f1;
  F=F1;
  m=10.^m1;
else
  xmin = x2;
  fmmin = f2;
  F=F2;
  m=10.^m2;
end
if fmmin<rmsdes
    xB=logspace(log10(xmin),log10(cx),20);
    for B=1:length(xB)
         AB=xB(B)*rgh2+son;
         mB=b/AB;
         [FB,~,~]=mt_fwd_occ(10.^mB,nl,w,thick);
         fB(B) = norm(W*FB'-W*d')/sqrt(ndG);
     end
    xx=logspace(log10(xmin),log10(cx),500);
    yy=spline(xB,fB,xx);
    %figure(2);   loglog(xx1,ff1,'*',xx2,ff2,'*');hold on;plot([ax cx],[rms rms]);hold on;plot([ax cx],[rmsdes rmsdes],'--r')
    %hold on;plot(xB,fB,'x','markersize',15); axis([0.0001 1000000 0.1 50])
    %hold on;plot(xx,yy); axis square
    %legend('Golden Search','Golden Search','X^2','Desired X^2','10 point search','spline interpolation')
    xinds=max(find(yy<rmsdes));
    [~, kk]=size(xinds);
    if kk==0
        op=find(min(yy));
        xmin=xx(op);
        A=xmin*rgh2+son;
        m=b/A;
        [F,~,~]=mt_fwd_occ(10.^m,nl,w,thick);
        fmmin = norm(W*F'-W*d')/sqrt(ndG);
        m=10.^m;
    else
        xmin=xx(xinds);
        A=xmin*rgh2+son;
        m=b/A;
        [F,~,~]=mt_fwd_occ(10.^m,nl,w,thick);
        fmmin = norm(W*F'-W*d')/sqrt(ndG);
        m=10.^m;
    end
    %hold on;plot(xmin,fmmin,'or','markersize',15); xlabel('Lagrange multiplayer'); ylabel('RMS Misfit');
    %title('Golden Search')
    %print('-djpeg',[finame,'_goldenS_iter_',num2str(iter-1)]);
end
% figure(3);
 x3=logspace(-5,4,30);
 for ii=1:30
     A3=x3(ii)*rgh2+son;
     m3=b/A3;
     [F3,~,~]=mt_fwd_occ(10.^m3,nl,w,thick);
     f4(ii) = norm(W*F3'-W*d')/sqrt(ndG);
 end
% xminr=floor(min(x3)); xmaxr=ceil(max(x3)); yminr=floor(min(f4)); ymaxr=ceil(max(f4));
% figure(3); semilogx(x3,f4,'.'); hold on; plot(xmin,fmmin,'p','markersize',20); xlabel('Lagrange Multiplier (mu)'); ylabel('RMS Misfit');
% title('Lagrange Multiplier versus RMS Misfit')
% print('-djpeg',[finame,'_muVSrms_iter_',num2str(iter-1)]);
% pause

%-----------------------manual-----------------------
% if fmmin<X2c/2;
% clf
%      loglog(xx1,ff1,'*',xx2,ff2,'*',x3,f4,'o');hold on;%plot(X2c,'o','markersize',15);hold on;plot(X2c/2,'or','markersize',15);ylabel('RMS')
%      [xmin,y]=ginput(1);
%      A=xmin*rgh2+son;
%      m=b/A;
%      F=mteturk(10.^m,choiceF);
%      fmmin = norm(W*F'-W*d')/sqrt(ndG);
% end  