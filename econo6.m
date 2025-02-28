%distance measurements-SIMILARITY MATRIX
%similarity matrices <|c_ij(T1)-c_ij(T2)|>
%correlation coef between correlation matrices of different windows (dist2)
%sqrt[2(1-dist2)]
%K-means clustering

clc; clear all; close all;
cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market/');
% cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1;
a = csvread('SP_194_new.csv',1,1); Input1=2;
%dates;
s = a;
dim=size(a);
totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
end
end
% avgC = corrcoef(return1);

% for USA total data 8068/8213, for JPN total data 7998
frame1=0; frame2=20;
wind=20; overlap=10;
I2=0; step=wind-overlap;
e=0.6;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
     I2=I2+1;
     return2=corrcoef(return1(t+1:t+wind,:));
     return2(isnan(return2))=0;
     return3(:,:,I2)=return2;
     
     coA1(:,:,I2)=sign(return3(:,:,I2));
     coB1(:,:,I2)=abs(return3(:,:,I2));
     coC1(:,:,I2)=coB1(:,:,I2).^(1.0+e);
     return4(:,:,I2)=(coA1(:,:,I2).*coC1(:,:,I2));
     
%      B=reshape(return2,[],1);
%      C(:,I2)=B;
end

for j=1:I2
   for k=1:I2
     dist1(j,k)=sum(sum(abs(return3(:,:,j) - return3(:,:,k))))/totalstock^2;
     nlin_dist1(j,k)=sum(sum(abs(return4(:,:,j) - return4(:,:,k))))/totalstock^2;
      A=return3(:,:,j);B=return3(:,:,k);
      dist2(j,k)=corr2(A(:),B(:));
%      dist3(j,k)=sqrt(2.0*(1.0-dist2(j,k)));
   end
end

figure(1);im1=imagesc(dist1);
 colorbar; colormap (jet);
figure(2);im1=imagesc(nlin_dist1);
 colorbar; colormap (jet);
  colorbar; colormap (redblue); colormap (flipud (redblue));
 figure(3);im2=imagesc(dist2);
%  colorbar; colormap (redblue); colormap (flipud (redblue));
% figure(4);im3=imagesc(dist3);
%  colorbar; colormap (redblue); colormap (flipud (redblue));


figure(5); plot (dist2(:,6));

%==============================================================
%distance points in eucledian space: to do K-means clustering
%==============================================================
  y=cmdscale(nlin_dist1,3); %3 for 3D plot
  
  [idx,C]=kmeans(y,4); %for number of K-means cluster
  
  figure(6);
  plot3(y(idx==1,1),y(idx==1,2),y(idx==1,3),'r.','MarkerSize',12);
  hold on;
  plot3(y(idx==2,1),y(idx==2,2),y(idx==2,3),'b.','MarkerSize',12);
  hold on;
  plot3(y(idx==3,1),y(idx==3,2),y(idx==3,3),'g.','MarkerSize',12);
  hold on;
  plot3(y(idx==4,1),y(idx==4,2),y(idx==4,3),'m.','MarkerSize',12);
  hold on;

  
%   plot3(y(idx==5,1),y(idx==5,2),y(idx==5,3),'c.','MarkerSize',12);
%   hold on;


%   n=length(dist1);
%   for k=1:n
%   text(y(k,1),y(k,2),y(k,3),num2str(k)); 
%   end