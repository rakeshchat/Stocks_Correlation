%all correlation matrices with overlapping windows
clc; clear all; close all;
cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1;
a = csvread('SP_194.csv',1,1); Input1=2;
dates;
s = a;
dim=size(a);
totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
end
end

frame1=0; frame2=50;
wind=20; overlap=10;
I2=0; step=wind-overlap;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
     I2=I2+1;
     return2=corrcoef(return1(t+1:t+wind,:));
     return2(isnan(return2))=0;
     return3(:,:,I2)=return2;
     B=reshape(return2,[],1);
     C(:,I2)=B;
end

% for j=1:I2
%    for k=1:I2
%      dist(j,k)=sum(sum(abs(return3(:,:,j) - return3(:,:,k))))/totalstock^2;
%    end
% end
% im=imagesc(dist);

for i=1:I2;
    ar=linspace(-1,1,90);
    r1(i,:)=hist(C(:,i),ar);
end
dim1=size(r1);
[X,Y] = meshgrid(linspace(-1,1,dim1(2)),linspace(-1,dim(1),dim1(1)));
figure(2);surf(X,Y,r1); shading interp;