%non-linear transformation of correlation matrices
clc; clear all; close all;
cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
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

frame1=20; frame2=50;

%select window, overlap sizes and e
wind=20; overlap=19; e=0.01;

h=0; step=wind-overlap; gap=wind-1;
start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    h=h+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,h)=return2;
    
A(:,:,h)=sign(return3(:,:,h));
B(:,:,h)=abs(return3(:,:,h));
C(:,:,h)=B(:,:,h).^(1.0+e);
D(:,:,h)=(A(:,:,h).*C(:,:,h));

E(:,h)=eig(D(:,:,h));
R(:,h)=(E(:,h)<0);
R1(h)=sum(R(:,h));
M1(h)=mean(mean(return3(:,:,h)));
end



N=(R1-mean(R1))/std(R1);
M=(M1-mean(M1))/std(M1);

plot(M,'b.-');
hold on;
plot(N,'r.-');legend('mean','# of -ve Eig val');
