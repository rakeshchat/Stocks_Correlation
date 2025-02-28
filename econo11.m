%non-linear transformation of correlation matrices
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

frame1=20; frame2=50;
%select window, overlap sizes and e
wind=20; overlap=19;
I2=0; step=wind-overlap; gap=wind-1;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
end

%selecting correlation matrix of a particular window number
corln1=return3(:,:,34);
e=0.01;

A1=sign(corln1);
B1=abs(corln1);
C1=B1.^(1.0+e);
D1=(A1.*C1);
E1=eig(D1);
figure(1);hist(E1(1:length(E1)-40,1),20);