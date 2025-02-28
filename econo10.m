%randomly select matrix from a window of critical period and their
%maximun eigenvalue distribution
clc; clear all; close all;
cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1);
a = csvread('SP_194.csv',1,1);
s = a;
dim=size(a);
%format long;
totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
end
end
I2=0;wind=40;
for t=0:wind:length(return1)-wind
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
end
lehman = return3(:,:,34);
maxmatrix=10000;
window=40; block=35;
%maxmatrix should be less than (totalstock_C_block)
%window=time window size, block=number of stocks randomly chosen
I1=0;
for k=1:maxmatrix
    I1=I1+1;
    raa=randperm(totalstock,block);
    ra12=lehman(1:window,raa);
    ra22=corrcoef(ra12);
    ra32(:,:,I1)=ra22;
    mEg(I1)=max(eig(ra32(:,:,I1)));
end
hist(mEg,60);