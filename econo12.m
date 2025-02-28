%correlation matrix analysis with RANDOM TIME SERIES
%non-linear transform analysis
%in a period, choosing randomly selected matrix and max e.v. dist
clc; clear all; close all;
%generating random time series of returns
length1=10000; column1=2500;
r = randn(length1,column1);
ret1 = r(41:80,:);
%==========================================
%non-linear transformation
%==========================================
% R = corrcoef(ret1);
% R1=corrcoef(r);
% A=sign(R);
% B=abs(R);
% e=0.5;
% C=B.^(1.0+e);
% D=(A.*C);
% E=eig(D);
%hist(D(1:length(D),1),100);
% E1=eig(R1);
% hist(E1(1:length(E1),1),60);
%[x,y]=hist(E(1:length(E),1),50);
%plot(y,x,'o');
%======================================
%max e.v. dist. at the critical period
%======================================
maxmatrix=10000;
window=40; block=35;
I1=0;
for k=1:maxmatrix
    I1=I1+1;
    ra=randperm(column1,block);
    ra1=ret1(1:window,ra);
    ra2=corrcoef(ra1);
    ra3(:,:,I1)=ra2;
    mE(I1)=max(eig(ra3(:,:,I1)));
end
figure(1);hist(mE,50);