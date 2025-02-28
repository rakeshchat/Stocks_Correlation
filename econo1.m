%return distribution of a single stock and normalized return dist
%s=load('/Users/rakeshchatterjee/1.codes/econophys/states_financial_market/Japan_375_stocks.csv');
clc; clear all; close all;
cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market/');
% a = csvread('Nikkei_165.csv',1,1);
a = csvread('SP_194.csv',1,1);
%('filename',1,1)-> last two 1,1 is to define the origin of file to read 
s = a(41:80,:);
for i=1:length(s)-1
    return1(i) = log(s(i+1)) - log(s(i));
%     return2(i) = (s(i+1) - s(i))/s(i);
%     return3(i) = (return2(i)-mean(return2))/(std(return2));
end
[x,y]=hist(return1,50);
% [x1,y1]=hist(return2,50);
% [x2,y2]=hist(return3,50);
%f1 = fit(y.',x.','gauss1');
%aa=15:101;
% plot(f1,y.',x.');
figure(1);plot(y,x,'*-');title('General Return distribution');
% figure(2);plot(y2,x2,'o-');title('Normalized Return distribution');
%semilogy(autocorr(return1.^2,100),'.-');title('Autocorrelation: Return^2:semilog');
%plot(autocorr(return1.^2,100),'.-');title('Autocorrelation: Return^2');
%auto=autocorr(return3,100);
%f2 = fit(aa.',auto(1,15:101).','power1');
%plot(f2,aa,auto(1,15:101));