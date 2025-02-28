%histogram of a sigle time window of correlation cefficients,
%NO difference between correlation of general returns and normalized returns
%eigenvalue distribution
clc;clear all; close all;
cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market/');
% cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
 a = csvread('Nikkei_165.csv',1,1); Input1=1;
%a = csvread('SP_194.csv',1,1); Input1=2;
%dates;
s = a;
dim=size(a);
%format long;
totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
    %return1a(i,k) = (return1(i,k)-mean(return1(:,k)))/(std(return1(:,k)));
end
end

I2=0;wind=40;
for t =0:wind:dim(1)-wind
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
end

% plot(a(:,34));
plot(return1(:,34));
t1=[1:length(a)-1];
str1=date_string(t1);
t2=[1:250:length(a)-1];
set(gca,'XTick',t2,'XTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;
set(gca,'XTickLabelRotation',90);

% C = corrcoef(return1);
% D = return3(:,:,34);
% figure(1);im=imagesc(C);caxis([-1 1]);
% B = reshape(D,[],1);
% [x,y]=hist(B,50);
% plot(y,x,'-');

% Y=pdf(B,x);

% E=eig(C);
% hist(E(1:length(E),1),70);
% [x,y]=hist(E(1:length(E)-1,1),100);
% plot(y,x,'o');

