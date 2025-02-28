%mean of correlation matrix over time and the Maximun of eigenvalue of each
%maximun eigenvalues captures the average correlation (when having few
%anti-correlations)
%measurement of skewness and kurtosis
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

wind=20; overlap=10;
I2=0; step=wind-overlap;
frame1=1; frame2=402;

start1=frame1*wind; stop1=frame2*wind;
for t =start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
    E(:,I2)=eig(return3(:,:,I2));
    Emax(:,I2)=max(eig(return3(:,:,I2)));
    M(:,I2)=mean(mean(return3(:,:,I2)));
    B=reshape(return2,[],1);
    C(:,I2)=B;
end

V=var(C);
S=skewness(C);
K=kurtosis(C);
    
%plot(1:139,M/max(M),'-*',1:139,Emax/max(Emax),'-o');
% figure(1);im1=plot(S,'-*');
% figure(2);im2=plot(K,'-o');

% title(['SP500(mean,var,skew,kurt)--Window=',num2str(wind),', Overlap=',num2str(overlap)]);


% subplot(411);plot(Emax,'g.-');legend('Max EV');grid on
subplot(411);plot(M,'r.-');legend('Mean');grid on;set(gca,'Xtick',[])
subplot(412);plot(V,'b.-');legend('Variance');grid on;set(gca,'Xtick',[])
subplot(413);plot(S,'c.-');legend('Skewness');grid on;set(gca,'Xtick',[])
subplot(414);plot(K,'k.-');legend('Kurtosis');grid on


t1=[start1:step:stop1];
str1=date_string(t1+wind);
t2=[1:40:I2];
set(gca,'XTick',t2,'XTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;
set(gca,'XTickLabelRotation',90);

