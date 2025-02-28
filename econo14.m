%emerging spectra: negative eigenvalues follows mean, for overlaping window
clc; clear all; close all;
cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1; 
a = csvread('SP_194.csv',1,1); Input1=2;
dates;
s = a;
dim=size(s);

totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
end
end

%for scanning a period of time, set frame1 and frame2
frame1=20; frame2=50;

%select window, overlap sizes and e
wind=20; overlap=19; e=0.01;
h=0; step=wind-overlap; gap=wind-1;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
%     t2(t)=t; disp(date_string(t)); fprintf('%d',t);
    h=h+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,h)=return2;
    
    A(:,:,h)=sign(return3(:,:,h));
    B(:,:,h)=abs(return3(:,:,h));
    C(:,:,h)=B(:,:,h).^(1.0+e);
    D(:,:,h)=(A(:,:,h).*C(:,:,h));

    E(:,h)=eig(D(:,:,h));
    Mxeg1(h)=max(E(1:end-gap,h));
    Mneg1(h)=min(E(1:end-gap,h));
    V1(h)=var(E(1:end-gap,h));
    S1(h)=skewness(E(1:end-gap,h));
    K1(h)=kurtosis(E(1:end-gap,h));
    R(:,h)=(E(:,h)<0);
    N1(h)=sum(R(:,h));
    M1(h)=mean(mean(return3(:,:,h)));
end

%========================================================
%Plot Subplot
%========================================================
% subplot(411);plot(M1,'r.-');legend('Mean');grid on
% title(['SP500(mean,var,skew,kurt)--Window=',num2str(wind),', Overlap=',num2str(overlap)]);
% subplot(412);plot(V1,'b.-');legend('Variance');grid on
% subplot(413);plot(S1,'y.-');legend('Skewness');grid on
% subplot(414);plot(K1,'k.-');legend('Kurtosis');grid on
%========================================================
subplot(411);plot(M1,'r.-');legend('Mean');grid on
% title(['SP500(mean,total-min-max)--Window=',num2str(wind),', Overlap=',num2str(overlap)]);
subplot(412);plot(Mneg1,'b.-');legend('Min -ve EV');grid on
subplot(413);plot(Mxeg1,'k.-');legend('Max -ve EV');grid on
subplot(414);plot(N1,'g.-');legend('Total # of -ve EV');grid on
%========================================================

t1=[start1:step:stop1];
str1=date_string(t1+wind);
t2=[1:20:length(M1)];
set(gca,'XTick',t2,'XTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;
set(gca,'XTickLabelRotation',90);

%correlation between two series:
% corr2(M1,N1);

%========================================================
%Plot with normalized
%========================================================
% N=(N1-mean(N1))/std(N1);
% M=(M1-mean(M1))/std(M1);
% V=(V1-mean(V1))/std(V1);
% S=(S1-mean(S1))/std(S1);
% Mxeg=(Mxeg1-mean(Mxeg1))/std(Mxeg1);
% Mneg=(Mneg1-mean(Mneg1))/std(Mneg1);

% plot(N,'b.-');
% hold on;
% plot(M,'r.-');
% hold on;
% plot(V,'g.-');
% hold on;
% plot(S,'k.-');
% hold on;
% plot(MegN,'y.-');
% legend('# of -ve EV','Mean','variance','skewness','Max of EV');

%=========================================



% function [F,c_v] = granger_cause(x,y,alpha,max_lag)
% alpha -- the significance level specified by the user
% max_lag -- the maximum number of lags to be considered
% if F > c_v we reject the null hypothesis 
% that y does not Granger Cause x

%for i=0.01:0.01:0.1
% i
% [F1,cv1]=granger_cause(N,M,0.04,2)
%end

% tau = 1; % time step 
% E   = 4; % dimension of reconstruction
% LMN = 5; % number of neigborhoods for L and M methods
% L=100:100:1000; % vector of computing points
% [ccm1]=CCM( M , N , tau , E , LMN , L )