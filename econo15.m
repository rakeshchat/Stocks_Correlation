%non-linear transformation and -ve eigenvalue analysis
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

avgcorr=corrcoef(return1);
Etotal=eig(avgcorr);

%for scanning a period of time, set frame1 and frame2
frame1=1; frame2=402;

%select window, overlap sizes and e
wind=20; overlap=0; e=0.002;
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
end

%========================================================
%Eigenvalue spectrum (for a particular window(w))
%========================================================

%without distortion:
%dim(2)-(wind-1) zero egvs.,(wind-1) non-zero egvs.
%normal period  ->MP dist. max egv
%critical period->MP dist., max egv shifts to high

%with distortion:
%from zero egvs, emerging spectra(ES) comes out.
%normal period  ->ES looks like distorted semicircle
%critical period->ES looks like lorentz curve

w=36;
Ev=E(:,w);
hist(Ev(1:end),50);
% hist(Ev(1:end-wind+1),50);
% Evn = Ev( Ev<0 ); %only -ve evs of ES
% hist(Evn,15);
% hist(Etotal(1:end),700);