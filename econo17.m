%Wishart matrices, WOE and CWOE analysis
%
clc; clear all; close all;
T=8; N=4;
% %create a random matrix
% rn = randn(r,cl);
% %do local normalization
% for k=1:cl
% for i=1:r
%     A(i,k) = ( rn(i,k)-mean(rn(:,k)) )/(std(rn(:,k)));
% end
% end
% %create WOE
% woe=(A'*A)/cl;
%create a correlation matrix with correlation "corl"
corl=0.6;
E1= corl*ones(N);
E1(logical(eye(size(E1)))) = 1;
E=sqrt(E1);
%create CWOE
A=rand(T,N);
B1=E*A';

% for k=1:N
% for i=1:T
%     M=mean(B1(:,k)) ;
%     std1=std(B1(:,k)); 
%     B1(i,k) = (B1(i,k)-M)/std1;
% end
% end

%cwoe=B1'*B1/T;



% [A1,B1]=correlatedGaussianNoise(A,N);
% A2=A1(1:500,1:1000);
% A3=corrcoef(A2);
% E3=eig(A3);
% D=sign(A3);
% D1=abs(A3);
% D2=D1.^(1.0+e);
% D3=D.*D2;
% 
% E=eig(D3);
% hist(E,20);
    
    