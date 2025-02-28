%Adjacency Matrix from correlation matrix, degree distribution, 
%forming network with a threshold value
%forming network with only highly correlated and highly anticorrelated
%stocks
clc; clear all; close all;
cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market');
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

frame1=20; frame2=50;


wind=20; overlap=19;
I2=0; step=wind-overlap;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
end

ret1 = return3(:,:,94);
R = corrcoef(ret1);
% R1 = R(100:105,100:105);

th=0.15;
A=(R < -0.6 | R > 0.8);
for i=1:totalstock
x(i)=sin(2*pi*i/totalstock);
y(i)=cos(2*pi*i/totalstock);
xy(i,:)=[x(i) y(i)];
end
gplot(A,xy);
S=sum(A);
% [x,y]=hist(S,60);
% plot(y,x,'o-');