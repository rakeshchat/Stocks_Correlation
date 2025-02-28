%Singular value decomposition
clc; clear all; close all;
cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1;
a = csvread('SP_194.csv',1,1); Input1=2;
dates;
s = a;
dim=size(a);
totalstock=dim(2);
stocksquare=totalstock^2;

for k=1:totalstock
    for i=1:dim(1)-1
        return1(i,k) = log(s(i+1,k)) - log(s(i,k));
    end
end

frame1=1; frame2=7;
wind=1000; overlap=0;
I2=0; step=wind-overlap;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
    E(:,I2)=eig(return3(:,:,I2));
end

E1=E(1:end-20,:);
% plot(E1(1:end-20,:));
dim=size(E1);
m=dim(1);
ns=dim(2);
E2 = zeros(m,ns);
[U,S,V]=svd(E1);
% semilogy(S,'o');

[r1 , r2] = size(S);
r         = min(r1,r2);
sigma     = diag(S);
pvar      = diag(S).^2;
% unfolded eigenvalues
ER = zeros(m,ns);
EG = zeros(m,ns);
EL = zeros(m,ns);
nt = 2;
% data adaptative unfolding
for i=1:m
  for j=1:ns
    sumg = 0;
    for k=1:nt
      xk   = sigma(k)*U(i,k)*V(j,k);
      sumg = sumg + xk;
    end
    EG(i,j) = sumg;
    suml = 0;
    for k=nt+1:r
      xk   = sigma(k)*U(i,k)*V(j,k);
      suml = suml + xk;
    end
    EL(i,j) = suml;
  end
end

ev = reshape(EG,1,m*ns);
% eig density
%dx = 0.01; % bin size
dx = 0.005; % bin size
%dx = 0.001; % bin size
[count,xe] = hist(ev,-2:dx:2);
% figure('Visible','Off');

figure();
for(i=1:ns)
    %loglog(EG(:,i),'o');
    plot(EG(:,i),'o');
    hold on;
end