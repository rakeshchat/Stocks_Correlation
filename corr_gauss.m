% %%%%%%%% WOE Wishart orthogonal ensemble %%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all;
% clc;
% N=1024;T=N*0.5;%% kappa = T/N=0.5 Vinayak paper
% en=200;
% C_woe=zeros(N,N,en);
% for k=1:en
%     M=rand(N,T);    size_M=size(M);
%     M_norm=(M-mean(M,2))./(std(M,0,2)); 
%     C_woe(:,:,k)=(M_norm*M_norm')/T;
%     k
% end
% eps=0.001;TT=512;
% clear nbin1 nbin2 a11 b11 a22 b22
%  for i=1:en
%      A11=C_woe(:,:,i);
%     r1=A11;r1(isnan(r1))=0;
%     r2=sign(r1);r22=abs(r1);mean(mean(r1));
%     r3=r22.^(1+eps);r4=r2.*r3;
%     Eig2=eig(r4); cut=N-TT+1;
%     Eig_ES=Eig2(1:cut);
%     %%%%%% Normal Spectra  %%%%%%%
%     nbin1=linspace(0,6.2,100);
%     [a1,a2]=hist(Eig2,nbin1); a11(i,:)=a1; a22(i,:)=a2;
%     %%%%%% Emerging spectra  %%%%%%%
%     nbin2=linspace(-.001,.006,100);
%     [b1,b2]=hist(Eig_ES,nbin2);b11(i,:)=b1;b22(i,:)=b2;
%  end
% a111=mean(a11,1);a222=mean(a22,1);
% figure();plot(a222,a111);%xlim([.1,10])
% xlabel('\lambda');ylabel('$\bar{\rho}(\lambda)$','Interpreter','Latex');
% set(gca,'LineWidth',1.5);figurepalette('show');
% set(gca, 'FontSize',14,'FontWeight','bold');
% title('N=194, T=3000, C=0.3, eps=0.001')
% b111=mean(b11,1);b222=mean(b22,1);
% figure();plot(b222,b111);%xlim([-.1,.1])
% xlabel('\lambda');ylabel('$\bar{\rho}(\lambda)$','Interpreter','Latex');
% set(gca,'LineWidth',1.5);figurepalette('show');
% set(gca, 'FontSize',14,'FontWeight','bold');
% title('N=194, T=3000, C=0.3, eps=0.001')
%%   CWOE sectorwise   %%%%%%%%%%%%%%%%2.0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100;r1=rand(4);sum_r1=sum(r1);i1=r1/sum_r1;
n1=30; n2=30; n3=20; n4=N-n1-n2-n3;
T=40;
C1=0.9;E = C1.*ones(N); E(logical(eye(size(E)))) = 1;

C1=0.3;E1 = C1.*ones(n1); E1(logical(eye(size(E1)))) = 1;
E(1:n1,1:n1)= E1;

C2=0.6;E2 = C2.*ones(n2);E2(logical(eye(size(E2)))) = 1;
E(n1+1:n1+n2,n1+1:n1+n2)= E2;

C3=0.9;E3 = C3.*ones(n3); E3(logical(eye(size(E3)))) = 1;
E(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3)= E3;

C4=0.9;E4 = C4.*ones(n4); E4(logical(eye(size(E4)))) = 1;
E(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4)= E4;
imagesc(E)
%% No power mapping eps=0
% clear nbin a11 a22
% for i=1:dim_woe(3)
%     A11=C_woe(:,:,i);    
%     r1=corrcoef(A11);r1(isnan(r1))=0;
%     Eig2=eig(r1);
%     nbin=linspace(0,5,100);
%     [a1,a2]=hist(Eig2,nbin);   
%     a11(i,:)=a1;a22(i,:)=a2;
% end
% a111=mean(a11);a222=mean(a22);
% figure();plot(a222,a111,'o-')
% xlabel('\lambda');ylabel('$\bar{\rho(\lambda)}$','Interpreter','Latex');
% set(gca,'LineWidth',1.5);figurepalette('show');
% set(gca, 'FontSize',14,'FontWeight','bold');
%% CWOE Correlated Wishart orthogonal ensemble   %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%   CWOE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;
N=1000;T=N*5;C1=0.9;
E = C1.*ones(N);E(logical(eye(size(E)))) = 1;
Rm=rand(N,T);B=E^(1/2)*Rm;dimB=size(B);
B_norm=(B-mean(B,2))./(std(B,0,2));
C=B_norm*B_norm'/T;
eps=0.001;TT=50;i1=0;
clc;clear nbin1 nbin2 a11 b11 a22 b22
mean(mean(C))
for i=1:TT:T
    i1=i1+1;
    A11=B_norm(:,i:i+TT-1);
    r1=corrcoef(A11');r1(isnan(r1))=0;
    r2=sign(r1);r22=abs(r1);mean(mean(r1));
    r3=r22.^(1+eps);r4=r2.*r3;
    Eig2=eig(r4); cut=N-TT+1;
    Eig_ES=Eig2(1:cut);
    fprintf('i1=%g, min=%g, max=%g, max=%g\n',i1, min(Eig_ES),max(Eig_ES),max(Eig2));
    %%%%%% Normal Spectra  %%%%%%%
    nbin1=linspace(-1,1000,300);
    [a1,a2]=hist(Eig2,nbin1); a11(i1,:)=a1; a22(i1,:)=a2;
    %%%%%% Emerging spectra  %%%%%%%
    nbin2=linspace(-0.00001,.0001,100000);
    [b1,b2]=hist(Eig_ES,nbin2);b11(i1,:)=b1;b22(i1,:)=b2;
end
a111=mean(a11,1);a222=mean(a22,1);
figure();bar(a222,a111);%xlim([.1,10])
xlabel('\lambda');ylabel('$\bar{\rho}(\lambda)$','Interpreter','Latex');
set(gca,'LineWidth',1.5);figurepalette('show');
set(gca, 'FontSize',14,'FontWeight','bold');
title(['N=', num2str(N), ', T=', num2str(TT),', C=', num2str(C1), ', eps=', num2str(eps)])

b111=mean(b11,1);b222=mean(b22,1);
figure();bar(b222,b111);%xlim([-.1,.1])
xlabel('\lambda');ylabel('$\bar{\rho}(\lambda)$','Interpreter','Latex');
set(gca,'LineWidth',1.5);figurepalette('show');
set(gca, 'FontSize',14,'FontWeight','bold');
title(['N=', num2str(N), ', T=', num2str(TT),', C=', num2str(C1), ', eps=', num2str(eps)])

% N=n1;
% C1=0.4;E = C1.*ones(N);
% E(logical(eye(size(E)))) = 1;
% imagesc(E)
% 
% M=rand(N,T);
% B=E^(1/2)*M;
% dimB=size(B);
% for i=1:dimB(1)
%     for j=1:dimB(2)
%         M1=mean(B(i,:));
%         Std1=std(B(i,:));
%         B_norm(i,j)=(B(i,j)-M1)/(Std1);
%     end
% end
% C=B_norm*B_norm'/T;
% imagesc(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% N=1024;C1=0.8;
% A = C1.*ones(N);
% A(logical(eye(size(A)))) = 1;
% [A1,B1] =(correlatedGaussianNoise(A,N));
% A11=A1(1:20,:);
% Eig1=eig(corrcoef(A11));
% figure(1);hist(Eig1,20);
% % plot(A1');
% eps=0.001;
% r1=corrcoef(A11);%corrcoef(M3d(:,:,t));r1(isnan(r1))=0;
% r2=sign(r1);r22=abs(r1);
% r3=r22.^(1+eps);
% r4=r2.*r3;nbin=linspace(-0.02,0.1,100);
% R2=eig(r4);
% figure(2);hist(R2(1:(end-50)),10)
% R11=R2(R2<0);neg_R(t)=length(R11);
% if Input1==1;cut=dim_r(2)-wind+1;end
% if Input1==2;cut=dim_r(2)-wind+1;end
% h1=figure(1);hist(R2(1:cut),nbin);
%title(['time=',num2str(t),'date',date_string(bb)]);set(gca,'FaceColor',[1 0.6 0.2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

