%surfaceplot of time evolution of correlation coefficients of return stocks
clc; clear all; close all;
%cd ('/Users/rakeshchatterjee/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1;
a = csvread('SP_194_new.csv',1,1); Input1=2;
%dates;
s = a;
dim=size(a);
totalstock=dim(2);
stocksquare=totalstock^2;
for k=1:totalstock
    for i=1:dim(1)-1
        return1(i,k) = log(s(i+1,k)) - log(s(i,k));
    end
end

frame1=1; frame2=200;
wind=40; overlap=0;
I2=0; step=wind-overlap;

% end1=floor(dim(1)/40);
start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    ut1=triu(return2,1);
    ut1(ut1==0)=[];
    ar=linspace(-1,1,90);
    h1(I2,:)=hist(ut1,ar);
    pdf1(I2,:)=hist(ut1,ar)/length(ut1);
end

%dimm2=size(pdf1);
%[X,Y] = meshgrid(linspace(-1,1,dimm2(2)),linspace(1,dimm2(1),dimm2(1)));

h1=figure(1);
surf(pdf1); shading interp; colormap('jet'); colorbar; hold on;



%surf(X,Y,pdf1); shading interp; colormap('jet'); colorbar; hold on; %view(0,90);
%xlim([-1 1]);ylim([1 dimm2(1)]);
% h1=surf(X,Y,pdf1);set(h1,'facecolor','none');hold on;

t1=[start1:step:stop1];
%str1=date_string(t1+wind);
t2=[1:50:length(Y)];
%set(gca,'YTick',t2,'YTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;
% set(gca,'XTickLabelRotation',90);


xlabel('C_{ij}','Fontsize',16,'FontWeight','b');
set(gca,'XTick',[-1 0 1],'XTickLabel',{'-1' '0' '1'})
set(gca,'Fontsize',16);set(gca,'FontWeight','b');