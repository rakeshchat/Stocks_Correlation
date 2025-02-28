%all correlation matrices
%average correlation of all the returns
%similarity matrix from all such correlation matrix 
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

avgcorr=corrcoef(return1);

frame1=1; frame2=200;
wind=20; overlap=10;
I2=0; step=wind-overlap;

start1=frame1*wind; stop1=frame2*wind;
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    return3(:,:,I2)=return2;
end

% for j=1:I2
%    for k=1:I2
%      dist(j,k)=sum(sum(abs(return3(:,:,j) - return3(:,:,k))))/totalstock^2;
%    end
% end

%======================
%overall average correlation
%======================
% im=imagesc(avgcorr);

%======================
%correlation matrix of a particular frame
%======================
%im=imagesc(return3(:,:,68)); caxis([-1,1]);


% if Input1==2;
%     t1=[22 41 57 90 110 142 160 173 178 194];
%     caxis([-1,1]);
%     set(gca,'XTick',t1,'XTickLabel',{'CD' 'CS' 'EG' 'FN' 'HC' 'ID' 'IT' 'MT' 'TC' 'UT'})
%     set(gca,'YTick',t1,'YTickLabel',{'CD' 'CS' 'EG' 'FN' 'HC' 'ID' 'IT' 'MT' 'TC' 'UT'})
%     set(gca,'Fontsize',15);set(gca,'FontWeight','b');
%     set(gca,'XTickLabelRotation',90)
% end
% 
% if Input1==1;
%     t1=[34 50 55 107 113 149 165];
%     caxis([-1,1]);
%     set(gca,'XTick',t1,'XTickLabel',{'CAP' 'CON' 'FIN' 'MAT' 'PHR' 'TECH' 'UTL'})
%     set(gca,'YTick',t1,'YTickLabel',{'CAP' 'CON' 'FIN' 'MAT' 'PHR' 'TECH' 'UTL'})
%     set(gca,'Fontsize',15);set(gca,'FontWeight','b');
%     set(gca,'XTickLabelRotation',90)
% end

%======================
%correlation with time
%======================
im=imagesc(return3);
colorbar; colormap (redblue); colormap (flipud (redblue));
t1=[0:wind:dim(1)-wind];
str1=date_string(t1+wind);
t2=[1:30:length(dist)];caxis([-1,1]);
set(gca,'XTick',t2,'XTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;
set(gca,'XTickLabelRotation',90);
set(gca,'YTick',t2,'YTickLabel',str1(t2)');
set(gca,'Fontsize',10);set(gca,'FontWeight','b');grid on;


colorbar; colormap (jet);
saveas(im,['1_fig'],'png');