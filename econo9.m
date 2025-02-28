%distance between correlation matrices |<c_ij(T1)-c_ij(T2)>|
%correlation coef between correlation matrices of different windows
clc; clear all; close all;
%cd ('/home/rakesh/Dropbox/econophys/states_financial_market');
% a = csvread('Nikkei_165.csv',1,1); Input1=1; 
a = csvread('SP_194.csv',1,1); Input1=2;
dates;
s = a;
dim=size(a);
totalstock=dim(2);
for k=1:totalstock
for i=1:dim(1)-1
    return1(i,k) = log(s(i+1,k)) - log(s(i,k));
end
end

frame1=1; frame2=30;
wind=20; overlap=10;
I2=0; step=wind-overlap;

start1=frame1*wind; stop1=frame2*wind;

v=VideoWriter('Window_40.avi');
open(v);
for t=start1:step:stop1
    I2=I2+1;
    return2=corrcoef(return1(t+1:t+wind,:));
    return2(isnan(return2))=0;
    im1=imagesc(return2);title(['S-Frame=',num2str(I2),' date:',date_string(t+wind)]);
    caxis([-1, 1]);colorbar;colormap(jet);
%     colorbar; colormap (redblue); colormap (flipud (redblue));
%     cd ('figures/corrmat_frames_S500/'); saveas(im1,['S-Frame_',num2str(I2)],'png'); cd ('../../');
    F=getframe(gcf);
    writeVideo(v,F);
    return3(:,:,I2)=return2;
end

close(v);
close(gcf);