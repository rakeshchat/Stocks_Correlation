%%% Intra cluster distance
%% hirdesh kumar pharasi
function [sum1, avg1] = intra_cluster_hkp(a,nc) 
%  a=[7 0;6 0;5 0;4 0; 3 0; 0 3;0 4;0 5;0 6;0 7;-3 0;-4 0; -5 0;-6 0;-7 0; 0 -7;0 -6;0 -5;0 -4; 0 -3];
%  nc=5;
% fprintf('nc=%g ',nc)
[idx,cen]=kmeans(a,nc);
for i=1:nc
    idx1{i}=find(idx==i);    
end
for i=1:nc
   UC{i}=a(idx1{i},:);
end
sum1=0;avg1=0;
for i=1:nc
    [r1]=size(UC{i});
    Ucc=UC{i};pdist1=0;
    for j=1:r1
        pdist1(j)=pdist([Ucc(j,:);cen(i,:)]);
    end  
    sum1(i)=sum(pdist1);%% sum_intra
    avg1(i)=mean(pdist1);%%avg_sum_intra
end
% %%%%%%%%%%%%%5
% y=a;rng(2);[idx,C] = kmeans(y,nc);
% Cl = {'[0 0 0]','[1 0 1]','[0 1 0]','[0 0 0.5]','[1 0 0]'...
%      '[0 0 1]','[0 0 .6]','[0.5 0 0.5]','[0.5 0 0]','[0.5 0 0.1]'};
% n=length(a);k = 1:n;cen
% h1=figure()
% for i=1:nc
%     text(y(k(idx==i),1),y(k(idx==i),2),num2str(find(idx==i)),'FontSize',12,'color',Cl{i});
% end
% xlim([min(y(:,1))-0.02, max(y(:,1))+0.02]); ylim([min(y(:,2))-0.02, max(y(:,2))+0.02]);
% % if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% % if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% set(gca, 'FontSize',12,'FontWeight','bold');
% xlabel('x-coordinates');ylabel('y-coordinates');figurepalette('show')
% 
% 
