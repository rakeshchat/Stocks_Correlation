clc;close all;clear all;cd ('/home/hirdesh/RMT_finance/MATLAB_prog/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
country=1;  wind=20;    shift=20;
a_country_data;save_plot=322;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:t11
    aa=shift*(t-1)+1;    bb=shift*(t-1)+wind;
    M(:,:)=return_all(aa:bb,:);
    M3d(:,:,t)=return_all(aa:bb,:);
    r1(:,:,t)=corrcoef(M);  
    r1(isnan(r1))=0;
    r2=corrcoef(M);
    r2(isnan(r2))=0;    
    E1(:,t)=eig(r2);    
    fprintf('%g\t%g\t%g\t',t,aa,bb);disp(date_string(bb))
end
% for eps=0.04:0.3:0.94
eps=0.65;
r2=sign(r1);r22=abs(r1);     r3=r22.^(1+eps);     r4=r2.*r3;
r1=r4;
% figure(1);imagesc(r1(:,:,19))
%%%%%%%%%%%%%%%%%%%%% Distance matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1=0;I3=0;dim3=size(r1);
for  i=1:dim3(3)
    for j=1:dim3(3)
        I3=dim3(1)*dim3(2);
         S1(i,j)=sum(sum(abs(r1(:,:,i)-r1(:,:,j))))/I3;%%
%        S1(i,j)=sqrt(sum(sum((r1(:,:,i)-r1(:,:,j)).^2)))/I3; fprintf('l2 distance')
        A=r1(:,:,i); B=r1(:,:,j);
        S2(i,j)=corr2(A(:),B(:));
        S3(i,j)=sqrt(2*(1-S2(i,j)));
    end
end
% h11=figure(11);imagesc(S1);colorbar;colormap1; save_input=61;saveplot;
% %% mds 2D     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y=cmdscale(S1,2);n=size(y,1);h5=figure(5);plot(y(1:n,1),y(1:n,2),'.');
% for k=1:n; text(y(k,1),y(k,2),num2str(k),'FontSize',10); end
% if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% set(gca, 'FontSize',10,'FontWeight','bold');
% xlabel('x coordinates');ylabel('y coordinates');
% xlim([-.3, .6]);ylim([-.1, .15]);
% saveas(h5,['../long_Paper/corr_coef_powermap/mds_S1=',num2str(t), '_100eps=', num2str(eps*100)],'png');
% end
% colormap(jet);title('C_{i,j}(t1)-C_{i,j}(t2)');
% % colorbar; colormap(redblue); colormap(flipud(redblue));
% if country==1 %% japan
%     t1=[1,16,30,43,57,69,82,94,106,119,131,144,156,168,180,193,205,217,230,242,254,267,279,291,303,316,328,340,352,365,377,389];
%     str1=[1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016];
% end
% if country==2 %% usa
%     t1=[1,14,27,39,52,65,77,90,103,115,128,140,153,166,178,191,204,216,229,241,254,266,279,292,305,317,329,342,354,367,380,392];
%     str1=[1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016];
% end
% % t1=[1:300:t11];
% % % str1=[(date_string(shift*(t1-1)+1))];
%  set(gca,'XTick',t1,'XTickLabel',str1(1:length(t1))');
%  set(gca,'XTickLabelRotation',90)
% % set(gca,'YTick',t1,'YTickLabel',str1(1:length(t1)))
% set(gca,'Fontsize',12);set(gca,'FontWeight','b');
%
%
% h22=figure(2);imagesc(S2);title('corr(A,B)');colorbar; colormap(redblue); colormap(flipud(redblue));
% h33=figure(3);imagesc(S3);title('d_{1,2}=sqrt{2(1-C(i,j)}');colorbar; colormap(redblue); colormap(flipud(redblue));
% 
% if Input==333;
%     saveas(h11,'plots/Sim_avg_dist/similarity_Matrix.jpg');
%     saveas(h22,'plots/Sim_avg_dist/correlation_Matrix.jpg');
%     saveas(h33,'plots/Sim_avg_dist/distance_Matrix.jpg');
% end
% %%%%%%%%%%%%%%%%%%%%%%%% MDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % d2=sqrt(2*(1-S2));
% y=cmdscale(S1,2);
% n=size(y,1);
% p1=figure(111);%p2=plot(y(1:n,1),y(1:n,2),'o');
% p2=plot3(y(1:n,1),y(1:n,2),y(1:n,3),'o');
% h = rotate3d;grid on;box on
% h.Enable = 'on';
% title(['Frame ',num2str(t)]);
% 
% for k=1:n
%     text(y(k,1)+0.009,y(k,2),y(k,3),num2str(k));
% end
% %%%%%%%%%%%%%%%%%%%% Kmean clustering 3-D %%%%%%%%%%%%%%5
% y=cmdscale(S1,3);X=y;nc=6; %% number of cluster
% [idx,C] = kmeans(y,nc);hold on;
% Cl = {'[0 0 0]','[1 0 1]','[1 1 0]','[0 0 0.5]','[1 0 0]'...
%     '[0 1 0]','[0 0 1]','[0.5 0 0.5]','[0.5 0 0]'};
% for i=1:nc
%     plot3(X(idx==i,1),X(idx==i,2),X(idx==i,3),'color',Cl{i}, ...
%         'marker','.','MarkerSize',4, 'LineStyle','none');
%     hold on;
% end
% grid on;  
% xlim([min(y(:,1)), max(y(:,1))]); ylim([min(y(:,2)), max(y(:,2))]);zlim([min(y(:,3)), max(y(:,3))]);
% set(gca, 'FontSize',10,'FontWeight','bold');
% xlabel('x coordinates');ylabel('y coordinates');zlabel('z coordinates');
% view([1,1,1]);hold on;n=length(S1);k = 1:n;
% for i=1:nc
%     text(X(k(idx==i),1),X(k(idx==i),2),X(k(idx==i),3),num2str(find(idx==i)),'FontSize',10,'color',Cl{i});
%     hold on;
% end
% h = rotate3d;grid on;box on;h.Enable='on';I=0;
% v = VideoWriter('peaks.avi');% Prepare the new file.
% open(v);
% for j=1:1:90
%     I=I+1;view(30,j);%pause(0.1);
%     Frms = getframe(gcf); writeVideo(v,Frms);
% end
% close(gcf);close(v);
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Kmean clustering 2-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=cmdscale(S1,2);
if country==1;nc=4;end
if country==2;nc=5;end %% number of cluster
%% Note: rng(SD) seeds the random number generator using the non-negative
%% integer SD so that RAND, RANDI, and RANDN produce a predictable
%% sequence of numbers.
rng(2);[idx,C] = kmeans(y,nc);
Cl = {'[0 0 0]','[1 0 1]','[0 1 0]','[0 0 0.5]','[1 0 0]'...
     '[0 0 1]'};%,'[0 0 1]','[0.5 0 0.5]','[0.5 0 0]'};
n=length(S1);k = 1:n;
h1=figure()
for i=1:nc
    text(y(k(idx==i),1),y(k(idx==i),2),num2str(find(idx==i)),'FontSize',12,'color',Cl{i});
end
xlim([min(y(:,1))-0.02, max(y(:,1))+0.02]); ylim([min(y(:,2))-0.02, max(y(:,2))+0.02]);
if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
set(gca, 'FontSize',12,'FontWeight','bold');
xlabel('x-coordinates');ylabel('y-coordinates');figurepalette('show')
%% sequence of numbers.
% rng(2);[idx,C] = kmedoids(y,nc);
% Cl = {'[0 0 0]','[1 0 1]','[0 1 0]','[0 0 0.5]','[1 0 0]'...
%      '[0 0 1]'};%,'[0 0 1]','[0.5 0 0.5]','[0.5 0 0]'};
% n=length(S1);k = 1:n;
% h41=figure()
% for i=1:nc
%     text(y(k(idx==i),1),y(k(idx==i),2),num2str(find(idx==i)),'FontSize',12,'color',Cl{i});
% end
% xlim([min(y(:,1))-0.02, max(y(:,1))+0.02]); ylim([min(y(:,2))-0.02, max(y(:,2))+0.02]);
% if country==1;title(['Medoids_USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% if country==2;title(['Medoids_MJPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% set(gca, 'FontSize',12,'FontWeight','bold');
% xlabel('x-coordinates');ylabel('y-coordinates');figurepalette('show')
% %%%%%%%%%%%%% cumulative frequeny %%%%%%%%%%%%%%%5555
%B1 = changem(idx1,[2 4 3 1],[1 2 3 4]);
if nc==4 && country==1 && eps==0.65 && wind==20 && shift==20
   idx1 = changem(idx,[1 2 3 4],[1 3 4 2]); 
   fprintf('fffkiff')
end
if nc==4 && country==1 && eps==0.65 && wind==20 && shift==10
   idx1 = changem(idx,[1 2 3 4],[1 4 3 2]); 
   fprintf('fsdfa')
end
if nc==4 && country==2 && eps==0.65 && wind==20 && shift==10
   idx1 = changem(idx,[1 2 3 4],[2 1 3 4]); 
   fprintf('fsdfa')
end
if nc==4 && country==2 && eps==0.65 && wind==20 && shift==20
   idx1 = changem(idx,[1 2 3 4],[2 1 4 3]); 
   fprintf('fffkiff')
end
if nc==5 && country==1 && eps==0.65 && wind==20 && shift==10
   idx1 = changem(idx,[1 2 3 4 5],[5 1 4  2]); 
end
if nc==5 && country==2 && eps==0.65 && wind==20 && shift==20
   idx1 = changem(idx,[1 2 3 4 5],[2 5 1 4 3]); 
end

 y1=idx1;%    y1=idx;
h2=figure();
plot(y1, '.','MarkerSize',10);ylim([0.5,6.5])
for i=1:length(idx)
    if mod(i,2)==0
        h=text(i,y1(i)+0.1,num2str(i),'FontSize',8);set(h,'Rotation',90);
    else
        h=text(i-0.3,y1(i)-0.2,num2str(i),'FontSize',8); set(h,'Rotation',90);
    end
end
set(gca, 'FontSize',12,'FontWeight','bold');
xlabel('Epoch');ylabel('State');
if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
figurepalette('show')
% set(gca, 'FontSize',8,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nbin=linspace(1,nc,nc);
% h3=figure();
% hist(y1,(1:nc));
% stem(h3,'b');xlim([0,nc+0.5]);
% set(gca, 'FontSize',12,'FontWeight','bold');
% xlabel('State');ylabel('Frequency');
% if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% figurepalette('show')
%  set(gca, 'FontSize',8,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % T11=find_patterns(idx);
% % for i=1:30
% %     fprintf('%g%g',T11{i,1},T11{i,2})
% %     fprintf('\n');
% % end
for i=1:nc
    for j=1:nc
        pattern=[i j];
        indices=findPattern2(y1,pattern);
        fprintf('%g\t%g\t%g\n',i,j,size(indices,1))
        freq1(i,j)=size(indices,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h4=figure();b=0;
b=bar3(freq1,0.6);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
for i=1:nc
    for j=1:nc
        text(j+0.4,i+0.4,freq1(i,j),num2str(freq1(i,j)),'HorizontalAlignment','left','Color','m','FontSize',8,'FontWeight','bold')
    end
end
set(gca, 'FontSize',12,'FontWeight','bold');
xlabel('Second MS');ylabel('First MS');zlabel('Frequency of paired states')
if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
view(50,40)
figurepalette('show')
%%%%%%%%%%%%%%%%%%%%%%%%%% probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h5=figure();b=0;
for i=1:size(freq1,1)
    for j=1:size(freq1,2)
        summ=sum(freq1(i,:));
        freq2(i,j)=freq1(i,j)/summ;
    end
end
b=bar3(freq2,0.6);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
for i=1:nc
    for j=1:nc
        text(j+0.4,i+0.4,freq2(i,j),num2str(round(freq2(i,j),2)),'HorizontalAlignment','left','Color','m','FontSize',8,'FontWeight','bold')
    end
end
set(gca, 'FontSize',12,'FontWeight','bold');
xlabel('Second MS');ylabel('First MS');zlabel('Probability of paired states')
if country==2;title(['JPN, ','Window= ',num2str(wind),', shift= ',num2str(shift)]);end
if country==1;title(['USA, ','Window= ',num2str(wind),', shift= ',num2str(shift)]);end
view(50,40)
figurepalette('show')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h6=figure();
% imagesc(freq1);
% set(gca, 'FontSize',12,'FontWeight','bold');
% xlabel('State');ylabel('State');colorbar
% if country==1;title(['USA (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% if country==2;title(['JPN (wind= ',num2str(wind),', shift=', num2str(shift),', eps=', num2str(eps),')']);end
% figurepalette('show')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=1:length(mean_all);
h7=figure(); plot(x1,mean_all,'-.');grid on;
xlabel('Frames');ylabel('Mean Correlation');
if country==2;title(['JPN, ','Window= ',num2str(wind),', shift= ',num2str(shift)]);end
if country==1;title(['USA, ','Window= ',num2str(wind),', shift= ',num2str(shift)]);end
for i=1:length(mean_all);text(x1(i),mean_all(i)+0.01,num2str(i),'FontSize',8);end
mean_all_eps=mean(mean(r4));hold on
plot(x1,mean_all_eps(:),'r');hold off
set(gca, 'FontSize',12,'FontWeight','bold');figurepalette('show')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% using idx probability of 5 epoch %%%%%%%%%%%
wind_len=5; clear idx22 Prob aa idx11;I1=0;
for i=1:wind_len:length(idx1)-wind_len
    I1=I1+1;
    %% frames of length wind_len(eg. 1-5,6-10,11-15...)
    aa= idx1(i:i-1+wind_len);
    idx22=0;
    %% calculating freq in length wind_len
    for j=1:nc
        idx11=sum(aa==j);
        idx22(j)=idx11/wind_len; 
    end
%     fprintf('%g\n%g\n', aa,idx22)
%% printing line at the last point of wind_len
    Prob(:,i-1+wind_len)=idx22;  
end
figure()
h8=bar(Prob',2, 'stacked');
if country==1;set(h8,{'FaceColor'},{'m';'g';'r';'b'});
legend('I','II','III','IV');end
if country==2;set(h8,{'FaceColor'},{'m';'g';'y';'r';'b'});
legend('I','II','III','IV','V');end

xlabel('Frames');ylabel('Probability');
if country==2;title(['JPN, ','Window= ',num2str(wind),', shift= ',num2str(shift), ', Number of epochs= ',num2str(wind_len)]);end
if country==1;title(['USA, ','Window= ',num2str(wind),', shift= ',num2str(shift), ', Number of epochs= ',num2str(wind_len)]);end
set(gca, 'FontSize',12,'FontWeight','bold');figurepalette('show')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h9=figure()
% subplot(221);plot(Prob(1,:),'m.-');
% xlabel('Frames');ylabel('Probability');legend('I state');
% set(gca, 'FontSize',8,'FontWeight','bold');figurepalette('show')
% subplot(222);plot(Prob(2,:),'g.-');
% xlabel('Frames');ylabel('Probability');legend('II state')
% set(gca, 'FontSize',8,'FontWeight','bold');figurepalette('show')
% subplot(223);plot(Prob(3,:),'r.-');
% xlabel('Frames');ylabel('Probability');legend('III state')
% set(gca, 'FontSize',8,'FontWeight','bold');figurepalette('show')
% subplot(224);plot(Prob(4,:),'b.-');
% xlabel('Frames');ylabel('Probability');legend('IV state')
% set(gca, 'FontSize',8,'FontWeight','bold');figurepalette('show')
% if country==2;title(['JPN, ','Window= ',num2str(wind),', shift= ',num2str(shift), ', Number of epochs= ',num2str(wind_len)]);end
% if country==1;title(['USA, ','Window= ',num2str(wind),', shift= ',num2str(shift), ', Number of epochs= ',num2str(wind_len)]);end
% set(gca, 'FontSize',8,'FontWeight','bold');figurepalette('show')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h8=figure();
% mean_all_eps=mean(mean(r4));
% plot(mean_all,'b')
% hold on
% plot(mean_all_eps(:),'r')
% hold off
% xlabel('Frames');ylabel('Mean Correlation(eps=0 "Blue")(eps=0.65 "Red")');
% if country==2;title(['JPN, ','Time window and shift= ',num2str(wind)]);end
% if country==1;title(['USA, ','Time window and shift= ',num2str(wind)]);end
% if Input==3;cd('plots/msd1/')
%     % saveas(h1,['msd_Frame',num2str(t)],'jpg');
%     % saveas(h21,['r1_Frame',num2str(t)],'jpg');
%     saveas(h22,['d1_Frame',num2str(t)],'jpg');
%     saveas(p1,['Kmeans',num2str(t)],'jpg');
%     cd ('../../')
% end
% cd('/..');
% plot(idx,'o');ylim([0,10]);
%   intracluster;
%  a_SVD;
save_input=212; saveplot;
if country==1;fprintf('\n USA stock list done\n');end
if country==2;fprintf('\n JPN stock list done\n');end

%%%%%%%%%%%%%%%%%%%%%%%%%% kiran code
% X = cmdscale(D,3);nc = 8; % numbe of clusters
% [idx,C] =  kmeans(X,nc);figure;
% %plot(X(:,1),X(:,2),'k.','MarkerSize',6);
% Cl = {'[0 0 0]','[1 0 1]','[1 1 0]','[0 0 0.5]','[1 0 0]'...
%     '[0 1 0]','[0 0 1]','[0.5 0 0.5]','[0.5 0 0]'};
% for i=1:nc
%     plot3(X(idx==i,1),X(idx==i,2),X(idx==i,3),'color',Cl{i}, ...
%         'marker','.','MarkerSize',4, 'LineStyle','none');
%     hold on;
% end
% grid on;
%   %xlim([-0.6 0.6]); ylim([-0.6 0.6]);
%  set(gca, 'FontSize',20,'FontWeight','bold');
%  xlabel('x coordinates');
%   ylabel('y coordinates');
%   hold on;
%   n=length(D);
%   k = 1:n;
%   for i=1:nc
%         text(X(k(idx==i),1),X(k(idx==i),2),X(k(idx==i),3),num2str(find(idx==i)),'FontSize',10,'color',Cl{i});
%         hold on;
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  Both for wind=10 and shift = 0
% if Input1==2
% t1=[1,27,52,77,103,128,153,178,204,229,254,279,306,330,355,381,406,431,456,481,506,531,556,582,607,632,657,682,707,733,758,783];
% str1=[1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016];
% end
% if Input1==1
% t1=[1,31,58,85,113,137,162,187,211,236,261,286,310,335,359,384,409,433,458,482,507,532,556,581,605,630,654,679,703,728,752,777];
% str1=[1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016];
% end