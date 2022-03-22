% preparing figures for the paper
%% fig 1
s1 = 19.24;
s2 = 7.242;
theta = 0;
length = 6;
scale = 1;
drawPolygon2fGrasp(PG,s1,s2,theta)
hold on
x = -2;
y = -0.9247;
xlim([-2.5,length+2])
quiver(x,y,0,2*scale,'k');
quiver(x,y,2*scale,0,'k');
text(x+1.6*scale,y+0.15*scale,'x','Interpreter','latex','fontsize',18)
text(x+0.12*scale,y+1.8*scale,'y','Interpreter','latex','fontsize',18)

rectangle('Position',[x+2,y-0.2*scale,length,0.2*scale]);

for i = 1:length/scale
    xstart = x+2+(i-1)*scale;
    ystart = y-0.2*scale;
    plot([xstart,xstart+scale],[ystart,ystart+0.2*scale],'k');
end
% text(0,0,'$\mathcal{B}$','Interpreter','latex','fontsize',18);
% set(gca,'visible','off')

%% fig 2
figure;
set(gcf,'renderer','Painters');
hold on
c_layers=20;cwidth=1.5;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers,'LineWidth',cwidth)
xlabel('$s_1$','interpreter','latex','fontsize',18)
ylabel('$s_2$','interpreter','latex','fontsize',18)
colorbar
print(gcf,'-depsc','-painters','fig2b.eps');




%% fig 2 (b)
figure;
set(gcf,'renderer','Painters');
hold on
c_layers=8;cwidth=1.5;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers,'LineWidth',cwidth)
xlabel('$s_1$','interpreter','latex','fontsize',20)
ylabel('$s_2$','interpreter','latex','fontsize',20)
h=colorbar % remove colorbar from non-rightmost part of figure
set(get(h,'label'),'string','$\:\sigma$','interpreter','latex','FontSize',20,'rotation',0);

col = [0 1 0; 0 0 0; 1 0 0];
s1ends = [5.481,8.103,1.01,2,2,3];
s2ends = [20.53,18.62,14.56,14.56,11.73,11.73];
% for i = 1:length(s1)
%     if status(i) ~=-1
% plot(s1(i),s2(i),'.','MarkerEdgeColor',col(status(i)+2,:),'MarkerFaceColor',col(status(i)+2,:),'markersize',8)
%      end
% end
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.g','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.g','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'g','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'g','lineWidth',4)
end
s1ends = [4,4,...
            4.5,8.103,...
            8.103,8.103,...
            0,2,...
            2,4,...
            4,4,...
            0,2];
s2ends = [17.74,20.53,...
            17.74,17.74,...
            18.62,20.53,...
            17.74,17.74,...
            8.103,8.103,...
            8.103,11.73,...
            4,4];
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.k','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.k','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'k','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'k','lineWidth',4)
end
set(gcf,'Position',[300,400,360,300])
% savefig(gcf,[subFolderName '\fig4a.fig']);
print(gcf,'-depsc','-painters',[subFolderName '\fig4a.eps']);
% epsclean('fig4a.eps');

%% fig 3 (a)
figure;
hold on
% Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers,'LineWidth',cwidth)
contour(S,S,Sigma,[6,6],'k','LineWidth',4)
xlabel('$s_1$','interpreter','latex','fontsize',20)
ylabel('$s_2$','interpreter','latex','fontsize',20)
h=colorbar % remove colorbar from non-rightmost part of figure
set(get(h,'label'),'string','$\:\sigma$','interpreter','latex','FontSize',20,'rotation',0);
col = [0 1 0; 0 0 0; 1 0 0];
s1ends = [7.1,8.103];
s2ends = [19.35,18.62];
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.r','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.r','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'r','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'r','lineWidth',4)
end
s1ends = [5.481,7.1,1.01,1.99,2.01,2.98];
s2ends = [20.53,19.35,14.56,14.56,11.73,11.73];
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.g','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.g','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'g','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'g','lineWidth',4)
end
set(gcf,'Position',[300,400,375,300])
% savefig(gcf,[subFolderName '\fig4b.fig']);
print(gcf,'-depsc','-painters','\fig4b.eps');

%% fig 3 (b)
figure;
hold on
Sigma=inter_finger_distance(X,X);

% rectangle('Position',[0 0 4 PG.S(end)],'FaceColor',[0.3,0.3,0.3],'lineStyle','none')
% rectangle('Position',[4 0 PG.S(end)-4 4],'FaceColor',[.3,.3,.3],'lineStyle','none')
contour(S,S,Sigma,c_layers,'LineWidth',cwidth)
xlabel('$s_1$','interpreter','latex','fontsize',20)
ylabel('$s_2$','interpreter','latex','fontsize',20)
h=colorbar % remove colorbar from non-rightmost part of figure
set(get(h,'label'),'string','$\:\sigma$','interpreter','latex','FontSize',20,'rotation',0);
col = [0 1 0; 0 0 0; 1 0 0];
s1ends = [5.481,8.103];
s2ends = [20.53,18.62];
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.g','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.g','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'g','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'g','lineWidth',4)
end
s1ends = [1.01,1.99,2.01,2.98];
s2ends = [14.56,14.56,11.73,11.73];
for i = 1:2:length(s1ends)
   plot(s1ends(i:i+1),s2ends(i:i+1),'.r','markerSize',8)
   plot(s2ends(i:i+1),s1ends(i:i+1),'.r','markerSize',8)
   plot(s1ends(i:i+1),s2ends(i:i+1),'r','lineWidth',4)
   plot(s2ends(i:i+1),s1ends(i:i+1),'r','lineWidth',4)
end
rectangle('Position',[0 0 4 PG.S(end)],'FaceColor',[0,0,0,0.4],'lineStyle','none')
rectangle('Position',[4 0 PG.S(end)-4 4],'FaceColor',[0,0,0,0.4],'lineStyle','none')
set(gcf,'Position',[300,400,375,300])
% savefig(gcf,[subFolderName '\fig4c.fig']);
print(gcf,'-depsc','-painters','\fig4c.eps');
%% fig 6


%% fig bottle
figure
f2 = [4;2];
f1 = [0;0];
s1 = 46.2041;
s2 = 23.4386;
theta = 0;
basepos = [f1,f2];
drawPolygon2fGrasp(PG,s1,s2,theta)
set(gca,'visible',0)
set(gcf,'renderer','painters')

%% fig gun
figure
s1 = 15.3010;
s2 = 37.2875;
theta = 0.0008;

drawPolygon2fGrasp(PG,s1,s2,theta)
set(gca,'visible',0)
set(gcf,'renderer','painters')