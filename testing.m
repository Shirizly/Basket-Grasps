clearvars
close all
%% Define the polygon
% PG = Polygon([1,0;2,0;2,3;-1,4;-2,3;-1,1;-0.5,1].',[0;1.5]);
PG = Polygon([-3,0;-1.5,-5;3,-4;3,4;-1.5,5].',[0;0]);
% PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].',[0;0]);
PG.drawPolygon();
%% Find equilibrium grasps
res = 50;
[PG,S,X,VL] = PG.findBdyVariable(res);
PG.res = res;
tic
[s1,s2,linesegs] = PG.Eqcheck();
toc
%% Check grasps for minimum/maximum
status = PG.Statuscheck(s1,s2,S,X,VL);
%% Full contact space figure
figure;
hold on
c_layers=50;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers)
xlabel('s_1')
ylabel('s_2')
colorbar
col = [0 1 0; 0 0 0; 1 0 0];
for i = 1:length(s1)
plot(s1(i),s2(i),'.','MarkerEdgeColor',col(status(i)+2,:),'MarkerFaceColor',col(status(i)+2,:))
end
% for i=1:size(linesegs,1) % plotting the areas defined by horizontal edges
%     plot(linesegs(i,[1,3]),linesegs(i,[2,4]),'k','linewidth',2)
%     plot(linesegs(i,[1,3]),linesegs(i,[2,4]),'k.')
%     plot(linesegs(i,[2,4]),linesegs(i,[1,3]),'k','linewidth',2)
%     plot(linesegs(i,[2,4]),linesegs(i,[1,3]),'k.')
% end
%%
% lineseg = [];
% PG.drawPolygonS(s1,s2,S,X,lineseg);

