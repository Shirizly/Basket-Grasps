%% Define the polygon
clear all
% PG = Polygon([1,0;2,0;2,3;-1,4;-2,3;-1,1;-0.5,1].',[0;1.5]);
PG = Polygon([-1,-1;2,-1;3,-3;3,3;1,1;0,2;-1,1;-2,3;0.2,3;-2,5;-3,1;-2,-4].',[0;0.5]);
% PG = Polygon([-3,0;-1.5,-5;3,-4;3,4;-1.5,5].',[0;0]);
% PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].',[0;0]);
PG.drawPolygon();
%% Find equilibrium grasps
res = 50;
[PG,S,X,VL] = PG.findBdyVariable(res);
figure
Sigma=inter_finger_distance(X,X);
c_layers=50;
contour(S,S,Sigma,c_layers)

%%
hold on
for i =1:length(PG.VL)-1
    
    x1 = [i,S(end)];
    y1=  [0,S(end)-i];
    plot(x1,y1,'k')
    y2 = [i,S(end)];
    x2 = [0,S(end)-i];
    plot(x2,y2,'k')
end
for i =1:3:round(S(end))
    x1 = [i,S(end)];
    y1=  [S(end),i];
    plot(x1,y1,'k')
    y2 = [i,0];
    x2 = [0,i];
    plot(x2,y2,'k')
end