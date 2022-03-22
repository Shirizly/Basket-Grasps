function [fig,dot] = drawPolygon2fGraspRange(PG,BGSegEnds,varargin)
% s1 = 1.184;
% s2 = [];
theta = 0;
% f2 = [0;0];
% f1 = [-0.2;-0.03];
switch nargin
    case 2
        fig = PG.drawPolygon();
        col = 'g';
    case 3
        fig = PG.drawPolygon();
        col = varargin{1};
    case 4
        fig = varargin{1};
        col = varargin{2};
end
    


dot = [];

hold on
for i=1:2
f1(:,i) = PG.get('1Pos',BGSegEnds(1,i));
f2(:,i) = PG.get('1Pos',BGSegEnds(2,i));
end
plot(f1(1,:),f1(2,:),'Color',col,'lineWidth',6)
plot(f2(1,:),f2(2,:),'Color',col,'lineWidth',6)
if norm(f1(:,1)-f1(:,2))<5E-2
    dot = plot(f1(1,:),f1(2,:),'.','Color',col,'lineWidth',20,'markersize',24);
end
if norm(f2(:,1)-f2(:,2))<5E-2
    dot = plot(f2(1,:),f2(2,:),'.','Color',col,'lineWidth',20,'markersize',24);
end
set(fig,'renderer','painters')
end