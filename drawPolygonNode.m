function fig = drawPolygonNode(PG,s1,theta,f1,f2)
% s1 = 5.441;
s2 = [];
% theta = 0.3218 -pi;
% f2 = [1;0];
% f1 = [-2;1];


if isempty(s2)
    axle = PG.get('1Pos',s1);
    d = f1-axle;
end
fig = drawPolyGraspMoved(PG,axle,theta,d,s1,s2,f1,f2);
end