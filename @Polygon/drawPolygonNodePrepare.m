function drawPolygonNodePrepare(PG,s1,s2,theta,f1,f2)
% s1 = 1.184;
% s2 = [];
% theta = -1.809;% -pi;
% f2 = [0;0];
% f1 = [-0.2;-0.03];


if isempty(s2)
    axle = PG.get('1Pos',s1);
    d = f1-axle;
else
    axle = PG.get('1Pos',s2);
    d = f2-axle;
end
PG.drawPolygonMoved(axle,theta,d,'k');
PG.drawPolygonNode(theta,s1,s2,f1,f2);
end