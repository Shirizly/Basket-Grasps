s1 = 24.7981;
s2 = [];
theta = -0.6550;% -pi;
f2 = [-4.25;0];
f1 = [4.25;0];

if isempty(s2)
    axle = PG.get('1Pos',s1);
    d = f1-axle;
else
    axle = PG.get('1Pos',s2);
    d = f2-axle;
end
fig = figure();
PG.drawPolygonMoved(axle,theta,d,'k');
PG.drawPolygonNode(theta,s1,s2,f1,f2);