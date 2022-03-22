theta = -pi()/2;
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
PG = Polygon(R*[-2,-2;2,-2;2,2;1,2;0,0;-3,2].'+[0;1],[0;0.5])
plot = PG.drawPolygonSS();
set(gca,'visible','off')
set(plot,'renderer','painters')