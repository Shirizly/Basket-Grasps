function fig = drawFingerObstacle(PG,origin,fingerToOrigin,varargin)
switch nargin
    case 3
        fig = figure();
        col = [0,0,1];
        xlabel('X');
        ylabel('\theta');
        zlabel('y');
    case 4
        fig = figure();
        col = varargin{1};
        xlabel('X');
        ylabel('\theta');
        zlabel('y');
    case 5
        fig = varargin{2};
        col = varargin{1};
end
deltath = 0.05;
tt = -pi()/2:deltath:pi/2;
verta = PG.vertex;
verta = [verta,verta(:,1)];
theta = wrapToPi(tt(1));
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
prevert = origin + R*fingerToOrigin+verta;
hold on
for th = tt(2:end)
    theta = wrapToPi(th);
    R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    vert = origin + R*fingerToOrigin+verta;
    for i = 1:(size(vert,2)-1)
        xx = [vert(1,i),vert(1,i+1),prevert(1,i+1),prevert(1,i)];
        yy = [vert(2,i),vert(2,i+1),prevert(2,i+1),prevert(2,i)];
        zz = [th,th,th-deltath,th-deltath];
         patch(xx,zz,yy,col,'linestyle','none')%,'edgeColor',col)
    end
%     plot3(vert(1,:),th*ones(1,length(vert)),vert(2,:),'k','linewidth',1)
    prevert = vert;
end
% [xx,tt,yy] = cobstacle(PG,finger);
% mesh(xx,tt,yy)
view(3)
camlight('headlight')

end
