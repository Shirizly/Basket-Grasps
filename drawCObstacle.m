function fig = drawCObstacle(PG,finger,varargin)
switch nargin
    case 2
        fig = figure();
        col = [0,0,1];
        xlabel('X');
        ylabel('\theta');
        zlabel('y');
    case 3
        fig = figure();
        col = varargin{1};
        xlabel('X');
        ylabel('\theta');
        zlabel('y');
    case 4
        fig = varargin{2};
        col = varargin{1};
end
deltath = 0.01;
tt = -pi():deltath:pi;
comTOvert = PG.vertex-PG.com;
comTOvert = [comTOvert,comTOvert(:,1)];
prevert = finger+comTOvert;
hold on
for th = tt(2:end)
    theta = wrapToPi(th+pi());
    R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    vert = finger+R*comTOvert;
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
camlight('headlight','infinite')

end
