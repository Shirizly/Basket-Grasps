%% Define the polygon
clearvars
clc
type = 4;
switch type
    case 0 % pentagon
        R = [0 1;-1 0];
        theta = -pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
%         PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        load([pwd '/' subFolderName '/AlgResult.mat'],'U','delU','pp','PG')
        f2 = [-4.4;1.8];
        f1 = [4.4;1.8];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,3];
        U = U*0.995;
    case 1 % cup
       
%         PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results
       load([pwd '/' subFolderName '/AlgResult.mat'],'U','delU','pp','PG','basepos')
       f2 = basepos(:,2);
       f1 = basepos(:,1);
       basepos = [f1,f2];
       res = 50;
       starting_node = [2,6];
       U = U*0.995;
    case 2 % bottle
        load object5.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [0.06;0.01];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 1500;
        starting_node = [2,21];
    case 3 % Elon
        load object.mat
        %     object = [object(:,7:8) object(:,1:6)];
        com = [0.045,0.07].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Elon'; % where to save results
    case 4 % Crown
        object = [-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 20;
         starting_node = [2,4];    
    case 5 % new bottle
        load NewBottle.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
%         theta =  1.515;
%         R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
%         object = R*object;
        com(2) = max(object(2,:))-com(2);
        object(2,:) = max(object(2,:))-object(2,:);
        object(2,1) = 0;
        object = round(object./20,1);
        com = round(com./20,0);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [4;2];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,24];
        case 6 % Gun
%         load gun2.mat
%         object = Round_object;
%         object = 0.8*[1,10;26,10;21,1;29,9;31,2;28,0;32,0;32,12;35,16;32,18;33,16;32,14;31,14;31,16;1,16;0.5,13.5;0,13;0.5,12.5].';
        
        object = 0.8*[1,10;26,10;21,1;29,9;32,0;32,12;36,12;40,8;41,9;35,14;31,16;1,16].';%32,14;31,14
        com = 0.8*[27;10];
        theta =  0.423;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
         object = R*object;
         com = R*com;
%          object = round(object./20,1);
%          com = round(com./20,0);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Gun'; % where to save results
        f2 = [18.2;16.09];
        f1 = [11.4;13.9];
        basepos = [f1,f2];
        res = 30;
        starting_node = [2,6];
end

%% draw the polygon (only once)
if ~exist('polyDrawing','var')
    polyDrawing = cell(0);
end
if isempty(polyDrawing)||~isgraphics(polyDrawing)
    polyDrawing = PG.drawPolygon();
end
%% draw the object configuration space c-obstacles:
cs = drawCObstacle(PG,f1,'b');
cs = drawCObstacle(PG,f2,'c',cs);
%% generate and combine the c-obstacle's ceiling:
deltath = 0.005;
tt = -pi():deltath:pi;
comTOvert = PG.vertex-PG.com;
comTOvertDist = comTOvert.'*comTOvert;
maxdist = sqrt(max(diag(comTOvertDist)));
fsort = sort([f1(1),f2(1)]);
xx = linspace(fsort(1)-maxdist,fsort(2)+maxdist,5000);

[~,~,yy1] = cobstacles(PG,f1,xx,tt);
[~,~,yy2] = cobstacles(PG,f2,xx,tt);
yy = yy1.*(yy1>=yy2) + yy2.*(yy2>yy1);
yy(yy==-500)=nan;
%% generate and draw the rim:
Rim = getRim(xx,tt,yy,U);
Rim = [Rim;1.01*U*ones(size(Rim(1,:)))];
plot3(Rim(1,:),Rim(2,:),Rim(3,:),'r','linewidth',3)
xlim([2*min(Rim(1,:)),2*max(Rim(1,:))])
ylim([2*min(Rim(2,:)),2*max(Rim(2,:))])
zlim([2,U*1.5])

%% plot the grasp and the puncture point:
plot3(0,0,(U/0.99-delU)*1.01,'+g','markersize',12,'linewidth',4)
pos = PG.get('comLocationDS',basepos,pp(1),pp(2));
plot3(pos(1),pos(3),pp(3)*1.05,'oy','markersize',4,'linewidth',4)
% surf(xx,tt,yy2)

%% draw the hand configuration space c-obstacles, 
% with the hand's origin at the relative position the CoM is when the object is in the basket grasp:
BG = [0;U-delU];
fr1 = BG-f1;
fr2 = BG-f2;
hs = drawFingerObstacle(PG,fr1,'b');
hs = drawFingerObstacle(PG,fr2,'c',hs);
% transform points in object c-space to their equivalent hand c-space
% positions:
for i = 1:length(Rim)
     HRim(:,i) = obpTohp(BG(2),Rim(:,i));
end
hold on
plot3(HRim(1,:),HRim(2,:),HRim(3,:),'r','linewidth',3)
HBG = obpTohp(BG(2),[BG(1);0;BG(2)]);
plot3(HBG(1),HBG(2),HBG(3),'+g','markersize',20,'linewidth',8)
HPP = obpTohp(BG(2),[pos(1);pos(3);pp(3)]);
plot3(HPP(1),HPP(2),HPP(3),'oy','markersize',8,'linewidth',4)
xlim([2*min(HRim(1,:)),2*max(HRim(1,:))])
ylim([2*min(HRim(2,:)),2*max(HRim(2,:))])
zlim([-3,1])

function hp = obpTohp(h,obp)
    xr = obp(1);
    thetar = obp(2);
    yr = obp(3);
    pos = [xr;yr-h];
    theta = -thetar;
    R =  [cos(theta) -sin(theta);sin(theta) cos(theta)];
    pos = R*pos;
    hp = [-pos(1),-thetar,-pos(2)];
end