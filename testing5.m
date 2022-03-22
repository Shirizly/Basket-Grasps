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
%         f2 = [-4.4;1.8];
%         f1 = [4.4;1.8];
        f2 = [8.3223;-1.3181];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,3];
        BG = [0;0;(U-delU)];
        U = U*0.995;
        
    case 1 % cup
       
%         PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results
       load([pwd '/' subFolderName '/AlgResult.mat'],'U','delU','pp','PG','basepos')
       f2 = basepos(:,2);
       f1 = basepos(:,1);
%        basepos = [f1,f2];
       res = 50;
       starting_node = [2,6];
       BG = [0;0;(U-delU)];
%        U = (U-delU)*0.995+delU;
%        pp = [0.2516;1.159;5.045];
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
% if ~exist('polyDrawing','var')
%     polyDrawing = cell(0);
% end
% if isempty(polyDrawing)||~isgraphics(polyDrawing)
%     polyDrawing = PG.drawPolygon();
% end
%% draw the object configuration space c-obstacles:

% generate and combine the c-obstacle's ceiling:
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
% surf(xx,tt,yy)
% generate and draw the rim:
Rim = getRim(xx,tt,yy,U*1);
Rim = [Rim;1*U*ones(size(Rim(1,:)))];
if length(pp)>4
    pos = PG.get('comLocationDS',basepos,pp(1),pp(2));
    pp = [pos(1);pos(3);pp(3)];
end
if length(pp)>3
pos = PG.get('comLocationSS',basepos(:,pp(end)),pp(1),pp(3));
pp = [pos(1);pp(3);pos(2)];
end
for i=1:1
hold on
cs = drawCObstacle(PG,f1,'b');
cs = drawCObstacle(PG,f2,'c',cs);
plot3(Rim(1,:),Rim(2,:),Rim(3,:),'r','linewidth',3)
xlim([2*min(Rim(1,:)),2*max(Rim(1,:))])
ylim([2*min(Rim(2,:)),2*max(Rim(2,:))])
zlim([2,U*1.5])
% Urange = linspace(U-delU,U,5);
% for Ul = Urange(1:end)
%     Riml = getRim(xx,tt,yy,Ul);
%     Riml = [Riml;U*ones(size(Riml(1,:)))];
%     plot3(Riml(1,:),Riml(2,:),Riml(3,:),'.k','linewidth',3)
% end
% plot the grasp and the puncture point:

plot3(BG(1),BG(2),BG(3)*1.01,'+g','markersize',20,'linewidth',6)

plot3(pp(1),pp(2),pp(3)*1.01,'oy','markersize',4,'linewidth',4)
% surf(xx,tt,yy2)
end
%% draw the hand configuration space c-obstacles, 
% with the hand's origin at either of the fingers:
% HO = BG([1,3]);
HO = f1;
fr1 = HO-f1;
fr2 = HO-f2;

% transform points in object c-space to their equivalent hand c-space
% positions:
for i = 1:length(Rim)
     HRim(:,i) = obpTohp(HO,BG,Rim(:,i));
end
for i=1:0
hold on
hs = drawFingerObstacle(PG,[0;0],fr1,'b');
hs = drawFingerObstacle(PG,[0;0],fr2,'c',hs);
plot3(HRim(1,:),HRim(2,:),HRim(3,:),'.r','linewidth',3)
HBG = obpTohp(HO,BG,BG);
plot3(HBG(1),HBG(2),HBG(3),'+g','markersize',20,'linewidth',8)
HPP = obpTohp(HO,BG,pp);
plot3(HPP(1),HPP(2),HPP(3),'oy','markersize',8,'linewidth',4)
% xlim([2*min(HRim(1,:)),2*max(HRim(1,:))])
% ylim([2*min(HRim(2,:)),2*max(HRim(2,:))])
% zlim([BG(2)-1,U+1])
% Urange = linspace((U-delU)*1.01,U,10);
% for Ul = Urange(1:end)
%     Riml = getRim(xx,tt,yy,Ul);
%     Riml = [Riml;U*ones(size(Riml(1,:)))];
%     for i = 1:length(Riml)
%      HRiml(:,i) = obpTohp(HO,BG,Riml(:,i));
%     end
% plot3(HRiml(1,:),HRiml(2,:),HRiml(3,:),'.k','linewidth',3)
% end
end
%%

fig = PG.drawPolygon();
hold on
HO = f1;
fr1 = HO-f1;
fr2 = HO-f2;

% transform points in object c-space to their equivalent hand c-space
% positions:
for i = 1:length(Rim)
     HRim1(:,i) = obpTohp(HO,BG,Rim(:,i));
end
hold on
HBG1 = obpTohp(HO,BG,BG);

plot(HRim1(1,:),HRim1(3,:),'r','linewidth',4)
plot(HBG1(1),HBG1(3),'ok','markersize',6,'linewidth',3)

HO = f2;
fr1 = HO-f1;
fr2 = HO-f2;
for i = 1:length(Rim)
     HRim2(:,i) = obpTohp(HO,BG,Rim(:,i));
end
HBG2 = obpTohp(HO,BG,BG);
plot(HRim2(1,:),HRim2(3,:),'r','linewidth',4)
plot(HBG2(1),HBG2(3),'ok','markersize',6,'linewidth',3)
plot([HBG1(1),HBG2(1)],[HBG1(3),HBG2(3)],'k')
%%
% theta = pi()/5;
% fn1 = [-1.5;0.6];
% fn2 = fn1+4.5432*[cos(theta);sin(theta)];
% plot([fn1(1),fn2(1)],[fn1(2),fn2(2)],'b')
% plot(fn1(1),fn1(2),'ob','markersize',6,'linewidth',3);
% plot(fn2(1),fn2(2),'ob','markersize',6,'linewidth',3);
%%
% HO = f1;
% fl = f2-f1;
% theta = theta-atan2(fl(2),fl(1));
% hp = [fn1(1),theta,fn1(2)];
% obp = hpToobp(HO,BG,hp);
% drawPolygonMoved(PG,PG.com,obp(2),obp([1,3]).'-PG.com,'b')
% range = round(linspace(1,4201,5));
% range = [1250,1100,1400];
range = [700,906,1000];
for i=1:length(range)
    fig = PG.drawPolygon();hold on;
    plot(HBG1(1),HBG1(3),'ok','markersize',6,'linewidth',3)
    plot(HBG2(1),HBG2(3),'ok','markersize',6,'linewidth',3);
    plot([HBG1(1),HBG2(1)],[HBG1(3),HBG2(3)],'k');
    obp = Rim(:,range(i));
    drawPolygonMoved(PG,PG.com,obp(2),obp([1,3])-PG.com,'b');
    axle = PG.com;
    d = obp([1,3])-PG.com;
    theta = obp(2);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    movedRim1 = axle + R*(HRim1([1,3],:)-axle)+d;
    plot(movedRim1(1,:),movedRim1(2,:),'r','linewidth',4)
    movedRim2 = axle + R*(HRim2([1,3],:)-axle)+d;
    plot(movedRim2(1,:),movedRim2(2,:),'r','linewidth',4)
    saveas(fig,['positionRimDS' num2str(i) '.bmp'])
    close all
end
function hp = obpTohp(HO,BG,obp)
    d0 = [BG(1);BG(3)];
%     theta0 = BG(2);
    X0 = HO;
    thetaO = obp(2);
    d = [obp(1);obp(3)];
    thetar = -thetaO;
    R =  [cos(thetar) -sin(thetar);sin(thetar) cos(thetar)];
    pos = R*(X0-d)+d0;
    hp = [pos(1),thetar,pos(2)];
end
function obp = hpToobp(HO,BG,hp)
    d0 = [BG(1);BG(3)];
    theta0 = BG(2);
    X0 = HO;
    thetar = hp(2);
    d = [hp(1);hp(3)];
    thetao = 0*theta0-thetar;
    R =  [cos(thetao) -sin(thetao);sin(thetao) cos(thetao)];
    pos = R*(d0-d)+X0;
    obp = [pos(1),thetao,pos(2)];
end