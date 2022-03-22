clearvars
close all
%% Define the polygon
% %% Define the polygon
clear all
clc
object = 2;
switch object
    case 0 % pentagon
        R = [0 1;-1 0];
        PG = Polygon(R*[-3,0;-1.5,-5;3,-4;3,4;-1.5,5].',[0;0]);
        subFolderName = 'Examples for research proposal\Convex polygon'; % where to save results
        f2 = [-4.25;0];
        f1 = [4.25;0];
        basepos = [f1,f2];
        res = 50;
    case 1 % cup
        % PG = Polygon([1,0;2,0;2,3;-1,4;-2,3;-1,1;-0.5,1].',[0;1.5]);
        % PG = Polygon([-1,-1;2,-1;3,-3;3,3;1,1;0,2;-1,1;-2,3;0.2,3;-2,5;-3,1;-2,-4].',[0;0.5]);
        PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].',[0;0]);
       subFolderName = 'Examples for research proposal\cup'; % where to save results  
       f2 = [0.45;0.0];
       f1 = [0;0];
       basepos = [f1,f2];
       res = 100;
    case 2 % bottle
        load NewBottle.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
%         theta =  1.515;
%         R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
%         object = R*object;
        com(2) = max(object(2,:))-com(2);
        object(2,:) = max(object(2,:))-object(2,:);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [0.06;0.01];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 1500;
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
        f2 = [2.5;1];
        f1 = [-2.5;1];
        basepos = [f1,f2];
        res = 50;
end

%% Find equilibrium grasps

polyDrawing = PG.drawPolygon();

% hold on
% plot([-0.05 0.3],-0.091*[1 1],'k','lineWidth',3)
% xlim([-0.05 0.3])
% ylim([-0.2 0.3])
hold off
%%
tic
[PG,S,X,VL] = PG.findBdyVariable(res);
PG.res = res;
%%
[s1,s2,linesegs] = PG.Eqcheck();
toc
%% Check grasps for minimum/maximum
status = PG.Statuscheck(s1,s2,S,X,VL);
%% Full contact space figure
figure;
hold on
c_layers=50;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers)
xlabel('$s_1$','interpreter','latex','fontsize',18)
ylabel('$s_2$','interpreter','latex','fontsize',18)
colorbar
col = [0 1 0; 0 0 0; 1 0 0];
for i = 1:length(s1)
%     if status(i) ==-1
plot(s1(i),s2(i),'.','MarkerEdgeColor',col(status(i)+2,:),'MarkerFaceColor',col(status(i)+2,:),'markersize',8)
%     end
end
% for i=1:size(linesegs,1) % plotting the areas defined by horizontal edges
%     plot(linesegs(i,[1,3]),linesegs(i,[2,4]),'k','linewidth',2)
%     plot(linesegs(i,[1,3]),linesegs(i,[2,4]),'k.')
%     plot(linesegs(i,[2,4]),linesegs(i,[1,3]),'k','linewidth',2)
%     plot(linesegs(i,[2,4]),linesegs(i,[1,3]),'k.')
% end
% rectangle('Position',[10 0 7.83 S(end)],'FaceColor',[0,0,0,0.2])
%%
% lineseg = [];
% PG.drawPolygonS(s1,s2,S,X,lineseg);
% rectangle('Position',[0 0 0.3015 S(end)],'FaceColor',[0,0,0,0.2])
% rectangle('Position',[0 0 S(end) 0.3015],'FaceColor',[0,0,0,0.2])
