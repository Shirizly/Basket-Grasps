clearvars
close all
%% Define the polygon
% %% Define the polygon
clear all
clc
object = 1;
switch object
    case 0 % pentagon
        R = [0 1;-1 0];
        theta = -pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        res = 50;
        
    case 1 % cup
       
       PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results

%        basepos = [f1,f2];
       res = 50;
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
        res = 50;
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
        
        case 7 % Crown with more edges
        object = [-2,-1;2,-1;2.25,1;3.25,3;1,0;0,1;-1,0;-4,1;-2.75,0].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 50;
         starting_node = [2,4];  
end

%

polyDrawing = PG.drawPolygon();
hold off
%% Find equilibrium grasps

[PG,S,X,VL] = PG.findBdyVariable(res);
[s1,s2,~] = PG.Eqcheck();

%% Check grasps for minimum/maximum
[BG,NBG] = PG.StatusSeparate(s1,s2,X);
%% Full contact space figure
contspace = figure;
hold on
c_layers=20;
Sigma=inter_finger_distance(X,X);
contour(S,S,Sigma,c_layers)
xlabel('$s_1$','interpreter','latex','fontsize',18)
ylabel('$s_2$','interpreter','latex','fontsize',18)
colorbar
col = [0 1 0; 0 0 0; 1 0 0];
plot(BG(1,:),BG(2,:),'.','MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'markersize',8)
plot(NBG(1,:),NBG(2,:),'.','MarkerEdgeColor',col(3,:),'MarkerFaceColor',col(3,:),'markersize',8)
if ~exist('BGSegEnds','var')
% Identify segment of basket grasps by pointing at end points:
e1 = ginput(1);
e2 = ginput(1);
distfrome1 = (BG(1,:)-e1(1)).^2+(BG(2,:)-e1(2)).^2;
inde1 = find(distfrome1==min(distfrome1),1,'first');
distfrome2 = (BG(1,:)-e2(1)).^2+(BG(2,:)-e2(2)).^2;
inde2 = find(distfrome2==min(distfrome2),1,'first');
BGSegEnds = [BG(1,inde1),BG(1,inde2);BG(2,inde1),BG(2,inde2)];
end
plot(BGSegEnds(1,:),BGSegEnds(2,:),'m+','markerSize',14,'lineWidth',5)
hold off
%%
drawPolygon2fGraspRange(PG,BGSegEnds)

%% Find basket grasp depth at discrete points along the segments:

N = 20;
BGs1 = linspace(BGSegEnds(1,1),BGSegEnds(1,2),N);
BGs2 = linspace(BGSegEnds(2,1),BGSegEnds(2,2),N);
BGslist = [BGs1;BGs2];


i = 15;
BGS = BGslist(:,i);
draw2fGrasp(PG,BGS(1),BGS(2),0)
%%


[~,sigma,iphi] = BasketDepth(PG,BGS);
N=50;

phirange = linspace(iphi-pi()/4,iphi+pi()/4,N);
BGdepthlistj = zeros(size(phirange));
for i=1:N
    phi = phirange(i);
    [BGdepthlistj(i),sigj(i),phij(i)] = BasketDepthDifTheta(PG,sigma,phi,BGS);
end
%%
figure
plot(phirange-iphi,BGdepthlistj)
xlabel('$\Delta \phi$','interpreter','latex','fontsize',25)
ylabel('$\Delta U$','interpreter','latex','fontsize',25)
title('Basket depth vs. inter-finger angle','interpreter','latex','fontsize',20)