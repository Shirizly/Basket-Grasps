clearvars
close all
%% Define the polygon
% %% Define the polygon
clear all
clc
object = 0;
switch object
    case 0 % pentagon
        R = [0 1;-1 0];
        theta = 0;%-pi()/2;
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
         
             case 8 % modified cup
       
       PG = Polygon([-3,2;-2,-2;2,-2;2,-0.5;4,1.5;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results

%        basepos = [f1,f2];
       res = 50;
end

%

polyDrawing = PG.drawPolygon();
hold off
%% Find equilibrium grasps

[PG,S,X,VL] = PG.findBdyVariable(res);
[s1,s2,~] = PG.Eqcheck();

%% Check grasps for minimum/maximum
[BG,NBG] = PG.StatusSeparate(s1,s2,X);
%% Full contact space figure, choose the basket grasp curve
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
% Identify segments of basket grasps by pointing at end points:
BGSegEnds = [];
rectangle('Position',[0,-S(end)/5,S(end),S(end)/8],'FaceColor',[0 .5 .5])
e1 = ginput(1);
while e1(2)>=-0.0750*S(end)
    distfrome1 = (BG(1,:)-e1(1)).^2+(BG(2,:)-e1(2)).^2;
    inde1 = find(distfrome1==min(distfrome1),1,'first');
    BGSegEnds(:,end+1) = [BG(1,inde1);BG(2,inde1)]; %#ok<SAGROW>
    e1 = ginput(1);
end
end
plot(BGSegEnds(1,:),BGSegEnds(2,:),'m+','markerSize',10,'lineWidth',5)
hold off
%% draw the grasp segments on the object
nseg = size(BGSegEnds,2)-1;
col = [0,1/nseg,1-1/nseg];
[polyfig,dot] = drawPolygon2fGraspRange(PG,BGSegEnds(:,[1,2]),col);
dots = cell(0);
if ~isempty(dot)
    dots{end+1} = dot;
end
for i=2:nseg
    if mod(i,2)==0
        col = [0.8,1-i/nseg,i/nseg];
    else
        col = [0,i/nseg,1-i/nseg];
    end
    [~,dot] = drawPolygon2fGraspRange(PG,BGSegEnds(:,[i,i+1]),polyfig,col);
    if ~isempty(dot)
        dots{end+1} = dot;
    end
end
for i = 1:numel(dots)
    uistack(dots{i},'top')
end
%% Find basket grasp depth at discrete points along the segments:
BGdepthlist = []; sig = []; phi = []; BGslist = []; N = 30;absforce1 = [];pplist = [];
absforce2 = [];
for j=2:nseg+1
    BGs1j = linspace(BGSegEnds(1,j-1),BGSegEnds(1,j),N);
    BGs2j = linspace(BGSegEnds(2,j-1),BGSegEnds(2,j),N);
    BGslistj = [BGs1j;BGs2j];
    % BGslist(:,[1,N+2]) = [];
    BGdepthlistj = zeros(1,size(BGslistj,2));
    sigj = zeros(1,size(BGslistj,2));
    phij = zeros(1,size(BGslistj,2));
    pplistj = zeros(size(BGslistj,2),7);
    forces = [];
    for i=1:N
        BGS = BGslistj(:,i);
        %     [BGdepthlistj(i),sigj(i),phij(i),forces] = BasketDepthForces(PG,BGS);
        % %     check = forces(:,1)+forces(:,2);
        %     absforce1j(i) = norm(forces(:,1));
        %     absforce2j(i) = norm(forces(:,2));
        [BGdepthlistj(i),sigj(i),phij(i),pp] = BasketDepth(PG,BGS);
        pplistj(i,1:length(pp)) = pp;
    end
    BGdepthlist = [BGdepthlist,BGdepthlistj];
    sig = [sig,sigj];
    phi = [phi,phij];
    BGslist = [BGslist,BGslistj];
    % absforce1 = [absforce1,absforce1j];
    % absforce2 = [absforce2,absforce2j];
    pplist = [pplist;pplistj];
end
N = N*nseg;
%% Plot the discrete points along the S parameterization
% figure;hold on;
% plot(1:N,BGdepthlist);
% plot(1:N,sig);hold off;
figure;
hold on
plot(sig,BGdepthlist,'.k','markersize',10);
xlabel('$\sigma$','interpreter','latex','fontsize',25)
ylabel('$\Delta U$','interpreter','latex','fontsize',25)
title('Basket depth vs. inter-finger distance','interpreter','latex','fontsize',20)
hold off
% figure;
% hold on
% plot(sig,absforce1,'.b','markersize',10);
% plot(sig,absforce2,'.r','markersize',10);
% legend('magnitude1','magnitude2');
% xlabel('$\sigma$','interpreter','latex','fontsize',25)
% ylabel('$F$','interpreter','latex','fontsize',25)
% title('Reaction forces vs. inter-finger distance','interpreter','latex','fontsize',20)
% hold off







%% plot phi vs sigma of the basket grasp curve
% figure;
% hold on
% plot(sig,wrapToPi(phi+pi()),'.k','markersize',10);
% xlabel('$\sigma$','interpreter','latex','fontsize',25)
% ylabel('$\phi$','interpreter','latex','fontsize',25)
% title('Hand orientation vs. inter-finger distance','interpreter','latex','fontsize',20)
% % ylim([-pi()/4,pi()/4])
% xlim([min(sig),max(sig)])
% hold off
%%
% i = 15;
% BGS = BGslist(:,i);
% [BGdepthlist(i),sig(i)] = BasketDepthPath(PG,BGS);
% %% draw the caging regions on the object
% cagedpoly = PG.drawPolygon();
% hold on;
% finger1 = cell(1,N*nseg);
% finger2 = cell(1,N*nseg);
% for i=1:1:N*nseg-1
%     BGS = BGslist(:,i);
%     dep = 0.99*BGdepthlist(i);
%     if dep<1E-2
%         continue
%     end
%     sigi = sig(i);
%     red = 2*floor(i/N)/nseg - 2*floor(i/N)/nseg*(floor(i/N)>2);
%     col = [red,max(0,1-3*floor(i/N)/nseg),max(0,3*floor(i/N)/nseg-2)];
%     [finger1{i},finger2{i}] = drawCaging(PG,BGS,dep,col);
% end
% %% draw the caging regions in 3-D, with an inscribed circle within each
% figure
% hold on;
% for i=1:1:N*nseg-1
%     sigi = sig(i);
%     red = 2*floor(i/N)/nseg - 2*floor(i/N)/nseg*(floor(i/N)>2);
%     col = [red,max(0,1-3*floor(i/N)/nseg),max(0,3*floor(i/N)/nseg-2)];
%     fing1 = finger1{i};
%     fing2 = finger2{i};
%     if isempty(fing1)||isempty(fing2)
%         continue
%     end
%     plot3(fing1(1,:),fing1(2,:),sigi*ones(size(fing1(1,:))),'Color',col,'lineWidth',1)
%     plot3(fing2(1,:),fing2(2,:),sigi*ones(size(fing2(1,:))),'Color',col,'lineWidth',1)
%     x1 = fing1(1,:);
%     y1 = fing1(2,:);
%     [xCenter1,yCenter1,radius1] = InscribedDisc(x1,y1);
%     x2 = fing2(1,:);
%     y2 = fing2(2,:);
%     [xCenter2,yCenter2,radius2] = InscribedDisc(x2,y2);
%     tet = 0:0.01:2*pi();
%     plot3(xCenter1+radius1*cos(tet),yCenter1+radius1*sin(tet),sigi*ones(size(tet)),'Color',0.5*col,'lineWidth',3);
%     plot3(xCenter2+radius2*cos(tet),yCenter2+radius2*sin(tet),sigi*ones(size(tet)),'Color',0.5*col,'lineWidth',3);
% end

%%
% cagedpoly2 = PG.drawPolygon();
% hold on;
% for i=1:1:N*nseg
%     BGS = BGslist(:,i);
%     dep = BGdepthlist(i)-0.2;
%     sigi = sig(i);
%     col = 'g';
%     drawCaging(PG,BGS,dep,col);
% end
%%
function [xCenter,yCenter,radius] = InscribedDisc(x,y)
% Make data into a 1000x1000 image.
xMin = min(x);
xMax = max(x);
yMin = min(y);
yMax = max(y);
scalingFactor = 1000 / min([xMax-xMin, yMax-yMin]);
x2 = (x - xMin) * scalingFactor + 1;
y2 = (y - yMin) * scalingFactor + 1;
mask = poly2mask(x2, y2, ceil(max(y2)), ceil(max(x2)));
% Compute the Euclidean Distance Transform
edtImage = bwdist(~mask);
% Find the max
radius = max(edtImage(:));
% Find the center
[yCenter, xCenter] = find(edtImage == radius);
% Scale and shift the center back to the original coordinates.
xCenter = (xCenter - 1)/ scalingFactor + xMin;
yCenter = (yCenter - 1)/ scalingFactor + yMin;
radius = radius / scalingFactor;
% rectangle('Position',[xCenter-radius, yCenter-radius, 2*radius, 2*radius],'Curvature',[1,1]);
end

% source for function code:

% https://www.mathworks.com/matlabcentral/answers/377838-please-how-can-i-find-the-center-and-the-radius-of-the-inscribed-circle-of-a-set-of-points

% full code:

% % Initialization / clean-up code.
% clc;    % Clear the command window.
% close all;  % Close all figures (except those of imtool.)
% clear;  % Erase all existing variables. Or clearvars if you want.
% workspace;  % Make sure the workspace panel is showing.
% format long g;
% format compact;
% fontSize = 25;
% x=[2.6381318;-18.66985584;-39.97784349;-69.8090262;-99.64020891;-92.131167993;-84.62315095;-57.4301004;-30.23704914;-3.65279788;22.93145337;34.09278024;45.2541071;23.94611945]
% y=[75.28822303;63.72102973;52.15383644;43.02184173;33.88984702;19.48158871;5.07333039;5.07333039;5.07333039;5.07333039;5.07333039;25.77251893;46.4717064;60.87996471]
% subplot(2, 2, 1);
% plot(x, y, 'b.-', 'MarkerSize', 30);
% grid on;
% title('Original Points', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % Make data into a 1000x1000 image.
% xMin = min(x)
% xMax = max(x)
% yMin = min(y)
% yMax = max(y)
% scalingFactor = 1000 / min([xMax-xMin, yMax-yMin])
% x2 = (x - xMin) * scalingFactor + 1;
% y2 = (y - yMin) * scalingFactor + 1;
% mask = poly2mask(x2, y2, ceil(max(y2)), ceil(max(x2)));
% % Display the image.
% p2 = subplot(2, 2, 2);
% imshow(mask);
% axis(p2, 'on', 'xy');
% title('Mask Image', 'FontSize', fontSize);
% % Compute the Euclidean Distance Transform
% edtImage = bwdist(~mask);
% % Display the image.
% p3 = subplot(2, 2, 3);
% imshow(edtImage, []);
% axis(p3, 'on', 'xy');
% % Find the max
% radius = max(edtImage(:))
% % Find the center
% [yCenter, xCenter] = find(edtImage == radius)
% % Display circles over edt image.
% viscircles(p3, [xCenter, yCenter], radius);
% % Display polygon over image also.
% hold on;
% plot(x2, y2, 'r.-', 'MarkerSize', 30, 'LineWidth', 2);
% title('Euclidean Distance Transform with Circle on it', 'FontSize', fontSize);
% % Display the plot again.
% subplot(2, 2, 4);
% plot(x, y, 'b.-', 'MarkerSize', 30);
% grid on;
% % Show the circle on it.
% hold on;
% % Scale and shift the center back to the original coordinates.
% xCenter = (xCenter - 1)/ scalingFactor + xMin
% yCenter = (yCenter - 1)/ scalingFactor + yMin
% radius = radius / scalingFactor
% rectangle('Position',[xCenter-radius, yCenter-radius, 2*radius, 2*radius],'Curvature',[1,1]);
% title('Original Points with Inscribed Circle', 'FontSize', fontSize);