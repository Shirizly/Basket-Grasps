clearvars
close all
%% Define the polygon
% %% Define the polygon
clc
object = 4;
switch object
    case 0 % pentagon
        theta = 0;%-pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        res = 50;
        
    case 1 % cup
       
       PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;0],[0;0]);
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
       
    case 9
        load DrawnObject.mat
        object = Round_object;
        PG = Polygon(object,[0;0]);
        res = 50;
end
[PG,S,X,VL] = PG.findBdyVariable(res);
%%

polyDrawing = PG.drawPolygon();
hold off

%% find the line segment of basket grasps curve (for the positions of the fingers)
% Find equilibrium grasps


[s1,s2,~] = PG.Eqcheck();

% Check grasps for minimum/maximum
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
%
drawPolygon2fGraspRange(PG,BGSegEnds)
%% Identify the fingerline vector over the segment as a function of S
[p11,p12] = PG.get('2Pos',BGSegEnds(1,1),BGSegEnds(1,2));
[p21,p22] = PG.get('2Pos',BGSegEnds(2,1),BGSegEnds(2,2));
FLa0 = p21-p11;
dels1 = BGSegEnds(1,2)-BGSegEnds(1,1); % s difference for finger 1 (over the segment)
dels2 = BGSegEnds(2,2)-BGSegEnds(2,1); % s difference for finger 2
FLa1 = ((p22-p21)-(p12-p11))/dels1;
%dp(S) = FLa0+S*FLa1; This is the function of the fingerline vector
Sstart = BGSegEnds(1,1);
Send = BGSegEnds(1,2);
%% Generate the solutions of edge-edge equilibrium for the range of S
curves = cell(0);
for edge1 = 1:PG.nv
    for edge2 = [1:edge1-1 edge1+1:PG.nv]
%         edge1 = 3;
%         edge2 = 6;
        dels1 = PG.S(PG.VL(edge1+1))-PG.S(PG.VL(edge1));
        dels2 = PG.S(PG.VL(edge2+1))-PG.S(PG.VL(edge2));
        
        nj = PG.normal(:,edge1);
        nj = nj/norm(nj);
        nk = PG.normal(:,edge2);
        nk = nk/norm(nk);
        epsi = 1E-7;
        J = [0 -1;1 0];
        if abs(nj.'*J*nk)<epsi
            continue
        end
        vj = PG.vertex(:,edge1);
        vk = PG.vertex(:,edge2);
        %
        N = 500;
        Sv = linspace(Sstart,Send,N);
        sol1 = zeros(3,size(Sv,2));
        sol2 = zeros(3,size(Sv,2));
        sol3 = zeros(3,size(Sv,2));
        sol4 = zeros(3,size(Sv,2));
        dp = zeros(2,size(Sv,2));
        
        for i=1:N
            dp(:,i) = FLa0+(Sv(i)-Sstart)*FLa1;
            [sol1(:,i),sol2(:,i),sol3(:,i),sol4(:,i)] = EEheight(dp(:,i),nj,nk,vj,vk,dels1,dels2);
        end
        %
        if any([~isnan(sol1(3,:)),~isnan(sol2(3,:)),~isnan(sol3(3,:)),~isnan(sol4(3,:))])
            curve = [];
            figure
            for i=1:N
                hold on
                if any([~isnan(sol1(3,i)),~isnan(sol2(3,i)),~isnan(sol3(3,i)),~isnan(sol4(3,i))])
                    plot(0,0,'og','linewidth',4)
                    plot(dp(1,i),dp(2,i),'og','linewidth',4)
                end
                if ~isnan(sol1(3,i))
                    PG.drawPolygonMoved(PG.com,sol1(3,i),sol1(1:2,i),'k');
                    curve = [curve,[sol1(2:3,i);Sv(i)]];
                end
                if ~isnan(sol2(3,i))
                    PG.drawPolygonMoved(PG.com,sol2(3,i),sol2(1:2,i),'k');
                    curve = [curve,[sol2(2:3,i);Sv(i)]];
                end
                if ~isnan(sol3(3,i))
                    PG.drawPolygonMoved(PG.com,sol3(3,i),sol3(1:2,i),'k');
                    curve = [curve,[sol3(2:3,i);Sv(i)]];
                end
                if ~isnan(sol4(3,i))
                    PG.drawPolygonMoved(PG.com,sol4(3,i),sol4(1:2,i),'k');
                    curve = [curve,[sol4(2:3,i);Sv(i)]];
                end
                
                hold off
            end
            curves{end+1} = curve;
        end
    end
end
%%
figure
hold on
for i=1:numel(curves)
    curve = curves{i};
    plot(curve(3,:),curve(1,:),'.')
end
title('Edge-edge escape height')
xlabel('S - parameter')
ylabel('Height')
hold off
%%
% tic
curves = cell(0);
g = 0;
for edge1 = 1:PG.nv
    for edge2 = [1:edge1-1 edge1+1:PG.nv]
        %         edge1 = 3;
        %         edge2 = 6;
        dels1 = PG.S(PG.VL(edge1+1))-PG.S(PG.VL(edge1));
        dels2 = PG.S(PG.VL(edge2+1))-PG.S(PG.VL(edge2));
        nj = PG.normal(:,edge1);
        nj = nj/norm(nj);
        nk = PG.normal(:,edge2);
        nk = nk/norm(nk);
        epsi = 1E-7;
        J = [0 -1;1 0];
        if abs(nj.'*J*nk)<epsi
            continue
        end
        vj = PG.vertex(:,edge1);
        vk = PG.vertex(:,edge2);
        %
        N = 500;
        Sv = linspace(Sstart,Send,N);
        sol1 = zeros(3,size(Sv,2));
        sol2 = zeros(3,size(Sv,2));
        sol3 = zeros(3,size(Sv,2));
        sol4 = zeros(3,size(Sv,2));
        dp = zeros(2,size(Sv,2));
        tic;
        for i=1:N
            dp(:,i) = FLa0+(Sv(i)-Sstart)*FLa1;
            [sol1(:,i),sol2(:,i),sol3(:,i),sol4(:,i)] = EEheight(dp(:,i),nj,nk,vj,vk,dels1,dels2);
        end
        toc
        %
        %         if any([~isnan(sol1(3,:)),~isnan(sol2(3,:)),~isnan(sol3(3,:)),~isnan(sol4(3,:))])
        %             curve = [];
        %             for i=1:N
        %                 if ~isnan(sol1(3,i))
        %                     curve = [curve,[sol1(2:3,i);Sv(i)]];
        %                 end
        %                 if ~isnan(sol2(3,i))
        %                     curve = [curve,[sol2(2:3,i);Sv(i)]];
        %                 end
        %                 if ~isnan(sol3(3,i))
        %                     curve = [curve,[sol3(2:3,i);Sv(i)]];
        %                 end
        %                 if ~isnan(sol4(3,i))
        %                     curve = [curve,[sol4(2:3,i);Sv(i)]];
        %                 end
        %             end
        %             curves{end+1} = curve;
        %         end
        
        if any(~isnan(sol1(3,:)))
            curveSt = find(~isnan(sol1(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol1(3,curveSt+1:end)),1,'first')-1;
            curve = sol1(:,curveSt:curveEn);
            if ~isempty(find(~isnan(sol1(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
            curves{end+1} = curve;
        end
        if any(~isnan(sol2(3,:)))
            curveSt = find(~isnan(sol2(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol2(3,curveSt+1:end)),1,'first')-1;
            curve = sol2(:,curveSt:curveEn);
            if ~isempty(find(~isnan(sol2(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
            curves{end+1} = curve;
        end
        if any(~isnan(sol3(3,:)))
            curveSt = find(~isnan(sol3(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol3(3,curveSt+1:end)),1,'first')-1;
            curve = sol3(:,curveSt:curveEn);
            if ~isempty(find(~isnan(sol3(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
            curves{end+1} = curve;
        end
        if any(~isnan(sol4(3,:)))
            curveSt = find(~isnan(sol4(3,:)),1,'first');
            curveEn = curveSt + find(isnan(sol4(3,curveSt+1:end)),1,'first')-1;
            curve = sol4(:,curveSt:curveEn);
            if ~isempty(find(~isnan(sol4(3,curveEn+1:end)),1))
                disp('Discontinuity!')
            end
            curves{end+1} = curve;
        end
        
    end
end
% toc