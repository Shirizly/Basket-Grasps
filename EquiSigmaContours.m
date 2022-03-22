clearvars
close all
%% Define the polygon
% %% Define the polygon
clc
object = 4;
symmetry = 0;
switch object
    case 0 % pentagon
        theta = 0;%-pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        res = 50;
        symmetry = [1,2,3,3,2;1,5,4,3,2;-1,-1,-5,-4,-5;-2,-1,-5,-4,-3];
    case 1 % cup
       PG = Polygon([-2,-2;2,-2;2,1;4,1;4,2;-3,2].'+[0;3],[0;3]);
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
        object = Round_object;
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
    case 8
                load new.mat
        object = Round_object(:,1:end-1);
        PG = Polygon(object,com);
        res = 50;
end

polyDrawing = PG.drawPolygon();
hold off
%% Find equilibrium grasps
tic
[PG,S,X,VL] = PG.findBdyVariable(res);
[s1,s2,~] = PG.Eqcheck();
toc
%% Check grasps for minimum/maximum
[BG,NBG] = PG.StatusSeparate(s1,s2,X);
%% Full contact space figure
contspace = figure;
hold on
for i=1:PG.nv
    plot([0 PG.S(end)],PG.S(PG.VL(i))*ones(1,2),'k--','linewidth',1)
    plot(PG.S(PG.VL(i))*ones(1,2),[0 PG.S(end)],'k--','linewidth',1)
end
c_layers=linspace(6.25,6.35,400);
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
% drawPolygon2fGraspRange(PG,BGSegEnds)

%%
% figure
clc
N = 50;
BGs1 = linspace(BGSegEnds(1,1),BGSegEnds(1,2),N);
BGs2 = linspace(BGSegEnds(2,1),BGSegEnds(2,2),N);
BGslist = [BGs1;BGs2];
tic
EquiSigmaNumericCell = cell(0);
pathpointsCell = cell(0);
for i=1:N
BGS = BGslist(:,i);
syms t
s1init = BGS(1);
s2init = BGS(2);
[O1,O2] = PG.get('2Pos',s1init,s2init);
sigma = norm(O2-O1);
relvector = O2-O1;
Fullcircuit = 0;
EquiSigmaContour = cell(0);
EquiSigmaNumeric = [];
s10 = s1init;
s20 = s2init;
s1prev = s10; s2prev = s20;
theta0 = 0;
R = [cos(t) sin(t);-sin(t) cos(t)];
J = [0 -1;1 0];
thetavector = 0:-0.001:-pi();
thetadir = 1;
prevedge = [0;0];
startIn = 1;
path = [];
pathpoints = [];
thetapoints = [];
s1vertex = [];
s2vertex = [];
s1vertexalt = [];
s2vertexalt = [];

while ~Fullcircuit %construct the contour from segments in each rectangle
% the contour passes through
try
    % Identify the edge-edge rectangle being explored
    Rt0 = [cos(theta0) -sin(theta0);sin(theta0) cos(theta0)];
    noRotation = 0;
    edge1 = PG.get('edgeNum',s10);
    edge2 = PG.get('edgeNum',s20);
    if all([edge1;edge2]==prevedge) % check if the default rectangle is the one already previously visited by the loop
        %fix the direction of the contour by picking the neighboring
        %rectangle - setting the edge supported by the vert
        isvertex = any(PG.S(PG.VL)==s10)*1 + any(PG.S(PG.VL)==s20)*2;
        switch isvertex
            case 1
                edge1 = edge1-1;
                if edge1 == 0
                    edge1 = PG.nv;
                    s10 = PG.S(end);
                end
            case 2
                edge2 = edge2-1;
                if edge2 == 0
                    edge2 = PG.nv;
                    s20 = PG.S(end);
                end
            case 3 
                disp('vertex-vertex case')
                edge1 = edge1-1; %to prevent crossing a corner, we move s2 a tiny bit
                if edge1 == 0
                    edge1 = PG.nv;
                    s10 = PG.S(end);
                end
                if s20>0
                    s20 = s20-1E-8;
                else
                    s20 = s20+1E-8;
                end
        end
    end 
    if edge1 == 1 && s10 == PG.S(end)
        s10 = 0;
    end
    if edge2 == 1 && s20 == PG.S(end)
        s20 = 0;
    end   
    % find the functions describing the contour
    s1left = PG.S(PG.VL(edge1+1))-s10;
    s1leftalt = PG.S(PG.VL(edge1))-s10;
    s2left = PG.S(PG.VL(edge2+1))-s20;
    s2leftalt = PG.S(PG.VL(edge2))-s20;
    if any([s1left,s2left,s1leftalt,s2leftalt]==0)
        startIn = 0;
    end  
    n1 = Rt0*PG.normal(:,edge1);
    n2 = Rt0*PG.normal(:,edge2);
    if abs(det([J*n1 -J*n2]))>1E-10
        dsvector = [J*n1 -J*n2]\(R-eye(2))*relvector;
        Sv = thetavector;
        if abs(det([J*n1 -J*n2]))<0.5 %edges are close to parallel, small rotations cause high displacements
            % to have a high-resolution curve, add resolution to theta
            % accordingly
            Sv = thetavector(1):diff(thetavector(1:2))*10^(round(log10(abs(det([J*n1 -J*n2]))),0)):thetavector(end);
        end
        % find the end of the segment in the current rectangle
        dsf = matlabFunction(dsvector);
        dsmatrix = feval(dsf,Sv);
        if startIn %special case when starting inside a rectangle
            s1vertex = find(dsmatrix(1,1:end)>=s1left,1,'first');
            s1vertexalt = find(dsmatrix(1,1:end)<=s1leftalt,1,'first');
            s2vertex = find(dsmatrix(2,1:end)>=s2left,1,'first');
            s2vertexalt = find(dsmatrix(2,1:end)<=s2leftalt,1,'first');
        else
            [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,0);
            if all(isempty([s1vertex,s2vertex,s1vertexalt,s2vertexalt])) % need to rotate in the opposite direction
                Sv = -Sv;
                dsmatrix = feval(dsf,Sv);
                [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,0);
            end
        end
        if isempty(s1vertex)
            s1vertex = inf;
        end
        if isempty(s2vertex)
            s2vertex = inf;
        end
        if isempty(s1vertexalt)
            s1vertexalt = inf;
        end
        if isempty(s2vertexalt)
            s2vertexalt = inf;
        end
        s1vertex = min([s1vertex,s1vertexalt]);
        s2vertex = min([s2vertex,s2vertexalt]);
        if min(s1vertex,s2vertex)==inf
            disp('intercept error');
        end
        Sv = Sv(1:min([s1vertex,s2vertex]));
        n1 = PG.normal(:,edge1);
        n2 = PG.normal(:,edge2);

        while s1vertex==s2vertex && abs(Sv(end)-Sv(end-1))>1E-8 %curve segments exits rectangle too close to corner to be clear which edge is crossed
            %add checking which edge is crossed
            Svtag = linspace(Sv(end-1),Sv(end),10);
            dsmatrixtag = feval(dsf,Svtag);
            [s1vertextag,s1vertexalttag,s2vertextag,s2vertexalttag] = checkRectangle(dsmatrixtag,s1left,s1leftalt,s2left,s2leftalt,1);
            if isempty(s1vertextag)
                s1vertextag = inf;
            end
            if isempty(s2vertextag)
                s2vertextag = inf;
            end
            if isempty(s1vertexalttag)
                s1vertexalttag = inf;
            end
            if isempty(s2vertexalttag)
                s2vertexalttag = inf;
            end
            mins1v = min([s1vertextag,s1vertexalttag]);
            mins2v = min([s2vertextag,s2vertexalttag]);
            minsv = min(mins1v,mins2v);
            s1vertex = s1vertex-2+mins1v;
            s2vertex = s2vertex-2+mins2v;
            Sv = [Sv(1:end-2),Svtag(1:minsv)];
        end
        if s1vertex==s2vertex
            disp('corner!');
        end
        dsmatrix = feval(dsf,Sv);
        s1prev = s10;
        s2prev = s20;
        if s1vertex<s2vertex
            s10 = PG.S(PG.VL(edge1+(dsmatrix(1,min([s1vertex,s2vertex]))>0)));
            s21 = s20+dsmatrix(2,end);
            if s21>PG.S(PG.VL(edge2+1))
                s21 = PG.S(PG.VL(edge2+1))
            end
            if s21<PG.S(PG.VL(edge2))
                s21 = PG.S(PG.VL(edge2))
            end
            [O10,O20] = PG.get('2pos',s10,s21);
            sigmagal = norm(O20-O10);
            p = [n2 -J*n2]\(O10-O20);
            h = p(1);
            s20check = s21+sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2);
            s20 = s21+sign(p(2))*(sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2));
            O1check = O10;
            O2check = PG.get('1pos',s20);
            fingers= O2check-O1check;
            sigmacheck = norm(fingers)-sigma;
%             drawsigmaissue(PG,n2,sigma,s10,s21)
%             set(0, 'CurrentFigure', contspace);
            theta0prev = theta0;
            thetatag = atan2(relvector(2),relvector(1));
            theta0 = atan2(relvector(2),relvector(1))-atan2(fingers(2),fingers(1));
        else
            s20 = PG.S(PG.VL(edge2+(dsmatrix(2,min([s1vertex,s2vertex]))>0)));
            s11 = s10+dsmatrix(1,end);
            if s11>PG.S(PG.VL(edge1+1))
                s11 = PG.S(PG.VL(edge1+1))
            end
            if s11<PG.S(PG.VL(edge1))
                s11 = PG.S(PG.VL(edge1))
            end
            [O10,O20] = PG.get('2pos',s11,s20);
            sigmagal = norm(O20-O10);
            p = [n1 -J*n1]\(O20-O10);
            h = p(1);
            s10 = s11+sign(p(2))*(sqrt(sigmagal^2-h^2)-sqrt(sigma^2-h^2));
            O1check = PG.get('1pos',s10);
            fingers= O20-O1check;
            sigmacheck = norm(fingers)-sigma;
%             drawsigmaissue(PG,n1,sigma,s20,s11)
%             set(0, 'CurrentFigure', contspace);
%             if abs(sigmacheck)>abs(sigmagal-sigma)
%                 disp('error start');
%             end
            theta0prev = theta0;
            thetatag = atan2(relvector(2),relvector(1));
            theta0 = atan2(relvector(2),relvector(1))-atan2(fingers(2),fingers(1));
        end
        dsmatrix = [s1prev;s2prev]+dsmatrix;
    
    
    
    else % the supported edges are parallel
        dsvector = t*[1;sign(n1.'*n2)]; % the object isn't rotating, instead use ds1 for the parametrization of the segment
        if sign(n1.'*n2)>0
            Sv = [0:0.001:min([s1left,s2left]),min([s1left,s2left])];
        else
            Sv = [0:0.001:min([s1left,-s2leftalt]),min([s1left,-s2leftalt])];
        end
        if s1left==0
            if sign(n1.'*n2)>0
                Sv = [0:-0.001:max(s1leftalt,s2leftalt),max(s1leftalt,s2leftalt)];
            else
                Sv = [0:-0.001:max(s1leftalt,-s2left),max(s1leftalt,-s2left)];
            end
        end
        if s2left==0
            if sign(n1.'*n2)>0
                Sv = [0:-0.001:max(s1leftalt,s2leftalt),max(s1leftalt,s2leftalt)];
            else
                Sv = [0:0.001:min(s1left,-s2leftalt),min(s1left,-s2leftalt)];
            end
        end
        if s2leftalt==0
            if sign(n1.'*n2)>0
                Sv = [0:0.001:min(s1left,s2left),min(s1left,s2left)];
            else
                Sv = [0:-0.001:max(s1leftalt,-s2left),max(s1leftalt,-s2left)];
            end
        end
        if Sv(2)<0
            Sv = -Sv;
            dsvector = -dsvector;
        end
        dsf = matlabFunction(dsvector);
        dsmatrix = feval(dsf,Sv);
        noRotation = 1;
        s1prev = s10;
        s2prev = s20;
        s10 = s10+dsmatrix(1,end);
        s20 = s20+dsmatrix(2,end);
        dsmatrix = dsmatrix + [s1prev;s2prev];
    end
    %     if theta0>2*pi()||theta0<-2*pi()
    %         Fullcircuit = 1;
    %     end
    if ~isempty(path)
        if find(all(path==[edge1;edge2],1))
            if find(all(abs(pathpoints-[s10;s20])<1E-4),1)
                Fullcircuit = 1;
            end
        end
    end
    path = [path,[edge1;edge2]];
    pathpoints = [pathpoints,[s10;s20]];
    thetapoints = [thetapoints,[theta0]];
    
    if startIn
        startIn=0;
    end
    % add the resulting segment to the numerical and analytical
    % representations of the double-support contour
    if ~noRotation
        EquiSigmaNumeric = [EquiSigmaNumeric,[dsmatrix;theta0prev+Sv]];
    else
        EquiSigmaNumeric = [EquiSigmaNumeric,[dsmatrix;theta0*ones(size(Sv))]];
    end
    
    % prepare the next iteration of the loop
    if imag(s10)~=0
        disp('big error')
    end
    hold on
    plot(dsmatrix(1,:),dsmatrix(2,:),'k','linewidth',2)
%     if size(dsmatrix,2)<5
%          disp('weird length error');
%     end
    prevedge = [edge1;edge2];
    s1vertex = [];
    s2vertex = [];
    s1vertexalt = [];
    s2vertexalt = [];
catch errmsg
    errmsg.message
    errmsg.stack.line
    disp(errmsg)
end
end
pathpointsCell{end+1} = [pathpoints;thetapoints];
EquiSigmaNumericCell{end+1} = EquiSigmaNumeric;
end
toc
%%
s1st = BGSegEnds(1,2);
s2st = BGSegEnds(2,2);
[f1,f2] = PG.get('2pos',s1st,s2st);
basepos = [f1,f2];
sig = norm(f2-f1);
[cont_original] = PG.GetSigmaContours(Sigma,sig);


cont = PG.CleanContour(cont_original,basepos);

[ds_max,ds_min,ds_virtual] = PG.DSNodes(cont);
[ss_max,ss_min,ss_saddle] = PG.SSNodes();
nodes{1} = ds_max;
nodes{2} = ds_min;
nodes{3} = ds_virtual;
nodes{4} = ss_max;
nodes{5} = ss_min;
nodes{6} = ss_saddle;

nodes = check_SS_for_DS(PG,nodes,f1,f2); %function checks SS_nodes for penetration (of the other finger),
% each node appears once for each finger it is relevant for, with the relevant finger index at the end 
% save([pwd '/' subFolderName '/PreGraph.mat'],'PG','nodes','cont','starting_node','f1','f2');
%% draw the contact space
loadfromfile = 0;
if loadfromfile == 1
    for finger = 1:2
        if exist([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig'],'file')
            graphfig{finger} = openfig([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig']); 
%             graphfig{finger} = openfig([pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.fig']); 
        end
    end
end
if ~exist('graphfig','var')
    graphfig = cell(1,2);
end

for finger = 1:2
    if isempty(graphfig{finger})||~isgraphics(graphfig{finger})
            otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
            contactHeight = basepos(2,finger);
            graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight);
            hold on
            for i=1:numel(EquiSigmaNumericCell)
                EquiSigmaNumeric = EquiSigmaNumericCell{i};
                plot(pathpointsNumeric(finger,:),pathpointsNumeric(3,:),'.r','lineWidth',18,'markerSize',14);
                plot(EquiSigmaNumeric(finger,:),EquiSigmaNumeric(3,:),'.k');
            end
            for i=[1,numel(EquiSigmaNumericCell)]
                EquiSigmaNumeric = EquiSigmaNumericCell{i};
                plot(EquiSigmaNumeric(finger,:),EquiSigmaNumeric(3,:),'.b');
            end
    end
    set(graphfig{finger},'Position',[(finger-1)*960,41,960,963])
    set(gcf,'renderer','Painters')
end
%%
for finger = 1:2


otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
contactHeight = basepos(2,finger);
graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight)
hold on
sigma0 = [];
for i=1:numel(EquiSigmaNumericCell)
    EquiSigmaNumeric = EquiSigmaNumericCell{i};
    pathpointsNumeric = pathpointsCell{i};
    s1init = EquiSigmaNumeric(1,1);
    s2init = EquiSigmaNumeric(2,1);
    [O1,O2] = PG.get('2Pos',s1init,s2init);
    sigma = norm(O2-O1);
    if isempty(sigma0)
        sigma0 = sigma;
    end
    plot3(EquiSigmaNumeric(finger,:),EquiSigmaNumeric(3,:),abs(sigma-sigma0)*ones(size(EquiSigmaNumeric(3,:))),'.b','lineWidth',3);
    plot3(pathpointsNumeric(finger,:),pathpointsNumeric(3,:),abs(sigma-sigma0)*ones(size(pathpointsNumeric(3,:))),'.r','lineWidth',12,'markerSize',14);
end
hold off
end
%%
function [s1vertex,s1vertexalt,s2vertex,s2vertexalt] = checkRectangle(dsmatrix,s1left,s1leftalt,s2left,s2leftalt,startIn)
s1vertex = [];
s2vertex = [];
s1vertexalt = [];
s2vertexalt = [];
% find on which of the rectangle's edge does the curve begin:
edge = 0; %edge number represents the edge in right-top-left-bottom order
if s1left ==0
    edge = edge+1;
end
if s2left ==0
    edge = edge+2;
end
if s1leftalt ==0
    edge = edge+3;
end
if s2leftalt ==0
    edge = edge+4;
end
% check if the curve found enters the right rectangle:
if ~startIn % startIn represents whether the curve segment starts inside the rectangle or on the boundary
    tangent = diff(dsmatrix(:,1:2),1,2); % if started on the boundary, check the direction of the curve to make sure it enters the right rectangle
    switch edge
        case 1
            if tangent(1)>0
                return
            end
        case 2
            if tangent(2)>0
                return
            end
        case 3
            if tangent(1)<0
                return
            end
        case 4
            if tangent(2)<0
                return
            end
    end
end
% check crossing the boundary only if it is not the starting boundary
% if s1left~=0
    s1vertex = find(dsmatrix(1,1:end)>s1left,1,'first');
%     if ~isempty(s1vertex)
%         ds2 = dsmatrix(2,s1vertex);
%         if ds2>s2left ||ds2<s2leftalt
%             s1vertex=[];
%         end
%     end
% end
% if s1leftalt~=0
    s1vertexalt = find(dsmatrix(1,1:end)<s1leftalt,1,'first');
%     if ~isempty(s1vertexalt)
%         ds2 = dsmatrix(2,s1vertexalt);
%         if ds2>s2left ||ds2<s2leftalt
%             s1vertexalt=[];
%         end
%     end
% end
% if s2left~=0
    s2vertex = find(dsmatrix(2,1:end)>s2left,1,'first');
%     if ~isempty(s2vertex)
%         ds1 = dsmatrix(1,s2vertex);
%         if ds1>s1left ||ds1<s1leftalt
%             s2vertex=[];
%         end
%     end
% end
% if s2leftalt~=0
    s2vertexalt = find(dsmatrix(2,1:end)<s2leftalt,1,'first');
%     if ~isempty(s2vertexalt)
%         ds1 = dsmatrix(1,s2vertexalt);
%         if ds1>s1left ||ds1<s1leftalt
%             s2vertexalt=[];
%         end
%     end
% end
end


%%
% if path(:,end)==path(:,2) %making the circuit complete with no overlap
%     path(:,end)=path(:,1);
%     
% end
    
% function ds = getOtherFingerSigma(PG,sigma,relfinger,dsj)
% 
% end


