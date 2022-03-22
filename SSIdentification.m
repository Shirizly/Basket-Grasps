%% Define the polygon
clearvars
clc
type = 4;
switch type
    case 0 % pentagon
        R = [0 1;-1 0];
        theta = -pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[0;3]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
%         f2 = [-4.4;1.8];
%         f1 = [4.4;1.8];
f2 = [8.3223;-1.3181];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,4];
        res = 50;
%         starting_node = [2,3];
    case 1 % cup
        % PG = Polygon([1,0;2,0;2,3;-1,4;-2,3;-1,1;-0.5,1].',[0;1.5]);
        % PG = Polygon([-1,-1;2,-1;3,-3;3,3;1,1;0,2;-1,1;-2,3;0.2,3;-2,5;-3,1;-2,-4].',[0;0.5]);
        PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results
       f2 = [2;4];
       f1 = [-2.5;2.8];
       basepos = [f1,f2];
       res = 50;
       starting_node = [2,6];
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

% draw the polygon (only once)
% if ~exist('polyDrawing','var')
%     polyDrawing = cell(0);
% end
% if isempty(polyDrawing)||~isgraphics(polyDrawing)
%     polyDrawing = PG.drawPolygon();
% end

% Find equilibrium grasps

% [PG,S,X,THL,CVL,~] = PG.findBdyVariable3(res);
[PG,S,X] = PG.findBdyVariable(res);

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
%
drawPolygon2fGraspRange(PG,BGSegEnds)

%% identify nodes, create double-support contours
Sigma=inter_finger_distance(X,X);

f1 = [0;0];
f2 = PG.get('1Pos',BGSegEnds(2,2))-PG.get('1Pos',BGSegEnds(1,2));
f2tag = PG.get('1Pos',BGSegEnds(2,2))-PG.get('1Pos',BGSegEnds(1,2));

sig = norm(f2-f1);
[cont_original] = PG.GetSigmaContours(Sigma,sig);
basepos = [f1,f2];
cont = PG.CleanContour(cont_original,basepos);

sigtag = norm(f2tag-f1);
[cont_original] = PG.GetSigmaContours(Sigma,sigtag);
basepostag = [f1,f2tag];
conttag = PG.CleanContour(cont_original,basepostag);

[ds_max,ds_min,ds_virtual] = PG.DSNodes(cont);
[ds_maxtag,ds_mintag,ds_virtualtag] = PG.DSNodes(conttag);
[ss_max,ss_min,ss_saddle] = PG.SSNodes();
nodes{1} = [ds_max;ds_maxtag];
nodes{2} = [ds_min;ds_mintag];
nodes{3} = [ds_virtual;ds_virtualtag];
nodes{4} = ss_max;
nodes{5} = ss_min;
nodes{6} = ss_saddle;

nodes = check_SS_for_DS(PG,nodes,f1,f2); %function checks SS_nodes for penetration (of the other finger),
% each node appears once for each finger it is relevant for, with the relevant finger index at the end 
% save([pwd '/' subFolderName '/PreGraph.mat'],'PG','nodes','cont','starting_node','f1','f2');
%% draw the contact space
loadfromfile = 1;
if loadfromfile == 1
    for finger = 1:2
        if exist([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig'],'file')
            graphfig{finger} = openfig([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig']); 
%             graphfig{finger} = openfig([pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.fig']); 
            hold on
            plot(EquiSigmaNumeric(finger,:),EquiSigmaNumeric(3,:),'.k')
            hold off
        end
    end
end
if ~exist('graphfig','var')
    graphfig = cell(1,2);
end

for finger = 1:2
    if isempty(graphfig{finger})||~isgraphics(graphfig{finger})
            otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
            otherFingerRelativeTag = basepostag(:,3-finger)-basepostag(:,finger);
            contactHeight = basepos(2,finger);
%             graphfig{finger} = plotGraphNodesTag(PG,cont,conttag,finger,nodes,otherFingerRelative,otherFingerRelativeTag,contactHeight);
            
            graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight);
            hold on
            for i=1:numel(EquiSigmaNumericCell)
                EquiSigmaNumeric = EquiSigmaNumericCell{i};
            plot(EquiSigmaNumeric(finger,:),EquiSigmaNumeric(3,:),'.k');
            end
            hold off
    end
     set(graphfig{finger},'Position',[(finger-1)*960,41,960,920])
    set(gcf,'renderer','Painters')
end
%%
% tic
% [OpenList,ClosedList,A,index_list,type_list,path_list,FingerMatrix,pathwayMatrix] = GraphConstruction(PG,cont,nodes,f1,f2,graphfig);
% toc
% save([pwd '/' subFolderName '/completeGraph.mat'],'ClosedList','A','index_list','type_list','path_list','FingerMatrix','pathwayMatrix');
%% Graphical presentation of the results:
% load([pwd '/' subFolderName '/completeGraph.mat'])
% graphMarkers = plotResults_fullgraph(nodes,cont,ClosedList,index_list,type_list,A,FingerMatrix,pathwayMatrix,graphfig,subFolderName);
%%
for finger = 1:2

%             graphfig{finger} = openfig([pwd '/' subFolderName '/Fullgraph' num2str(finger) '.fig']);

%            graphfig{finger} = openfig([pwd '/' subFolderName '/Nodes s' num2str(finger) '.fig']);

%     set(graphfig{finger},'Position',[(finger-1)*960,41,960,963])
end

%% Algorithm Run:

[Open_list,Closed_list,A,index_list,type_list,path_list,FingerMatrix,pathwayMatrix] = GraphSearch(PG,cont,nodes,starting_node,f1,f2,graphfig);



%%
plotResults_graphSearch(nodes,cont,Closed_list,index_list,type_list,path_list,FingerMatrix,pathwayMatrix,graphfig,subFolderName)

%%
[U,delU,pp] = plotResults_pathPositionsNew(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,subFolderName)
%%
save([pwd '/' subFolderName '/AlgResult.mat'],'U','delU','pp','PG','Open_list','Closed_list','A','index_list','type_list','path_list','FingerMatrix','pathwayMatrix','nodes','cont','basepos','starting_node')
% plotResults_NotpathPositions(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,subFolderName)

%%
for drawResults=1:0 %for loop to allow folding
theta_index = [4,4,4,3,3,3];
if ~exist('graphMarkers','var')
graphMarkers = cell(0);
end
removegraphics = [];
for i=1:numel(graphMarkers)
    if ~isgraphics(graphMarkers{i})
        removegraphics = [removegraphics,i]; %#ok<*AGROW>
    end
end
graphMarkers(removegraphics) = [];
hold on
for i=1:length(Closed_list) % mark the current node of each step of the algorithm
    finger = [];
    cur_node_ind=Closed_list(i);
    connect_ind = find(path_list==cur_node_ind);
    cur_node_ind_intype = index_list(cur_node_ind);
    cur_node_type = type_list(Closed_list(i));
    node_par = nodes{cur_node_type}(cur_node_ind_intype,:);
    s_origin =node_par(1);
    theta_origin = node_par(theta_index(cur_node_type));
    if cur_node_type<4
        set(groot,'CurrentFigure',graphfig{1})
        hold on
        graphMarkers{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3); %#ok<*SAGROW>
        s_origin(2) =node_par(2);
        set(groot,'CurrentFigure',graphfig{2})
        hold on
        graphMarkers{end+1} = plot(s_origin(2),theta_origin,'ok','markerSize',14,'lineWidth',3);
    else
        finger = node_par(end);
        set(groot,'CurrentFigure',graphfig{finger})
        hold on
        graphMarkers{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3);
    end
    for j = 1:length(connect_ind) % mark each neighbor of the current node and connect them
        node_ind = connect_ind(j);
        pathway = pathwayMatrix{cur_node_ind,node_ind};
        node_ind_intype = index_list(node_ind);
        node_type = type_list(node_ind);
        node_par = nodes{node_type}(node_ind_intype,:);
        s =node_par(1);
        theta = node_par(theta_index(node_type));
        if node_type<4 % next node is DS_node
            s(2) =node_par(2);
            if isempty(finger) % current node is DS_node
%                 sMin1 = PG.S(PG.VL(find(PG.S(PG.VL)<s_origin(1),1,'last')));
%                 sMax1 = PG.S(PG.VL(find(PG.S(PG.VL)>s_origin(1),1,'first')));
%                 sMin2 = PG.S(PG.VL(find(PG.S(PG.VL)<s_origin(2),1,'last')));
%                 sMax2 = PG.S(PG.VL(find(PG.S(PG.VL)>s_origin(2),1,'first')));
                if any([1,3]==FingerMatrix(Closed_list(i),connect_ind(j)))
                    set(groot,'CurrentFigure',graphfig{1})
                    hold on
                    graphMarkers{end+1} = plot(s(1),theta,'or','markerSize',14);
                    graphMarkers = plotEdge(pathway,cont,graphMarkers,1);
%                     graphMarkers{end+1} = plot([s(1),s_origin(1)],[theta,theta_origin],'r-','lineWidth',4);
                    
                end
                if any([2,3]==FingerMatrix(Closed_list(i),connect_ind(j)))
                    set(groot,'CurrentFigure',graphfig{2})
                    hold on
                    graphMarkers{end+1} = plot(s(2),theta,'or','markerSize',14);
                    graphMarkers = plotEdge(pathway,cont,graphMarkers,2);
%                     graphMarkers{end+1} = plot([s(2),s_origin(2)],[theta,theta_origin],'r-','lineWidth',4);
                    
                end
            else % current node is SS_node
                set(groot,'CurrentFigure',graphfig{finger})
                hold on
                graphMarkers{end+1} = plot(s(finger),theta,'or','markerSize',14);
                graphMarkers = plotEdge(pathway,cont,graphMarkers,finger);
%                 graphMarkers{end+1} = plot([s(finger),s_origin],[theta,theta_origin],'r-','lineWidth',4);
                
            end
        else % next node is SS_node
            next_finger = node_par(end);
            set(groot,'CurrentFigure',graphfig{next_finger})
            hold on
            s_origin_temp = s_origin;
            if length(s_origin)>1
                s_origin_temp = s_origin(next_finger);
            end
            graphMarkers{end+1} = plot(s,theta,'or','markerSize',14);
            graphMarkers = plotEdge(pathway,cont,graphMarkers,next_finger);
%             graphMarkers{end+1} = plot([s,s_origin_temp],[theta,theta_origin],'r-','lineWidth',4);
            
        end
    end
end
hold off 
for finger = 1:2
saveas(graphfig{finger},[pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.bmp'])
saveas(graphfig{finger},[pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.fig'])
end
%% plot each position of the polygon along the path in separate figures
if 1
xpos = [(0:384:1536),(0:384:1536)];
ypos = [542*ones(1,5),41*ones(1,5)];
for i=1:length(Closed_list)
    f1t = f1;
    f2t = f2;
    connect_ind = find(path_list==Closed_list(i));
    node_ind = index_list(Closed_list(i));
    node_type = type_list(Closed_list(i));
    node_par = nodes{node_type}(node_ind,:);
    s_origin =node_par(1);
    theta_origin = node_par(theta_index(node_type));
    h_origin = node_par(theta_index(node_type)-1);
    if node_type>3 && node_par(end)==2
        f1t = f2;
        f2t = f1;
    end
    fig = drawPolygonNode(PG,s_origin,theta_origin,f1t,f2t);
    if i<=10
    set(fig,'Position',[xpos(i),ypos(i),384,461])
    end
    xlim([-6 6]);
    ylim([-10 6]);
    if i<length(Closed_list)
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin)]);
    else
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin) ', Final']);
    end
    saveas(fig,[pwd '/' subFolderName '/Node ' num2str(i) '.bmp']);
    saveas(fig,[pwd '/' subFolderName '/Node ' num2str(i) '.fig']);
end
end
%% Removing markings from the graphs to reuse them
if ~isempty(graphMarkers) && false
    for i = 1:numel(graphMarkers)
        delete(graphMarkers{i});
    end
end
end
%%
function graphMarkers = plotEdge(pathway,cont,graphMarkers,finger)
for i=1:numel(pathway)-1
    node1 = pathway{i}; %nodes here refer to points along the pathway, not necessarily nodes of the graph
    node2 = pathway{i+1};
    s1 = node1(1,1); theta1 = node1(1,3);
    s2 = node2(1,1); theta2 = node2(1,3);
    if length(node1)==3
        if abs(theta2-theta1)<pi() % if the pathway goes directly over the plane
            graphMarkers{end+1} = plot([s1,s2],[theta1,theta2],'r--','lineWidth',4);
        else %otherwise, it goes to +-pi, and comes from the other (over the cycle of theta)
            thetas = sort([theta1,theta2]);
            graphMarkers{end+1} = plot([s1,s2],[thetas(1),-pi()],'r--','lineWidth',4);
            graphMarkers{end+1} = plot([s1,s2],[thetas(2),pi()],'r--','lineWidth',4);
        end
        % need to add check for theta going the wrong direction
        % (delta_theta>pi)
    else
        seg = cont{node1(end-3)};
        if node1(end-2)>node2(end-2)
            graphMarkers{end+1} = plot(seg(finger,node2(end-2):8:node1(end-2)),seg(4,node2(end-2):8:node1(end-2)),'r--','lineWidth',6,'markerSize',6);
        else
            graphMarkers{end+1} = plot(seg(finger,node1(end-2):8:node2(end-2)),seg(4,node1(end-2):8:node2(end-2)),'r--','lineWidth',6,'markerSize',6);
        end
    end
    
end
end