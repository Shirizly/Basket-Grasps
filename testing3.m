% %% Define the polygon
% clear all
clc
% PG = Polygon([1,0;2,0;2,3;-1,4;-2,3;-1,1;-0.5,1].',[0;1.5]);
% PG = Polygon([-1,-1;2,-1;3,-3;3,3;1,1;0,2;-1,1;-2,3;0.2,3;-2,5;-3,1;-2,-4].',[0;0.5]);
% PG = Polygon([-3,0;-1.5,-5;3,-4;3,4;-1.5,5].',[0;0]);
% PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].',[0;0]);
load object.mat
PG = Polygon(object,com);
if ~exist('polyDrawing','var')
    polyDrawing = cell(0);
end
if isempty(polyDrawing)||~isgraphics(polyDrawing)
    polyDrawing = PG.drawPolygon();
end

%% Find equilibrium grasps
res = 50;
[PG,S,X,~] = PG.findBdyVariable(res);

%% draw s-theta plane
f2 = [1;0];
f1 = [-2;1];
basepos = [f1,f2];
%%
Sigma=inter_finger_distance(X,X);
sig = norm(f2-f1);
[cont_original] = PG.GetSigmaContours(Sigma,sig);

cont = PG.CleanContour(cont_original,basepos);
%%
[ds_max,ds_min,ds_virtual] = PG.DSNodes(cont);
[ss_max,ss_min,ss_saddle] = PG.SSNodes();
nodes{1} = ds_max;
nodes{2} = ds_min;
nodes{3} = ds_virtual;
nodes{4} = ss_max;
nodes{5} = ss_min;
nodes{6} = ss_saddle;
% nodes{7} = ss_hook;

nodes = check_SS_for_DS(PG,nodes,f1,f2); %function checks SS_nodes for penetration (of the other finger),
% each node appears once for each finger it is relevant for, with the relevant finger index at the end 
%%
if ~exist('graphfig','var')
    graphfig = cell(1,2);
end

for finger = 1:2
    if isempty(graphfig{finger})||~isgraphics(graphfig{finger})
        if exist(['S' num2str(finger) '_Theta.fig'],'file')
            graphfig{finger} = openfig(['S' num2str(finger) '_Theta.fig']);
        else
            otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
            contactHeight = basepos(2,finger);
            graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight);
        end
        set(graphfig{finger},'Position',[1921+(finger-1)*960,41,960,963])
    end
end


%% Algorithm Run:
starting_node = [2,1];

[Open_list,Closed_list,A,index_list,type_list,path_list,FingerMatrix] = GraphSearch(PG,cont,nodes,starting_node,f1,f2,graphfig);

%% Graphical presentation of the results:
for fakeCounter=1:1
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
for i=1:length(Closed_list)
    finger = [];
    connect_ind = find(path_list==Closed_list(i));
    node_ind = index_list(Closed_list(i));
    cur_node_type = type_list(Closed_list(i));
    node_par = nodes{cur_node_type}(node_ind,:);
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
    for j = 1:length(connect_ind)
        node_ind = index_list(connect_ind(j));
        node_type = type_list(connect_ind(j));
        node_par = nodes{node_type}(node_ind,:);
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
                    graphMarkers{end+1} = plot([s(1),s_origin(1)],[theta,theta_origin],'r-','lineWidth',4);
                    graphMarkers{end+1} = plot(s(1),theta,'or','markerSize',14);
                end
                if any([2,3]==FingerMatrix(Closed_list(i),connect_ind(j)))
                    set(groot,'CurrentFigure',graphfig{2})
                    hold on
                    graphMarkers{end+1} = plot([s(2),s_origin(2)],[theta,theta_origin],'r-','lineWidth',4);
                    graphMarkers{end+1} = plot(s(2),theta,'or','markerSize',14);
                end
            else % current node is SS_node
                set(groot,'CurrentFigure',graphfig{finger})
                hold on
                graphMarkers{end+1} = plot([s(finger),s_origin],[theta,theta_origin],'r-','lineWidth',4);
                graphMarkers{end+1} = plot(s(finger),theta,'or','markerSize',14);
            end
        else % next node is SS_node
            next_finger = node_par(end);
            set(groot,'CurrentFigure',graphfig{next_finger})
            hold on
            s_origin_temp = s_origin;
            if length(s_origin)>1
                s_origin_temp = s_origin(next_finger);
            end
            graphMarkers{end+1} = plot([s,s_origin_temp],[theta,theta_origin],'r-','lineWidth',4);
            graphMarkers{end+1} = plot(s,theta,'or','markerSize',14);
        end
    end
end
hold off 
for finger = 1:2
% set(graphfig{finger},'Position',[0,41,1400,600])

% saveas(graphfig{finger},['Polygon Escape Graph' num2str(finger) '.bmp'])
end
%%
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
    xlim([-4 4]);
    ylim([-5 5]);
    if i<length(Closed_list)
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin)]);
    else
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin) ', Final']);
    end
%     saveas(fig,['Node ' num2str(i) '.bmp']);
end
%%
if ~isempty(graphMarkers) && false
    for i = 1:numel(graphMarkers)
        delete(graphMarkers{i});
    end
end
end