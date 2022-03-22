function graphMarkers =  plotResults_fullgraph(nodes,cont,Closed_list,index_list,type_list,A,FingerMatrix,pathwayMatrix,graphfig,subFolderName)

% r = 0.4; c = 'w'; centeringx = 0.2; centeringy = -0.02;
% deltax = 0.3; deltay = 0.18;
theta_index = [4,4,4,3,3,3];
epsi = 1E-5;
graphMarkers = [];
topMarkersl = gobjects(0);
topMarkerst = gobjects(0);
for i=1:length(Closed_list) % mark the current node of each step of the algorithm
    finger = [];
    cur_node_ind=Closed_list(i);
    connect_ind = find(A(cur_node_ind,:)>0);
    cur_node_ind_intype = index_list(cur_node_ind);
    cur_node_type = type_list(Closed_list(i));
    cur_node_par = nodes{cur_node_type}(cur_node_ind_intype,:);
    theta_origin = cur_node_par(theta_index(cur_node_type));
    for j = 1:length(connect_ind) % mark each neighbor of the current node and connect them
        node_ind = connect_ind(j);
        pathway = pathwayMatrix{cur_node_ind,node_ind};
        node_ind_intype = index_list(node_ind);
        node_type = type_list(node_ind);
        node_par = nodes{node_type}(node_ind_intype,:);
        s =node_par(1);
        theta = node_par(theta_index(node_type));
        if ~(s == cur_node_par(1) && abs(theta-theta_origin)<epsi)
            
            if node_type<4 % next node is DS_node
                s(2) =node_par(2);
                if isempty(finger) % current node is DS_node
                    if any([1,3]==FingerMatrix(Closed_list(i),connect_ind(j)))
                        set(groot,'CurrentFigure',graphfig{1})
                        hold on
                        %                     graphMarkers{end+1} = plot(s(1),theta,'or','markerSize',14);
                        graphMarkers = plotEdge(pathway,cont,graphMarkers,1);
                        
                    end
                    if any([2,3]==FingerMatrix(Closed_list(i),connect_ind(j)))
                        set(groot,'CurrentFigure',graphfig{2})
                        hold on
                        %                     graphMarkers{end+1} = plot(s(2),theta,'or','markerSize',14);
                        graphMarkers = plotEdge(pathway,cont,graphMarkers,2);
                        
                    end
                else % current node is SS_node
                    set(groot,'CurrentFigure',graphfig{finger})
                    hold on
                    %                 graphMarkers{end+1} = plot(s(finger),theta,'or','markerSize',14);
                    graphMarkers = plotEdge(pathway,cont,graphMarkers,finger);
                    
                end
            else % next node is SS_node
                next_finger = node_par(end);
                set(groot,'CurrentFigure',graphfig{next_finger})
                hold on
                %             graphMarkers{end+1} = plot(s,theta,'or','markerSize',14);
                graphMarkers = plotEdge(pathway,cont,graphMarkers,next_finger);
                
            end
        end
    end
    
    s_origin =cur_node_par(1);
    if cur_node_type<4
        set(groot,'CurrentFigure',graphfig{1})
        hold on
%         graphMarkers{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3); %#ok<*SAGROW>
%         topMarkersl{end+1} = circle(s_origin+deltax+centeringx,theta_origin+deltay+centeringy,r,c);
%         topMarkerst{end+1} = text(s_origin+deltax,theta_origin+deltay,num2str(i),'FontSize',16,'FontWeight','bold');
        s_origin(2) =cur_node_par(2);
        set(groot,'CurrentFigure',graphfig{2})
        hold on
%         graphMarkers{end+1} = plot(s_origin(2),theta_origin,'ok','markerSize',14,'lineWidth',3);
%         topMarkersl{end+1} = circle(s_origin(2)+deltax+centeringx,theta_origin+deltay+centeringy,r,c);
%         topMarkerst{end+1} = text(s_origin(2)+deltax,theta_origin+deltay,num2str(i),'FontSize',16,'FontWeight','bold');
    else
        finger = cur_node_par(end);
        set(groot,'CurrentFigure',graphfig{finger})
        hold on
%         graphMarkers{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3);
%         topMarkersl{end+1} = circle(s_origin+deltax+centeringx,theta_origin+deltay+centeringy,r,c);
%         topMarkerst{end+1} = text(s_origin+deltax,theta_origin+deltay,num2str(i),'FontSize',16,'FontWeight','bold');
    end
end
% for i =1:numel(topMarkersl)
%     uistack(topMarkersl{i},'top')
%     uistack(topMarkerst{i},'top')
% end
hold off
for finger = 1:2
    saveas(graphfig{finger},[pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.bmp'])
    saveas(graphfig{finger},[pwd '/' subFolderName '/Polygon Escape Graph' num2str(finger) '.fig'])
end


%% Removing markings from the graphs to reuse them
if ~isempty(graphMarkers) && false
    for i = 1:numel(graphMarkers)
        delete(graphMarkers{i});
    end
end
end

function graphMarkers = plotEdge(pathway,cont,graphMarkers,finger)
col = 0.5*[1,1,1];
lw = 2;
for i=1:numel(pathway)-1
    node1 = pathway{i}; %nodes here refer to points along the pathway, not necessarily nodes of the graph
    node2 = pathway{i+1};
    s1 = node1(1,1); theta1 = node1(1,3);
    s2 = node2(1,1); theta2 = node2(1,3);
    if length(node1)<7||length(node2)<7
        if abs(theta2-theta1)<pi() % if the pathway goes directly over the plane
            graphMarkers{end+1} = plot([s1,s2],[theta1,theta2],'-','lineWidth',lw,'color',col);
        else %otherwise, it goes to +-pi, and comes from the other (over the cycle of theta)
            thetas = sort([theta1,theta2]);
            graphMarkers{end+1} = plot([s1,s2],[thetas(1),-pi()],'-','lineWidth',lw,'color',col);
            graphMarkers{end+1} = plot([s1,s2],[thetas(2),pi()],'-','lineWidth',lw,'color',col);
        end
        % need to add check for theta going the wrong direction
        % (delta_theta>pi)
    else
        seg = cont{node1(end-3)};
        if node1(end-2)>node2(end-2)
            graphMarkers{end+1} = plot(seg(finger,node2(end-2):1:node1(end-2)),seg(4,node2(end-2):1:node1(end-2)),'.','lineWidth',lw,'markerSize',2*lw,'color',col);
        else
            graphMarkers{end+1} = plot(seg(finger,node1(end-2):1:node2(end-2)),seg(4,node1(end-2):1:node2(end-2)),'.','lineWidth',lw,'markerSize',2*lw,'color',col);
        end
    end
    
end
end
function circles = circle(x,y,r,c)
hold on
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
deltax = xlim(2)-xlim(1);
deltay = ylim(2)-ylim(1);
rx = r;
ry = (deltay/deltax)*r;
th = 0:pi/50:2*pi;
x_circle = rx * cos(th) + x;
y_circle = ry * sin(th) + y;
circles = plot(x_circle, y_circle,c);
fill(x_circle, y_circle, c)

end