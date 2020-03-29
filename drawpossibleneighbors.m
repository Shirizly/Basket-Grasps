function labels = drawpossibleneighbors(PG,nodes,index_list,type_list,node,possible_neighbors,fingerIndex,graphfig)
theta_index = [4,4,4,3,3,3];
finger = [];
labels = cell(0);
node_ind = index_list(node);
    cur_node_type = type_list(node);
    node_par = nodes{cur_node_type}(node_ind,:);
    s_origin =node_par(1);
    theta_origin = node_par(theta_index(cur_node_type));
    if cur_node_type<4
        set(groot,'CurrentFigure',graphfig{1})
        hold on
        labels{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3);
        s_origin(2) =node_par(2);
        set(groot,'CurrentFigure',graphfig{2})
        hold on
        labels{end+1} = plot(s_origin(2),theta_origin,'ok','markerSize',14,'lineWidth',3);
    else
        finger = node_par(end);
        set(groot,'CurrentFigure',graphfig{finger})
        hold on
        labels{end+1} = plot(s_origin,theta_origin,'ok','markerSize',14,'lineWidth',3);
    end
for j = 1:length(possible_neighbors)
    node_ind = index_list(possible_neighbors(j));
    node_type = type_list(possible_neighbors(j));
    node_par = nodes{node_type}(node_ind,:);
    s =node_par(1);
    theta = node_par(theta_index(node_type));
    if node_type<4 % next node is DS_node
        s(2) =node_par(2); 
        if isempty(finger) % current node is DS_node
%             sMin1 = PG.S(PG.VL(find(PG.S(PG.VL)<s_origin(1),1,'last')));
%             sMax1 = PG.S(PG.VL(find(PG.S(PG.VL)>s_origin(1),1,'first')));
%             sMin2 = PG.S(PG.VL(find(PG.S(PG.VL)<s_origin(2),1,'last')));
%             sMax2 = PG.S(PG.VL(find(PG.S(PG.VL)>s_origin(2),1,'first')));
%             if isempty(sMin1)
                
            if fingerIndex(j)==1 %s(1)>=sMin1 && s(1)<=sMax1
                set(groot,'CurrentFigure',graphfig{1})
                hold on
                labels{end+1} = plot([s(1),s_origin(1)],[theta,theta_origin],'r-','lineWidth',4); %#ok<*AGROW>
                labels{end+1} = plot(s(1),theta,'or','markerSize',14);
            end
            if fingerIndex(j)==2 %s(2)>=sMin2 && s(2)<=sMax2
                set(groot,'CurrentFigure',graphfig{2})
                hold on
                labels{end+1} = plot([s(2),s_origin(2)],[theta,theta_origin],'r-','lineWidth',4);
                labels{end+1} = plot(s(2),theta,'or','markerSize',14);
            end
        else % current node is SS_node
            set(groot,'CurrentFigure',graphfig{finger})
            hold on
            labels{end+1} = plot([s(finger),s_origin],[theta,theta_origin],'r-','lineWidth',4);
            labels{end+1} = plot(s(finger),theta,'or','markerSize',14);
        end
    else % next node is SS_node
        next_finger = node_par(end);
        set(groot,'CurrentFigure',graphfig{next_finger})
        hold on
        s_origin_temp = s_origin;
        if length(s_origin)>1
            s_origin_temp = s_origin(next_finger);
        end
        labels{end+1} = plot([s,s_origin_temp],[theta,theta_origin],'r-','lineWidth',4);
        labels{end+1} = plot(s,theta,'or','markerSize',14);
    end 
end
end