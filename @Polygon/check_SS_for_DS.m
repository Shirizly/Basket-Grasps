function nodes = check_SS_for_DS(PG,nodes,f1,f2)
fingers = [f1,f2];
for i = 4:6 %types of ss_nodes
    ss_nodes = nodes{i};
    ss_nodes_updated = [];
    for j=1:size(ss_nodes,1) % number of nodes in type i
        node = ss_nodes(j,:);
        for k=1:2 % each finger
            if not(check_ss_penetration(PG,node,fingers(:,3-k)-fingers(:,k)))
                ss_nodes_updated(end+1,:) = [node,k];
                ss_nodes_updated(end,2) = ss_nodes_updated(end,2)+fingers(2,k); %updating the height of the node from the assumption of contact at height 0 to contact at finger location
            end
        end 
    end
    nodes{i} = ss_nodes_updated;
end   
end

function check = check_ss_penetration(PG,node,f2_rel)
    cur_s = node(1);
    cur_theta = node(3);
    axle = PG.get('1Pos',cur_s);
    R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
    vertex = R*(PG.vertex-axle);
    vx = vertex(1,:);
    vy = vertex(2,:);
    check = inpolygon(f2_rel(1),f2_rel(2),vx,vy);
end