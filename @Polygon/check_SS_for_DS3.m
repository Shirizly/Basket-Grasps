function nodes = check_SS_for_DS3(PG,nodes,f1,f2)
fingers = [f1,f2];
for i = 4:6 %types of ss_nodes
    ss_nodes = nodes{i};
    ss_nodes_updated = [];
    for j=1:size(ss_nodes,1) % number of nodes in type i
        node = ss_nodes(j,:);
        for k=1:2 % each finger
            if not(check_ss_penetration3(PG,node,fingers(:,3-k)-fingers(:,k),fingers(2,k)))
                ss_nodes_updated(end+1,:) = [node,k];
                ss_nodes_updated(end,2) = ss_nodes_updated(end,2)+fingers(2,k); %updating the height of the node from the assumption of contact at height 0 to contact at finger location
            end
        end 
    end
    nodes{i} = ss_nodes_updated;
end   
for i = 1:3 %types of ds_nodes
    ds_nodes = nodes{i};
    ds_nodes_updated = [];
    for j=1:size(ds_nodes,1) % number of nodes in type i
        node = ds_nodes(j,:);
        for k=1:2 % each finger
            if not(check_ds_penetration3(PG,node,fingers(2,k)))
                ds_nodes_updated(end+1,:) = [node];
            end
        end 
    end
    nodes{i} = ds_nodes_updated;
end 
end

function check = check_ss_penetration3(PG,node,f2_rel,ContactHeight)
    cur_s = node(1);
    cur_theta = node(3);
    axle = PG.get('1Pos',cur_s);
    R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
    vertex = R*(PG.vertex-axle);
    vx = vertex(1,:);
    vy = vertex(2,:);
    check = inpolygon(f2_rel(1),f2_rel(2),vx,vy);
    convertex = vertex(:,PG.CVL(1:end-1)); %relative to contact
    mindisttable = min(convertex(2,:))+ContactHeight;
    check = check || (mindisttable<0);
end

function check = check_ds_penetration3(PG,node,ContactHeight)
    cur_s = node(1);
    cur_theta = node(4);
    axle = PG.get('1Pos',cur_s);
    R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
    vertex = R*(PG.vertex-axle);
    convertex = vertex(:,PG.CVL(1:end-1)); %relative to contact
    mindisttable = min(convertex(2,:))+ContactHeight;
    check =  (mindisttable<0);
end