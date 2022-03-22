function plotResults_pathPositions(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,subFolderName)
% need to make this work on path_list instead of Closed_list, so it will
% only show nodes on the path
%% plot each position of the polygon along the path in separate figures
lCl = length(Closed_list)
xpos = [(0:384:1536),(0:384:1536)];
ypos = [542*ones(1,5),41*ones(1,5)];
f1 = basepos(:,1);
f2 = basepos(:,2);
theta_index = [4,4,4,3,3,3];

for i=1:lCl
    f1t = f1;
    f2t = f2;
    node_ind = index_list(Closed_list(i));
    node_type = type_list(Closed_list(i));
    node_par = nodes{node_type}(node_ind,:);
    s_origin =node_par(1);
    s2 = [];
    theta_origin = node_par(theta_index(node_type));
    h_origin = node_par(theta_index(node_type)-1);
    if node_type>3 && node_par(end)==2
        f1t = f2;
        f2t = f1;
        s2 = s_origin;
        s_origin = [];
    end
    fig = drawPolygonNode(PG,s_origin,s2,theta_origin,f1t,f2t);
    if i<=10
        set(fig,'Position',[xpos(i),ypos(i),384,461])
    end
    xlim([-6 6]);
    ylim([-8 8]);
    if i<length(Closed_list)
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin)]);
    else
        title(['Node ' num2str(i) ', height = ,' num2str(h_origin) ', Final']);
    end
    saveas(fig,[pwd '/' subFolderName '/Node ' num2str(i) '.bmp']);
    saveas(fig,[pwd '/' subFolderName '/Node ' num2str(i) '.fig']);
end
end