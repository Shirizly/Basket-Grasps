function plotResults_NotpathPositions(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,subFolderName)
% need to make this work on path_list instead of Closed_list, so it will
% only show nodes on the path
%% create a list of nodes on the escape path
current_node = Closed_list(end);
escape_path = current_node;
while path_list(current_node)~=0
    current_node = path_list(current_node);
    escape_path = [current_node escape_path];
end
not_path = setdiff(Closed_list,escape_path);
not_path = [Closed_list(1),not_path.'];
%% plot each position of the polygon along the not-path in separate subplots
maxlim = zeros(1,4);
lep = length(not_path);
f1 = basepos(:,1);
f2 = basepos(:,2);
theta_index = [4,4,4,3,3,3];
fig = figure();
set(fig,'position',[200,200,960,300]);
for i=1:lep
    
    f1t = f1;
    f2t = f2;
    node_ind = index_list(not_path(i));
    node_type = type_list(not_path(i));
    node_par = nodes{node_type}(node_ind,:);
    s_origin =node_par(1);
    s2 = [];
    theta_origin = node_par(theta_index(node_type));
    h_origin = node_par(theta_index(node_type)-1);
    if node_type>3 && node_par(end)==2
        s2 = s_origin;
        s_origin = [];
    end
    if node_type<4
        s2 = node_par(2);
    end
    ax(i) = subplot(1,lep,i);
    drawPolygonNodePrepare(PG,s_origin,s2,theta_origin,f1t,f2t);
    xlimits = get(ax(i),'xlim');
    ylimits = get(ax(i),'ylim');
    lim = [xlimits,ylimits];
    for j=1:4
        if mod(j,2) == 1
            maxlim(j) = min(maxlim(j),lim(j));
        else
            maxlim(j) = max(maxlim(j),lim(j));
        end
    end
    nodenum = find(Closed_list==not_path(i));
    c = newline;
    h_origin = round(h_origin,2);
    title(['Node ' num2str(nodenum) ',' c 'Height = ,' num2str(h_origin)]);

end
set(ax,'xlim',maxlim(1:2),'ylim',maxlim(3:4));

saveas(fig,[pwd '/' subFolderName '/not-Escape path'  '.bmp']);
saveas(fig,[pwd '/' subFolderName '/not-Escape path'  '.fig']);
end