function [U,delU,pp] = plotResults_pathPositionsNew(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,subFolderName)
% need to make this work on path_list instead of Closed_list, so it will
% only show nodes on the path
%% create a list of nodes on the escape path
current_node = Closed_list(end);
escape_path = current_node;
while path_list(current_node)~=0
    current_node = path_list(current_node);
    escape_path = [current_node escape_path];
end
% escape_path(3)=[];
%% plot each position of the polygon along the path in separate figures
maxlim = zeros(1,4);
lep = length(escape_path);
f1 = basepos(:,1);
f2 = basepos(:,2);
theta_index = [4,4,4,3,3,3];
% fig = figure();
% set(fig,'position',[200,200,960,300]);
deltaU = 0;
pp = [];
for i=1:lep
    
    f1t = f1;
    f2t = f2;
    node_ind = index_list(escape_path(i));
    node_type = type_list(escape_path(i));
    node_par = nodes{node_type}(node_ind,:);
    s_origin =node_par(1);
    s2 = [];
    theta_origin = node_par(theta_index(node_type));
    h_origin = node_par(theta_index(node_type)-1);
    if i==1
        U0 = h_origin;
        U = U0;
    end
    if h_origin>U
        U = h_origin;
        pp = node_par;
    end
    if node_type>3 && node_par(end)==2
        s2 = s_origin;
        s_origin = [];
    end
    if node_type<4
        s2 = node_par(2);
    end

%     ax(i) = subplot(1,lep,i);
    figure
    drawPolygonNodePrepare(PG,s_origin,s2,theta_origin,f1t,f2t);
%     xlimits = get(ax(i),'xlim');
%     ylimits = get(ax(i),'ylim');
%     lim = [xlimits,ylimits];
%     for j=1:4
%         if mod(j,2) == 1
%             maxlim(j) = min(maxlim(j),lim(j));
%         else
%             maxlim(j) = max(maxlim(j),lim(j));
%         end
%     end
    nodenum = find(Closed_list==escape_path(i));
    c = newline;
    h_origin = round(h_origin,4);
    set(gca,'visible',false)
    if i<length(escape_path)
        title(['$U = ' num2str(h_origin) '$'],'interpreter','latex','fontsize',26,'visible',1);
    else
        title(['$U = ' num2str(h_origin) '$, escape'],'interpreter','latex','fontsize',26,'visible',1);
    end
delU = U-U0;
if isempty(pp)
        pp = node_par;
end
end
%% print only the pp
if length(pp)<5
    node_type = 5;
else
    node_type = 2;
end
s_origin =pp(1);
s2 = [];
theta_origin = pp(theta_index(node_type));
h_origin = pp(theta_index(node_type)-1);
if node_type>3 && pp(end)==2
    s2 = s_origin;
    s_origin = [];
end
if node_type<4
    s2 = pp(2);
end
figure
drawPolygonNodePrepare(PG,s_origin,s2,theta_origin,f1t,f2t);
set(gca,'visible',false)
title(['point ' num2str(subFolderName)],'interpreter','latex','fontsize',26,'visible',1);
% set(ax,'xlim',maxlim(1:2)+[0,1],'ylim',maxlim(3:4)-[0,2]);
% set(gcf,'renderer','Painters')
% saveas(fig,[pwd '/' subFolderName '/Escape path'  '.bmp']);
% saveas(fig,[pwd '/' subFolderName '/Escape path'  '.fig']);
end