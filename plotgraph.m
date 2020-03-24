function plotgraph(contv,graph_nodes,A)
hold on
for i = 1:size(A,1)
    
    switch graph_nodes.type(i)
        case 1
            plota(contv(:,graph_nodes.index(i)),'+r','markersize',6,'linewidth',3);
            textpos = contv(:,graph_nodes.index(i));
            textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
            text(textpos(1),textpos(2),['h = ' num2str(round(graph_nodes.h(i),2))],'fontweight','bold');
        case 0
            plota(contv(:,graph_nodes.index(i)),'*r','markersize',6,'linewidth',3);
            textpos = contv(:,graph_nodes.index(i));
            textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
            text(textpos(1),textpos(2),['h = ' num2str(round(graph_nodes.h(i),2))],'fontweight','bold');
        case -1
            plota(contv(:,graph_nodes.index(i)),'+g','markersize',6,'linewidth',3);
            textpos = contv(:,graph_nodes.index(i));
            textpos = textpos + (0.5*[2;0]*(textpos(1)<textpos(2))+1*[-2;0]*(textpos(1)>textpos(2)));
            text(textpos(1),textpos(2),['h = ' num2str(round(graph_nodes.h(i),2))],'fontweight','bold');
    end
    for j = 1:i
        if A(i,j) == 1
            plota([contv(:,graph_nodes.index(i)) contv(:,graph_nodes.index(j))],'k');
        end
    end
end
hold off
end