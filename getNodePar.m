function nodePar = getNodePar(nodes,index_list,type_list,nodeList)
for i =1:length(nodeList)
    node = nodeList(i);
    nodeInd = index_list(node);
    type = type_list(node);
    temp = nodes{type}(nodeInd,:);
    nodePar(i,1:length(temp)) = temp;
end
% nodePar
end