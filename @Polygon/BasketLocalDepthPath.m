function [BGdepth,sig] = BasketLocalDepthPath(PG,BGS,j)
% input: PG - polygon set already at initial orientation, BGS - relative
% positions (s-parameter) of the fingers along the object boundary
% output: BGdepth - depth of the basket grasp at the given positions

f1 = PG.get('1Pos',BGS(1));
f2 = PG.get('1Pos',BGS(2));
basepos = [f1,f2];
% identify nodes, create double-support contours
Sigma=inter_finger_distance(PG.X,PG.X);

sig = norm(f2-f1);
[cont_original] = PG.GetSigmaContours(Sigma,sig);


cont = PG.CleanContour(cont_original,basepos);

[ds_max,ds_min,ds_virtual] = PG.DSNodes(cont);
[ss_max,ss_min,ss_saddle] = PG.SSNodes();
nodes{1} = ds_max;
nodes{2} = ds_min;
nodes{3} = ds_virtual;
nodes{4} = ss_max;
nodes{5} = ss_min;
nodes{6} = ss_saddle;

diffToBG = ds_min(:,1:2)-BGS.';
diffToBG2 = diffToBG.^2;
distToBG = diffToBG2(:,1)+diffToBG2(:,2);
start_node = [2,find(distToBG == min(distToBG),1,'first')];
nodes = check_SS_for_DS(PG,nodes,f1,f2); 
[Open_list,Closed_list,A,index_list,type_list,path_list,FingerMatrix,pathwayMatrix] = GraphLocalSearchNoGraphics(PG,cont,nodes,start_node,f1,f2);

% 
% for finger = 1:2
%             otherFingerRelative = basepos(:,3-finger)-basepos(:,finger);
%             contactHeight = basepos(2,finger);
%             graphfig{finger} = plotGraphNodes(PG,cont,finger,nodes,otherFingerRelative,contactHeight);
%     set(graphfig{finger},'Position',[(finger-1)*960,41,960,963])
%     set(gcf,'renderer','Painters')
% end
% plotResults_graphSearch(nodes,cont,Closed_list,index_list,type_list,path_list,FingerMatrix,pathwayMatrix,graphfig,[])

[~,BGdepth,~] = plotResults_pathPositionsNew(PG,nodes,Closed_list,index_list,type_list,path_list,basepos,j)
end
