function [BGdepth,sig,phi] = BasketDepthDifTheta(PG,sig,phi,BGS)
% input: PG - polygon set already at initial orientation, BGS - relative
% positions (s-parameter) of the fingers along the object boundary
% output: BGdepth - depth of the basket grasp at the given positions

f1 = [0;0];
f2 = [sig*cos(phi);sig*sin(phi)];
basepos = [f1,f2];
% identify nodes, create double-support contours
Sigma=inter_finger_distance(PG.X,PG.X);
interfinger = f2-f1;
phi = wrapToPi(atan2(interfinger(2),interfinger(1)));
sig = norm(interfinger);
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
[~,Closed_list,~,index_list,type_list,path_list,~,~] = GraphSearchNoGraphics(PG,cont,nodes,start_node,f1,f2);
%% create a list of nodes on the escape path
current_node = Closed_list(end);
if length(Closed_list)==1
    BGdepth = 0;
else
    theta_index = [4,4,4,3,3,3];
    U = -inf;
    while path_list(current_node)~=0
        current_node = path_list(current_node);
        node_ind = index_list(current_node);
        node_type = type_list(current_node);
        node_par = nodes{node_type}(node_ind,:);
        h_origin = node_par(theta_index(node_type)-1);
        if path_list(current_node)==0
            U0 = h_origin;
        end
        if h_origin>U
            U = h_origin;
        end
    end
    try
        BGdepth = U-U0;
    catch me
        disp('problem with U0')
    end
end
end
