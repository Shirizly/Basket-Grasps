function [BGdepth,sig,phi,forces] = BasketDepthForces(PG,BGS)
% input: PG - polygon set already at initial orientation, BGS - relative
% positions (s-parameter) of the fingers along the object boundary
% output: BGdepth - depth of the basket grasp at the given positions

f1 = PG.get('1Pos',BGS(1));
f2 = PG.get('1Pos',BGS(2));
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
%% find the reaction forces at the basket grasp
BG = ds_min(start_node(2),:);
p1stat = PG.get('1Pos',BG(1));
p2stat = PG.get('1Pos',BG(2));
theta = BG(4);
n1 = PG.get('normal',BG(1));
n2 = PG.get('normal',BG(2));        
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
p1 = R*(p1stat-PG.com);
p2 = R*(p2stat-PG.com);
% check if either finger supports a vertex
vertexflag = 0;
if any(PG.S(PG.VL)==BG(1))
    vertexflag = vertexflag+1;
end
if any(PG.S(PG.VL)==BG(2))
    vertexflag = vertexflag+2;
end
switch vertexflag % calculate the forces, changes if vertices are supported
    case 0 % E-E BG
        eq = [n1 n2]\[0;1]; % sum of forces equals zero, finds forces mag.
        forces = [eq(1)*n1,eq(2)*n2];
    case 1 % V-E BG
        inter = [0,1]*([n2 -[0;1]]\(-p2)); % intersection point of forces above COM
        n1s = [0;inter]-p1;
        n1s = n1s/norm(n1s);
        eq = [n1s n2]\[0;1]; % sum of forces equals zero, finds forces mag.
        forces = [eq(1)*n1s,eq(2)*n2];
    case 2 % E-V BG
        inter = [0,1]*([n1 -[0;1]]\(-p1)); % intersection point of forces above COM
        n2s = [0;inter]-p2;
        n2s = n2s/norm(n2s);
        eq = [n1 n2s]\[0;1]; % sum of forces equals zero, finds forces mag.
        forces = [eq(1)*n1,eq(2)*n2s];
    case 3 % V-V BG
        % find ranges of heights that can include the force intersection
        % point, for each support seperately
        % first support:
        n11 = n1;
        supn11 = [0,1]*([n11 -[0;1]]\(-p1)); % intersection point of forces above COM
        n12 = PG.get('normal',BG(1)+E-3);
        supn12 = [0,1]*([n12 -[0;1]]\(-p1)); % intersection point of forces above COM
        supn1 = sort([supn11,supn12]);
        % second support:       
        n21 = n2;
        supn21 = [0,1]*([n21 -[0;1]]\(-p2)); % intersection point of forces above COM
        n22 = PG.get('normal',BG(2)+E-3);
        supn22 = [0,1]*([n22 -[0;1]]\(-p2)); % intersection point of forces above COM
        supn2 = sort([supn21,supn22]);
        % find the highest point that intersects both ranges
        inter = min([supn1(2),supn2(2)]);
        % find the normals and complete the equation
        n1s = [0;inter]-p1;
        n1s = n1s/norm(n1s);
        n2s = [0;inter]-p2;
        n2s = n2s/norm(n2s);
        eq = [n1s n2s]\[0;1]; % sum of forces equals zero, finds forces mag.
        forces = [eq(1)*n1s,eq(2)*n2s];
end


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
