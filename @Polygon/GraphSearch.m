function [OpenList,ClosedList,A,index_list,type_list,path_list,FingerMatrix] = GraphSearch(PG,cont,nodes,start_node,f1,f2,graphfig)
% Inputs:
% PG is the polygon
% cont is a cell array of the various segments of the equisigma contour
% nodes is a cell array of the various node matrices, divided by node type
% start_node = [node_type,Node_index]
% f1 and f2 are spatial locations of the fingers
% graphfig contains handles for the figures of the s-theta planes

% Outputs:
% OpenList and ClosedList are lists of nodes, closed (checked) and open,
% each node defined by their index in index_list and type_list(same number)
% A is a matrix of neighbors (for each two nodes, a boolean variable
% determining if they are neighbors)
% type_list contains the type of each node
% index_list contains the index of each node in the type's matrix (inside
% nodes)

% Further clarifications:
% DS_nodes are [s1,s2,h,theta,index of contour segment, index inside
% contour segment].
% SS_nodes are [s,h,theta], with SS_min being [s,h,theta,type], type = 1 if
% node is hook, 0 if not.
%     labels = {'ds_max','ds_min','ds_virtual','ss_max','ss_min','ss_saddle'}
%some programming and debugging variables:
debug = 1;
labels = [];
epsi = 1E-5;
%% Generate an array of all graph nodes:
% each node is represented by a number, from 1 to #nodes, which serves as
% its index in index_list and type_list (and h_list).
% type_index_list contains the index for the first node of each type.
index_list = []; %index of every node relative to its type
type_list = []; %type of every node (parallel to index)
type_index_list = []; %start of every group of nodes in index_list
h_list = []; % list of all heights of nodes.
h_index = [3,3,3,2,2,2]; %the index of height in the node\s parameters (depending on node type).
for i = 1:numel(nodes)
    nodes_num = size(nodes{i},1);
    type_list = [type_list,i*ones(1,nodes_num)];
    type_index_list = [type_index_list,size(index_list,2)+1];
    index_list = [index_list,1:1:nodes_num];
    heights = nodes{i}(:,h_index(i)).';
    if i>3
        heights = heights+f1(2);
    end
    h_list = [h_list,heights];
end

%% Finding potential neighbors (all nodes within the neighboring half-strips):
stripNum = length(PG.VL)-1;
maxnb = size(index_list,2);
path_list = zeros(size(index_list,2),1);
A = false(maxnb);
FingerMatrix = zeros(maxnb);
cur_type = start_node(1);
cur_node = start_node(2)-1+type_index_list(cur_type);
OpenList=cur_node;
ClosedList=[];
cur_type = type_list(cur_node);
node_par = nodes{cur_type}(index_list(cur_node),:);
escape = escape_check(PG,node_par,cur_type);
while ~escape
    possible_neighbors = []; fingerIndex = [];
    switch cur_type
        case {1,2,3} %DS_max%DS_min%DS_virtual
            node_par = nodes{cur_type}(index_list(cur_node),:);
            for finger = 1:2 % find neighbors in both s-theta planes
                stripIndV = find(PG.S(PG.VL)<=node_par(finger),1,'last');
                if stripIndV>stripNum %the node is at the right edge of the plane
                    stripIndV = stripNum;
%                     stripIndV(2) = 1; %%%%%%%%current correction for code to run!!!
                end
                if PG.S(PG.VL(stripIndV))==node_par(finger) %the node is on a strip's edge
                    stripIndV(2) = stripIndV(1)-1; % add the previous strip (left neighboring to the node)
                    if stripIndV(2) < 1 % the node is on the left edge of the plane
                        stripIndV(2) = [];%stripNum; %%%%%%%%current correction for code to run!!!
                    end
                end
                currentTheta = node_par(4);
                for stripInd = stripIndV
                    halfIndDelta = findHalfInd(PG,stripInd,currentTheta);
                    halfInd = sign(halfIndDelta);
                    possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,halfInd,type_index_list,finger)]; %#ok<*AGROW>
                    fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                    if abs(halfIndDelta)<epsi
                        possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,-halfInd,type_index_list,finger)];
                        fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                    end
                end
            end
            fingerIndex(possible_neighbors==cur_node)=[];
            possible_neighbors(possible_neighbors==cur_node)=[];
            %             if debug
            %                 labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels);
            %             end
        case {4,5,6} %SS_max %SS_min %SS_virtual
            node_par = nodes{cur_type}(index_list(cur_node),:);
            finger = node_par(end);
            stripIndV = find(PG.S(PG.VL)<=node_par(1),1,'last');
            if stripIndV>stripNum %the node is at the right edge of the plane
                stripIndV = stripNum;
                stripIndV(2) = 1;
            end
            if PG.S(PG.VL(stripIndV))==node_par(1) %the node is on a strip's edge
                stripIndV(2) = stripIndV(1)-1; % add the previous strip (left neighboring to the node)
                if stripIndV(2) < 1 % the node is on the left edge of the plane
                    stripIndV(2) = stripNum;
                end
            end
            currentTheta = node_par(4);
            for stripInd = stripIndV
                halfIndDelta = findHalfInd(PG,stripInd,currentTheta);
                halfInd = sign(halfIndDelta);
                possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,halfInd,type_index_list,finger)]; %#ok<*AGROW>
                fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                if abs(halfIndDelta)<epsi
                    possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,-halfInd,type_index_list,finger)];
                    fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                end
            end
            fingerIndex(possible_neighbors==cur_node)=[];
            possible_neighbors(possible_neighbors==cur_node)=[];
    end
    if ~isempty(labels)
        for i = 1:numel(labels)
            delete(labels{i});
        end
    end
    if debug == 1
        labels = drawpossibleneighbors(PG,nodes,index_list,type_list,cur_node,possible_neighbors,fingerIndex,graphfig);
    end
    %% identify neighbors among potential neighbors by constructing a legal path:
    neighbors = []; fingerIndexNeighbors = [];
    for i = 1:length(possible_neighbors)
        %         disp(i);
        %         relevantFingers = findRelevantFingers(nodes,cur_node,possible_neighbors(i),type_list,index_list);
        connect = false;
        %         for finger = relevantFingers
        finger = fingerIndex(i);
        connect = connect || BuildEdge(PG,nodes,f1,f2,cur_node,possible_neighbors(i),type_list,index_list,cont,finger);
        %         end
        if connect
            neighbors = [neighbors,possible_neighbors(i)];
            fingerIndexNeighbors(end+1) = fingerIndex(i);
        end
    end
    %% add edges to the graph:
    A(cur_node,neighbors) = true;
    A(neighbors,cur_node) = true; % add connection on both directions
    for i = 1:length(neighbors)
        if FingerMatrix(cur_node,neighbors(i))<3
            FingerMatrix(cur_node,neighbors(i)) = FingerMatrix(cur_node,neighbors(i))+fingerIndexNeighbors(i);
            FingerMatrix(neighbors(i),cur_node) = FingerMatrix(neighbors(i),cur_node)+fingerIndexNeighbors(i);
        end
        if path_list(neighbors(i),1)==0
            % keep info of the first node leading to each new node, to help
            % reconstruct the path
            path_list(neighbors(i),1) = cur_node;
        end
    end
    %% prepare the data structures for the algorithm continuation:
    % Add the neighbors of the current node to the open list:
    OpenList=union(OpenList,setdiff(find(A(cur_node,:)).',[OpenList;ClosedList]));
    
    % Remove current node from open list:
    OpenList=setdiff(OpenList,cur_node);
    
    % Add the current node to the end of the closed list:
    ClosedList=[ClosedList;cur_node];
    
    %% Find the next node with smallest height and define it as the new current node:
    [~,min_ind] = min(h_list(OpenList));
    cur_node=OpenList(min_ind);
    cur_type = type_list(cur_node);
    node_par = nodes{cur_type}(index_list(cur_node),:);
    escape = escape_check(PG,node_par,cur_type);
end
% Remove current node from open list:
OpenList=setdiff(OpenList,cur_node);

% Add the current node to the end of the closed list:
ClosedList=[ClosedList;cur_node];

if ~isempty(labels)
    for i = 1:numel(labels)
        delete(labels{i});
    end
end
end

%%
function halfInd = findHalfInd(PG,stripInd,currentTheta)
normal = PG.normal(:,stripInd);
thetaSaddle = pi()/2-atan2(normal(2),normal(1)); %angle of the higher saddle line
halfInd = sign(wrapToPi(currentTheta-thetaSaddle));
end

function possible_neighbors = SearchHalfStrip(PG,nodes,stripInd,halfInd,type_index_list,finger)
%% preprocessing of the inputs, preparing the ranges for the search
possible_neighbors = [];
sMin = PG.S(PG.VL(stripInd));
sMax = PG.S(PG.VL(stripInd+1));
normal = PG.normal(:,stripInd);
thetaSaddle = pi()/2-atan2(normal(2),normal(1)); %angle of the higher saddle line
thetaNext = thetaSaddle+halfInd*pi();
thetaNextLimited = max(min(thetaNext,pi()),-pi());
thetaRange1 = sort([thetaSaddle, thetaNextLimited]);
thetaRange2 = [-2*pi(),-2*pi()];
if thetaSaddle*halfInd>0
    ThetaEdge = -pi()*halfInd;
    thetaRange2 = sort([ThetaEdge, thetaSaddle-halfInd*pi()]);
end
%% checking DS nodes
for i = 1:3
    sCond = nodes{i}(:,finger)>=sMin & nodes{i}(:,finger)<=sMax;  %the nodes are in the strip
    thetaCond1 = (nodes{i}(:,4)>=thetaRange1(1) & nodes{i}(:,4)<=thetaRange1(2)); %the nodes are in the first part of the right theta-range
    thetaCond2 = (nodes{i}(:,4)>=thetaRange2(1) & nodes{i}(:,4)<=thetaRange2(2)); % the nodes are in the second part of the right theta-range
    thetaCond = thetaCond1 | thetaCond2; %the nodes are in the right theta-range
    pos = find(sCond & thetaCond); % indexes of nodes that are in the right half-strip
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index];
end
%% checking SS nodes
for i = 4:6
    fingerCond = nodes{i}(:,end)==finger; % these are nodes on the right s-theta plane
    sCond = nodes{i}(:,1)>=sMin & nodes{i}(:,1)<=sMax;  %the nodes are in the strip
    thetaCond1 = (nodes{i}(:,3)>=thetaRange1(1) & nodes{i}(:,3)<=thetaRange1(2)); %the nodes are in the first part of the right theta-range
    thetaCond2 = (nodes{i}(:,3)>=thetaRange2(1) & nodes{i}(:,3)<=thetaRange2(2)); % the nodes are in the second part of the right theta-range
    thetaCond = thetaCond1 | thetaCond2; %the nodes are in the right theta-range
    pos = find(sCond & thetaCond & fingerCond); % indexes of nodes that are in the right half-strip
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index];
end
end

function connect = BuildEdge(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont,finger)
connect = false;
%% define A and B as the top and bottom nodes
cur_node_par = findNodePar(nodes,cur_node,type_list,index_list,finger);
next_node_par = findNodePar(nodes,next_node,type_list,index_list,finger);
if finger == 2 % if were working around finger 2, switch the order (since all functions that work with the fingers assume primary finger is f1
    f1temp = f1;
    f1 = f2;
    f2 = f1temp;
end
% at this point nodes have parameters: 4(SS not min),5(SS_min),7(DS)
if cur_node_par(2)>next_node_par(2)
    A = cur_node_par;
    B = next_node_par;
else
    A = next_node_par;
    B = cur_node_par;
end
% if (abs(A(1)-11.31)<0.05&&abs(A(3)+1.016)<0.05) % debugging stop
%     disp('here')
% end

flag_dsm=false;
if A(end)==1
    flag_dsm=true;
else
    %% define deltaTheta and deltaS
    deltaTheta = B(3)-A(3);
    deltaS = B(1)-A(1);
    %% generate list of possible height-decreasing directions from A (dirList)
    % directions are 1-4 in order: s+,theta+,s-,theta- (going CCW from s+)
    possibleDirList = genPosDirList(PG,A);
    %% choose the initial direction from dirList (that closes down delta)
    requiredDirList = genReqDirList(deltaTheta,deltaS);
    dirList = intersect(requiredDirList,possibleDirList);
    if isempty(dirList)
        return
    end
    if length(dirList)>1
        dirList = dirList(1);
    end
    %% perform the checks
    checkRot = false; checkTrans = false;
    if any([1,3]==dirList)
        [checkTrans,intercept] = check_translation(PG,f1,f2,A,B);
        Atemp = A;
        Atemp(1) = B(1); % "moving" Atemp to the corner of the L movement
        if checkTrans
            [checkRot,intercept] = check_rotation(PG,f1,f2,Atemp,B);
        end
    else
        [checkRot,intercept] = check_rotation(PG,f1,f2,A,B);
        Atemp = A;
        Atemp(3) = B(3); % "moving" Atemp to the corner of the L movement
        if checkRot
            [checkTrans,intercept] = check_translation(PG,f1,f2,Atemp,B);
        end
    end
    if checkRot&&checkTrans
        connect = true;
        return
    end
end
%% If one of the checks return false, or A is a DS_max, and if B is DS, check if it is the nearest DS node to the intercept point
% first, define the intercept point:
if flag_dsm %the intercept is node A
    intercept = A;
else % find the segment and index on the contour
    intercept = findOnCont(intercept,cont,finger);
end
if length(intercept)<7
    disp('error')
end
connect = checkContNeighbor(nodes,intercept,B);
end

function possibleDirList = genPosDirList(PG,A)
% function finds cardinal directions that are height-decreasing and allowed.
% at this point nodes have parameters: 4(SS not min),5(SS_min),6(DS)
possibleDirList = [];
%% First, find height decreasing directions in the relevant s-theta plane
edgeNum = PG.get('edgeNum',A(1));
tangent = PG.tangent(:,edgeNum);
cur_theta = A(3);
R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
tangent = R*tangent;
tangent2 = [];
if any(PG.S(PG.VL)==A(1)) % if A (first contact) is on a vertex - add the tangent of the previous edge
    if edgeNum==1 % if A is the first vertex, the previous edge is the last edge
        edgeNum = length(PG.VL);
    end
    tangent2 = PG.tangent(:,edgeNum-1);
    tangent2 = R*tangent2;
end
for i=1:4 %go over the directions:
    if i==1 %for s-positive parallel direction (translation):
        if tangent(2)>=0 %edge pointing upwards, moving along it in the positive direction decreases com height
            possibleDirList = [possibleDirList,i];
        end
    end
    if i==3 %for s-negative parallel direction (translation):
        if isempty(tangent2) % if the node isn't on a vertex - use the same tangent as (i==1)
            if tangent(2)<=0 %edge pointing downwards, moving along it in the negative direction decreases com height
                possibleDirList = [possibleDirList,i];
            end
        else
            if tangent2(2)<=0 %edge pointing downwards, moving along it in the negative direction decreases com height
                possibleDirList = [possibleDirList,i];
            end
        end
    end
    if any([2,4]==i) %for theta parallel directions (rotation):
        rCOM = R*(PG.com-PG.get('1Pos',A(1)));
        rotDir = (3-i); %negative and positive according to right-hand rule;
        if rotDir*rCOM(1)<0 || (rotDir*rCOM(1)==0 && rCOM(2)>=0) %rotating the object decreases com height,
            % first condition is a reduction of rCOM cross rotDir(as
            % vector), taking only the vertical direction of the
            % instantaneous movement of com.
            % second condition to diff. upper equilibrium (where both
            % directions are h-dec) from lower (where neither are)
            possibleDirList = [possibleDirList,i];
        end
    end
end
switch length(A)
    case {4,5} % SS_nodes - only need to find height-decreasing directions:
        % no changes needed, the free possible directions list are indeed
        % free
    case 7 % DS_nodes - need to remove directions that go through the other finger:
        fRel = R*(PG.get('1Pos',A(6))-PG.get('1Pos',A(1))); %vector from first contact to second
        edgeNum2 = PG.get('edgeNum',A(6));
        normal = PG.normal(:,edgeNum2);
        
        
        normal = R*normal;
        normal2 = [];
        if any(PG.S(PG.VL)==A(6)) % if A (second contact) is on a vertex - add the tangent of the previous edge
            if edgeNum2==1 % if A is the first vertex, the previous edge is the last edge
                edgeNum2 = length(PG.VL);
            end
            normal2 = PG.normal(:,edgeNum2-1);
            normal2 = R*normal2;
        end
        remove = [];
        for i = 1:length(possibleDirList) % go over the height decreasing directions:
            dir = possibleDirList(i);
            % translation directions:
            if dir==1 && normal.'*tangent>0 % this is condition to remove, since the movement direction is -tangent
                remove = [remove,i];
            else
                if dir==1 && ~isempty(normal2) % take into account the limit imposed by the second edge
                    if normal2.'*tangent>0
                        remove = [remove,i];
                    end
                end
            end
            if dir==3 && (normal.'*tangent<0 && isempty(tangent2))% here the movement direction is tangent
                remove = [remove,i];
            else
                if dir==3 && ~isempty(normal2) % take into account the limit imposed by the second edge
                    if normal2.'*tangent<0
                        remove = [remove,i];
                    end
                end
            end
            if ~isempty(tangent2) % the movement is tangent to the previous edge
                if dir==3 && (normal.'*tangent2<0) %the normal allows it
                    remove = [remove,i];
                else
                    if dir==3 && ~isempty(normal2) % take into account the limit imposed by the second edge
                        if normal2.'*tangent2<0
                            remove = [remove,i];
                        end
                    end
                end
            end
            J = [0 -1;1 0];
            if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<0
                remove = [remove,i];
            else
                if any([2,4]==dir) && ~isempty(normal2) % take into account the limit imposed by the second edge
                    if ((3-dir)*J*fRel).'*normal2<0
                        remove = [remove,i];
                    end
                end
            end
        end
        possibleDirList(remove) = [];
end
end

function requiredDirList = genReqDirList(deltaTheta,deltaS)
% function finds cardinal directions that are needed in order to close
% deltas
requiredDirList = [];
switch sign(deltaS)
    case -1
        requiredDirList = [requiredDirList,3];
    case 1
        requiredDirList = [requiredDirList,1];
end
switch sign(deltaTheta)
    case -1
        requiredDirList = [requiredDirList,4];
    case 1
        requiredDirList = [requiredDirList,2];
end
end

function node_par = findNodePar(nodes,node,type_list,index_list,finger)
type = type_list(node);
node_par = nodes{type}(index_list(node),:);
if length(node_par)==6 % move the s of the other finger to the end.
    temp = node_par;
    node_par(1) = temp(finger);
    node_par(2:end-1) = temp(3:end);
    node_par(end) = temp(3-finger);
end
node_par(end+1) = type; % add the type at the end
end

function intercept = findOnCont(intercept,cont,finger)
epsi = 1E-4;
for i =1:numel(cont)
    seg = cont{i};
    delta = seg([finger,4],:)-repmat(intercept.',1,size(seg,2));
    delta2 = delta.*delta;
    dist = delta2(1,:)+delta2(2,:);
    if min(dist)<epsi
        index = find(dist==min(dist),1,'first');
        temp = [seg(1:4,index).',i,index];
        intercept = [temp(finger),temp(3:6),temp(3-finger),7];
        return
    end
end
end
%%
function connect = checkContNeighbor(nodes,intercept,B)
%this function finds whether there are other nodes between intercept and B
connect = false;
cur_seg=intercept(end-3);
next_seg = B(end-3);

if cur_seg==next_seg
    contourind = [intercept(end-2),B(end-2)];
    maxind = max(contourind);
    minind = min(contourind);
    for i = 1:3
        pos_between = nodes{i}(:,5:6);
        pos_between(pos_between(:,1)~=cur_seg,:)=[];
        in_between = find(pos_between(:,2)<maxind & pos_between(:,2)>minind); %#ok<*EFIND>
        if ~isempty(in_between)
            return
        end
    end
    connect = true;
end
end

function [connect,intercept] = check_rotation(PG,f1,f2,A,B)
connect = true;
intercept = [];
delta_theta = wrapToPi(B(3)-A(3));
if delta_theta==0
    return
end
cur_s = A(1);
cur_theta = A(3);
vertex = PG.vertex;
f1O = PG.get('1Pos',cur_s);
tangent = PG.tangent;
f2_rel = f2-f1;
hand_theta = atan2(f2_rel(2),f2_rel(1));
sigma = norm(f2_rel);
for i=1:size(PG.vertex,2)
    % first find the point on the edge i that will intercept f2
    v1 = vertex(:,i);
    t1 = tangent(:,i);
    vertexTof1O = v1-f1O;
    % construct the polynom in alpha  - the relative distance on the
    % edge from the vertex to the point
    a = norm(t1)^2;
    b = 2*t1.'*vertexTof1O;
    c = norm(vertexTof1O)^2-sigma^2;
    if b==0
        alpha = sqrt(-c/a);
    else
        alpha = (-b+sign(b)*sqrt(b^2-4*a*c))/(2*a); %take the option that comes closer to v1
        % if the closest point on the edge to the axle falls inside the
        % edge (-1<b<0), both intersections need to be checked.
        if (b<0 && b>-1) && (alpha<0 || alpha>=1)
            % if the automatic choice doesn't fit, try the other.
            alpha = (-b+sqrt(b^2-4*a*c))/(2*a);
        end
    end
    if alpha>=0 && alpha<1 && isreal(alpha) % if the point is inside the edge
        int_point=v1+alpha*t1;
        rel_int_point = int_point-f1O;
        int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point relative to the inter-finger angle
        int_theta = wrapToPi(hand_theta-(int_theta+cur_theta)); % angle of the intercept relative to cur_theta
        if int_theta*delta_theta>0 && abs(int_theta)<abs(delta_theta) % if the rotation would intercept the edge before completion
            if abs(delta_theta-pi())<1E-8||abs(delta_theta+pi())<1E-8 % if we're connecting nodes on the parallel equi-height line (the one representing horizontal orientation of the edge), than both half-strip/directions need to be checked
                [connect,intercept] = check_rotation_Exactly_PI(PG,f1,f2,A,B);
                return
            else
                int_theta = wrapToPi(int_theta+cur_theta);
                intercept = [cur_s,int_theta];
                connect = false;
                return
            end
        end
    end
end
end

function [connect,intercept] = check_rotation_Exactly_PI(PG,f1,f2,A,B)
connect = true;
intercept = [];
delta_theta = -wrapToPi(B(3)-A(3));
if delta_theta==0
    return
end
cur_s = A(1);
cur_theta = A(3);
vertex = PG.vertex;
f1O = PG.get('1Pos',cur_s);
tangent = PG.tangent;
f2_rel = f2-f1;
hand_theta = atan2(f2_rel(2),f2_rel(1));
sigma = norm(f2_rel);
for i=1:size(PG.vertex,2)
    % first find the point on the edge i that will intercept f2
    v1 = vertex(:,i);
    t1 = tangent(:,i);
    vertexTof1O = v1-f1O;
    % construct the polynom in alpha  - the relative distance on the
    % edge from the vertex to the point
    a = norm(t1)^2;
    b = 2*t1.'*vertexTof1O;
    c = norm(vertexTof1O)^2-sigma^2;
    if b==0
        alpha = sqrt(-c/a);
    else
        alpha = (-b+sign(b)*sqrt(b^2-4*a*c))/(2*a); %take the option that comes closer to v1
        % if the closest point on the edge to the axle falls inside the
        % edge (-1<b<0), both intersections need to be checked.
        if (b<0 && b>-1) && (alpha<0 || alpha>=1)
            % if the automatic choice doesn't fit, try the other.
            alpha = (-b+sqrt(b^2-4*a*c))/(2*a);
        end
    end
    if alpha>=0 && alpha<1 && isreal(alpha) % if the point is inside the edge
        int_point=v1+alpha*t1;
        rel_int_point = int_point-f1O;
        int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point relative to the inter-finger angle
        int_theta = wrapToPi(hand_theta-(int_theta+cur_theta)); % angle of the intercept relative to cur_theta
        if int_theta*delta_theta>0 && abs(int_theta)<abs(delta_theta)
            int_theta = wrapToPi(int_theta+cur_theta);
            intercept = [cur_s,int_theta];
            connect = false;
            return
        end
    end
end
end

function [connect,intercept] = check_translation(PG,f1,f2,A,B)
connect = true;
intercept = [];
delta_s = B(1)-A(1);
if delta_s==0
    return
end
cur_theta = A(3);
R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
cur_edge = PG.get('edgeNum',A(1));
trans = R*PG.tangent(:,cur_edge);
trans = trans/norm(trans);
cur_position = PG.get('1Pos',A(1));
vertex = R*(PG.vertex-cur_position);
tangent = R*PG.tangent;
f2_rel = f2-f1;
for i=1:size(PG.vertex,2)
    v1 = vertex(:,i);
    t1 = tangent(:,i);
    mat = [t1 -trans];
    if abs(det(mat))>1E-8 % don't check the current edge or any edges parallel to it (there will be other inersections if any exist)
        p = mat\(f2_rel-v1);
        if (p(1)>=0&&p(1)<=1)&&(abs(p(2))<=abs(delta_s)&&p(2)*delta_s>0) % find if the edge intercepts the finger after start of movement and before end.
            intercept = [A(1)+p(2),cur_theta];
            connect = false;
            return
        end
    end
end
end

function escape = escape_check(PG,cur_node_par,cur_type)
escape = true;
theta_index = [4,4,4,3,3,3];
cur_theta = cur_node_par(theta_index(cur_type));
R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
cur_position = PG.get('1Pos',cur_node_par(1));
vertex = R*(PG.vertex-cur_position);
if any(vertex(2,:)>0)
    escape = false;
end
end

%%
% function labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels)
% theta_index = [4,4,4,3,3,3];
% if ~isempty(labels)
%     delete(labels)
% end
% hold on
% for i=1:length(possible_neighbors)
%     node = possible_neighbors(i);
%     node_ind = index_list(node);
%     node_type = type_list(node);
%     node_par = nodes{node_type}(node_ind,:);
%     s =node_par(1);
%     theta = node_par(theta_index(node_type));
%     h = node_par(theta_index(node_type)-1);
%     labels(2*i-1) = plot(s,theta,'ko','markerSize',14,'lineWidth',3);
%     labels(2*i) = text(s,theta+0.2,[num2str(i) ',' num2str(node) ',' num2str(h)] ,'fontSize',10);
% end
% end
function relevantFingers = findRelevantFingers(nodes,cur_node,next_node,type_list,index_list)
node = [cur_node,next_node];
types = type_list(node);
if all(types<4) % both nodes are DS, both planes need to be checked
    relevantFingers = [1,2];
else % one of the nodes (at least) is SS, need to find the relevant plane
    ssind = find(types>3,1,'first');
    ssnode = node(ssind);
    sstype = types(ssind);
    ssnode_par = nodes{sstype}(index_list(ssnode),:);
    relevantFingers = ssnode_par(end);
end
end


function check = inBetweenTheta(theta1,theta2,betweentheta)
% this functions checks if betweentheta is an angle inside the smallest
% range connecting theta1 and theta2, not including theta1 and theta2
% themselves
check = false;
theta2 = wrapToPi(theta2-theta1);
between = wrapToPi(betweentheta-theta1);
if between*theta2>0 && abs(between)<abs(theta2)
    check = true;
end
end


