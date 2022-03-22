function [OpenList,ClosedList,A,index_list,type_list,path_list,FingerMatrix,pathwayMatrix] = GraphLocalSearchNoGraphics(PG,cont,nodes,start_node,f1,f2)
% Inputs:
% PG is the polygon
% cont is a cell array of the various segments of the equisigma contour
% nodes is a cell array of the various node matrices, divided by node type
% start_node = [node_type,Node_index]
% f1 and f2 are spatial locations of the fingers
% Outputs:
% OpenList and ClosedList are lists of nodes, closed (checked) and open (unchecked),
% each node defined by their index in index_list and type_list (same number)
% A is a matrix of neighbors (for each two nodes, a boolean variable
% determining if they are neighbors)
% type_list contains the type of each node
% index_list contains the index of each node in the type's matrix (inside
% nodes)
% Further clarifications:
% DS_nodes are [s1,s2,h,theta,index of contour segment, index inside
% contour segment].
% SS_nodes are [s,h,theta,finger], with SS_min being [s,h,theta,type,finger], type = 1 if
% node is hook, 0 if not. finger is the relevant finger for the node
%     labels = {'ds_max','ds_min','ds_virtual','ss_max','ss_min','ss_saddle'}
% path_list is a vector indexed similarly to index_list with describing for
% each visited node the first node that led to it.
%some programming and debugging variables:

epsi = 1E-5;

%% Generate an array of all graph nodes:
% each node is represented by a number, from 1 to #nodes, which serves as
% its index in index_list, type_list, and h_list.
% type_index_list contains the index for the first node of each type.
index_list = []; %index of every node relative to its type
type_list = []; %type of every node
type_index_list = []; %start of every type of nodes in index_list
h_list = []; % list of all heights of nodes.
h_index = [3,3,3,2,2,2]; %the index of height in the node's parameters (depending on node type).
for i = 1:numel(nodes)
    nodes_num = size(nodes{i},1);
    type_list = [type_list,i*ones(1,nodes_num)];
    type_index_list = [type_index_list,size(index_list,2)+1];
    index_list = [index_list,1:1:nodes_num];
    heights = nodes{i}(:,h_index(i)).';
%     if i>3
%         heights = heights+f1(nodes{i}(:,end));
%     end
    h_list = [h_list,heights];
end

%% Finding potential neighbors (all nodes within the neighboring half-strips):
sEnd = PG.S(PG.VL(end));
stripNum = length(PG.VL)-1;
maxnb = size(index_list,2);
path_list = zeros(size(index_list,2),1);
A = false(maxnb);
FingerMatrix = zeros(maxnb);
cur_type = start_node(1);
cur_node = start_node(2)-1+type_index_list(cur_type);
OpenList=cur_node;
ClosedList=[];
pathwayMatrix = cell(maxnb,maxnb);
cur_type = type_list(cur_node);
node_par = nodes{cur_type}(index_list(cur_node),:);
escape = escape_check(PG,node_par,cur_type,f1,f2);
caged = false;
while ~(escape||caged)
    possible_neighbors = []; fingerIndex = [];
    switch cur_type
        case {1,2,3} %DS_max%DS_min%DS_virtual
            node_par = nodes{cur_type}(index_list(cur_node),:);
            for finger = 1:2 % find neighbors in both s-theta planes
                stripIndV = find(PG.S(PG.VL)<=node_par(finger),1,'last');
                if stripIndV>stripNum %the node is at the right edge of the plane
                    stripIndV = stripNum;
                end
                if PG.S(PG.VL(stripIndV))==node_par(finger) %the node is on a strip's edge
                    stripIndV(2) = stripIndV(1)-1; % add the previous strip (left neighboring to the node)
                    if stripIndV(2) < 1 % the node is on the left edge of the plane
                        stripIndV(2) = []; %don't connect over the edge
                    end
                end
                currentTheta = node_par(4);
                for stripInd = stripIndV
                    halfIndDelta = findHalfInd(PG,stripInd,currentTheta);
                    halfInd = sign(halfIndDelta);
                    possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,halfInd,type_index_list,finger)]; %#ok<*AGROW>
                    fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                    if abs(halfIndDelta)<1E-2 || abs(abs(halfIndDelta)-pi())<1E-2
                        possible_neighbors = [possible_neighbors;SearchHalfStrip(PG,nodes,stripInd,-halfInd,type_index_list,finger)];
                        fingerIndex = [fingerIndex, finger*ones(1,length(possible_neighbors)-length(fingerIndex))];
                    end
                end
            end
            fingerIndex(possible_neighbors==cur_node)=[];
            possible_neighbors(possible_neighbors==cur_node)=[];

        case {4,5,6} %SS_max %SS_min %SS_virtual
            node_par = nodes{cur_type}(index_list(cur_node),:);
            finger = node_par(end);
            stripIndV = find(PG.S(PG.VL)<=node_par(1),1,'last');
            if stripIndV>stripNum %the node is at the right edge of the plane
                stripIndV = stripNum;
            end
            if PG.S(PG.VL(stripIndV))==node_par(1) %the node is on a strip's edge
                stripIndV(2) = stripIndV(1)-1; % add the previous strip (left neighboring to the node)
                if stripIndV(2) < 1 % the node is on the left edge of the plane
                    stripIndV(2) = [];
                end
            end
            currentTheta = node_par(3);
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
    neighborsAndfinger = unique([possible_neighbors,fingerIndex.'],'rows');
    possible_neighbors = neighborsAndfinger(:,1);
    fingerIndex = neighborsAndfinger(:,2);

    %% identify neighbors among potential neighbors by constructing a legal path:
    neighbors = []; fingerIndexNeighbors = []; pathways = cell(0);
    for i = 1:length(possible_neighbors)
        finger = fingerIndex(i);
        [connect,pathway] = BuildEdge(PG,nodes,f1,f2,cur_node,possible_neighbors(i),type_list,index_list,cont,finger);
        %         end
        if connect && ~A(cur_node,possible_neighbors(i)) % don't connect again nodes that are already connected
            neighbors = [neighbors,possible_neighbors(i)];
            pathways{end+1} = pathway;
            fingerIndexNeighbors(end+1) = fingerIndex(i);
        end
    end
    
    
    
    %% add neighbors for current nodes on the s-edges of the graph (the duplicates from the other end)
    dupeInd = []; dupeFingerInd = [];
    if cur_type<=3 % DS_nodes
    if node_par(1)==0
        s = sEnd;
        theta = currentTheta;
        finger = 1;
        dupeInd = type_index_list(cur_type)-1+findDSNodeSTheta(nodes{cur_type},s,theta,finger);
        dupeFingerInd = [dupeFingerInd,finger];
    end
    if node_par(1)==sEnd
        s = 0;
        theta = currentTheta;
        finger = 1;
        dupeInd = type_index_list(cur_type)-1+findDSNodeSTheta(nodes{cur_type},s,theta,finger);
        dupeFingerInd = [dupeFingerInd,finger];
    end
    if node_par(2)==0
        s = sEnd;
        theta = currentTheta;
        finger = 2;
        dupeInd = [dupeInd,type_index_list(cur_type)-1+findDSNodeSTheta(nodes{cur_type},s,theta,finger)];
        dupeFingerInd = [dupeFingerInd,finger];
    end
    if node_par(2)==sEnd
        s = 0;
        theta = currentTheta;
        finger = 2;
        dupeInd = [dupeInd,type_index_list(cur_type)-1+findDSNodeSTheta(nodes{cur_type},s,theta,finger)];
        dupeFingerInd = [dupeFingerInd,finger];
    end
    else % SS_nodes
    if node_par(1)==0
        s = sEnd;
        theta = currentTheta;
        finger = node_par(end);
        dupeInd = type_index_list(cur_type)-1+findSSNodeSTheta(nodes{cur_type},s,theta,finger);
        dupeFingerInd = [dupeFingerInd,finger];
    end
    if node_par(1)==sEnd
        s = 0;
        theta = currentTheta;
        finger = node_par(end);
        dupeInd = type_index_list(cur_type)-1+findSSNodeSTheta(nodes{cur_type},s,theta,finger);
        dupeFingerInd = [dupeFingerInd,finger];
    end
    end
    neighbors = [neighbors,dupeInd];
    pathways = [pathways,cell(size(dupeInd))];
    fingerIndexNeighbors = [fingerIndexNeighbors,dupeFingerInd];
    %% add edges to the graph:
    A(cur_node,neighbors) = true;
    A(neighbors,cur_node) = true; % add connection on both directions
    for i = 1:length(neighbors)
        if FingerMatrix(cur_node,neighbors(i))<3
            FingerMatrix(cur_node,neighbors(i)) = FingerMatrix(cur_node,neighbors(i))+fingerIndexNeighbors(i);
            FingerMatrix(neighbors(i),cur_node) = FingerMatrix(neighbors(i),cur_node)+fingerIndexNeighbors(i);
        end
        pathwayMatrix(cur_node,neighbors(i)) = pathways(i);
        pathwayMatrix(neighbors(i),cur_node) = pathways(i);
        if path_list(neighbors(i),1)==0
            % keep info of the first node leading to each new node, to help
            % reconstruct the path
            path_list(neighbors(i),1) = cur_node;
        end
    end
    %% prepare the data structures for the algorithm continuation:
    % Add the neighbors of the current node to the open list:
    OpenList=union(OpenList,setdiff(neighbors.',[OpenList;ClosedList]));
    
    % Remove current node from open list:
    OpenList=setdiff(OpenList,cur_node);
    
    % Add the current node to the end of the closed list:
    ClosedList=[ClosedList;cur_node];
    

    %% Find the next node with smallest height and define it as the new current node:
    [~,min_ind] = min(h_list(OpenList));
    cur_node = OpenList(min_ind);
    if ~isempty(cur_node)
        cur_type = type_list(cur_node);
        node_par = nodes{cur_type}(index_list(cur_node),:);
        if any(cur_type == [1,3,6]) && node_par(3)>2
            escape = true;
        end
    else
        caged = true;
    end
end
% Remove current node from open list:
OpenList=setdiff(OpenList,cur_node);

% Add the current node to the end of the closed list:
ClosedList=[ClosedList;cur_node];


end
%%
function halfIndDelta = findHalfInd(PG,stripInd,currentTheta)
normal = PG.normal(:,stripInd);
thetaSaddle = pi()/2-atan2(normal(2),normal(1)); %angle of the higher saddle line
halfIndDelta = wrapToPi(currentTheta-thetaSaddle)+1E-7;
end

function possible_neighbors = SearchHalfStrip(PG,nodes,stripInd,halfInd,type_index_list,finger)
%% preprocessing of the inputs, preparing the ranges for the search
possible_neighbors = [];
epsi = 3E-3;
sMin = PG.S(PG.VL(stripInd));
sMax = PG.S(PG.VL(stripInd+1));
normal = PG.normal(:,stripInd);
thetaSaddle = wrapToPi(pi()/2-atan2(normal(2),normal(1))); %angle of the higher saddle line
thetaNext = thetaSaddle+halfInd*pi();
thetaNextLimited = max(min(thetaNext,pi()),-pi()); % take either the theta+-pi or the angle limit
thetaRange1 = sort([thetaSaddle, thetaNextLimited]);
thetaRange2 = [-2*pi(),-2*pi()];
if thetaSaddle*halfInd>0
    ThetaEdge = -pi()*halfInd;
    thetaRange2 = sort([ThetaEdge, thetaSaddle-halfInd*pi()]);
end
%% checking DS nodes
for i = 1:3
    sCond = nodes{i}(:,finger)>=sMin & nodes{i}(:,finger)<=sMax;  %the nodes are in the strip
    thetaCond1 = ((nodes{i}(:,4)-thetaRange1(1))>=-epsi & (nodes{i}(:,4)-thetaRange1(2))<=epsi); %the nodes are in the first part of the right theta-range
    thetaCond2 = ((nodes{i}(:,4)-thetaRange2(1))>=-epsi & (nodes{i}(:,4)-thetaRange2(2))<=epsi); % the nodes are in the second part of the right theta-range    thetaCond = thetaCond1 | thetaCond2; %the nodes are in the right theta-range
    thetaCond = thetaCond1 | thetaCond2; %the nodes are in the right theta-range
    pos = find(sCond & thetaCond); % indexes of nodes that are in the right half-strip
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index];
end
%% checking SS nodes
for i = 4:6
    fingerCond = nodes{i}(:,end)==finger; % these are nodes on the right s-theta plane
    sCond = nodes{i}(:,1)>=sMin & nodes{i}(:,1)<=sMax;  %the nodes are in the strip
    thetaCond1 = ((nodes{i}(:,3)-thetaRange1(1))>=-epsi & (nodes{i}(:,3)-thetaRange1(2))<=epsi); %the nodes are in the first part of the right theta-range
    thetaCond2 = ((nodes{i}(:,3)-thetaRange2(1))>=-epsi & (nodes{i}(:,3)-thetaRange2(2))<=epsi); % the nodes are in the second part of the right theta-range
    thetaCond = thetaCond1 | thetaCond2; %the nodes are in the right theta-range
    pos = find(sCond & thetaCond & fingerCond); % indexes of nodes that are in the right half-strip
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index];
end
end


%% building edges
function [connect,pathway] = BuildEdge(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont,finger)
connect = false;intercept = [];
%% define A and B as the top and bottom nodes
cur_node_par = findNodePar(nodes,cur_node,type_list,index_list,finger);
next_node_par = findNodePar(nodes,next_node,type_list,index_list,finger);
if finger == 2 % if we're working around finger 2, switch the order (since all functions that work with the fingers assume primary finger is f1
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
pathway = {A}; %define the first point in the path as A
flag_ADS = A(end)<=3; flag_BDS = B(end)<=3;
if flag_ADS && flag_BDS
    connect = checkContNeighbor(nodes,A,B);
    if connect
        pathway = [pathway,B];
        return
    end
end

%% define deltaTheta and deltaS
deltaTheta = wrapToPi(B(3)-A(3));
deltaS = B(1)-A(1);
if deltaTheta == 0 && deltaS == 0
    connect = 1;
    pathway = [pathway,B];
end
%% generate list of possible height-decreasing directions from A (dirList)
% directions are 1-4 in order: s+,theta+,s-,theta- (going CCW from s+)
possibleDirList = genPosDirList(PG,A);
%% choose the initial direction from dirList (that closes down delta)
requiredDirList = genReqDirList(deltaTheta,deltaS);
dirList = intersect(requiredDirList,possibleDirList);
if isempty(dirList)
    return % no possible path under constraints of the algorithm
end
%% Attempt to build the cardinal direction segments, perform the checks:
if length(dirList)>1 %if there is more than one open direction
    if A(end)==6 %if the A node is a saddle, prefer rotation
        temp = find(mod(dirList,2)==0); %are there any open rotation directions?
        if ~isempty(temp) % if there are, prefer them
            dirList = dirList(temp);
        end
    end
    temp = find(mod(dirList,2)~=0); %are there any open translation directions?
    if ~isempty(temp) % if there are, prefer them
        dirList = dirList(temp);
    end
    dirList = dirList(1); %make sure to have eventually only one starting direction
end
checkRot = false; checkTrans = false;
if any([1,3]==dirList) %starting direction is translation
    [checkTrans,intercept] = check_translation(PG,f1,f2,A,B);
    Atemp = A;
    Atemp(1) = B(1); % "moving" Atemp to the corner of the L movement
%     Atemp(end) = 6; % this is a naive assumption that each middlepoint will be SS and not DS, necessary for genPosDirList()
    if checkTrans
        reqRotDir = requiredDirList(find(mod(requiredDirList,2)==0)); % finds the required rotational direction
        if ~isempty(reqRotDir)
            pathway = [pathway,Atemp(1:3)];
            middlepoint_PosDirList = genPosDirListMiddlePoint(PG,Atemp,finger,cont); % find free downwards directions at the middlepoint
            if any(middlepoint_PosDirList==reqRotDir)
                [checkWater,thetaWater] = checkwatershed(PG,A,Atemp,B); % not crossing the high watershed curve
                [checkRot,intercept] = check_rotation(PG,f1,f2,Atemp,B);
                
                if ~checkRot && ~checkWater % if there's both a crossing and an intercept
                    if ~checkrange(Atemp(3),thetaWater,intercept(2)) % if the intercept is not between the middlepoint and the watershed
                        intercept = []; % the path is not admissible
                    end
                end
                checkRot = checkRot && checkWater;
            else
                return % no possible path under constraints of the algorithm
            end
        else
            checkRot = true;
        end
    end
else % starting direction is rotation
    [checkRot,intercept] = check_rotation(PG,f1,f2,A,B);
    Atemp = A;
    Atemp(3) = B(3); % "moving" Atemp to the corner of the L movement
    if checkRot
        pathway = [pathway,Atemp(1:3)];
        [checkTrans,intercept] = check_translation(PG,f1,f2,Atemp,B);
    end
end
if checkRot&&checkTrans
    connect = true;
    pathway = [pathway,B(1:3)];
    return
end

%% If one of the checks return false, and there was an intercept point, and if B is DS, it may still connect
if flag_BDS && ~isempty(intercept)
    % first, find the intercept point as a point on the contour:
    epsi = 1E-5;
    while length(intercept)<7 % in case the distance to the nearest contour data point is too large
        intercept = findOnCont(intercept,cont,finger,epsi);
        epsi = epsi*10;
    end
    if intercept(2)<=A(2) % the intercept point is not higher than the starting point
        connect = checkContNeighbor(nodes,intercept,B);
    end
    if connect
        pathway = [pathway,intercept,B];
    end
end
end

function possibleDirList = genPosDirList(PG,A)
% function finds cardinal directions that are height-decreasing and allowed.
% at this point nodes have parameters: 4(SS not min),5(SS_min),7(DS)
possibleDirList = [];
epsi = 1E-5;
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
for m=1:4 %go over the directions:
    if m==1 %for s-positive parallel direction (translation):
        if tangent(2)>=-epsi %edge pointing upwards, moving along it in the positive direction decreases com height
            possibleDirList = [possibleDirList,m];
        end
    end
    if m==3 %for s-negative parallel direction (translation):
        if isempty(tangent2) % if the node isn't on a vertex - use the same tangent as (i==1)
            if tangent(2)<=epsi %edge pointing downwards, moving along it in the negative direction decreases com height
                possibleDirList = [possibleDirList,m];
            end
        else
            if tangent2(2)<=epsi %edge pointing downwards, moving along it in the negative direction decreases com height
                possibleDirList = [possibleDirList,m];
            end
        end
    end
    if any([2,4]==m) %for theta parallel directions (rotation):
        rCOM = R*(PG.com-PG.get('1Pos',A(1)));
        rotDir = (3-m); %negative and positive according to right-hand rule;
        lowersaddle = (rCOM(2)<0 && abs(rCOM(1))<0.02); %checking if the node is a lower saddle, no rotation decreases height
        if (rotDir*rCOM(1)<0 || (abs(rCOM(1))<=0.02 && rCOM(2)>=-epsi)) && ~lowersaddle %rotating the object decreases com height,
            % first condition is a reduction of rCOM cross rotDir(as
            % vector), taking only the vertical direction of the
            % instantaneous movement of com.
            % second condition to diff. upper equilibrium (where both
            % directions are h-dec) from lower (where neither are)
            possibleDirList = [possibleDirList,m];
        end
    end
end
switch A(end)
    case {4,5,6} % SS_nodes - only need to find height-decreasing directions:
        % no changes needed, the free possible directions list are indeed
        % free
    case {1,2,3} % DS_nodes - need to remove directions that go through the other finger:
        fRel = R*(PG.get('1Pos',A(6))-PG.get('1Pos',A(1))); %vector from first contact to second
        edgeNum2 = PG.get('edgeNum',A(6));
        normal = PG.normal(:,edgeNum2);
        
        J = [0 -1;1 0]; % rotation matrix for 90 degrees (CCW)
        normal = R*normal;
        normal2 = [];
        if any(PG.S(PG.VL)==A(6)) % if A (second contact) is on a vertex - add the tangent of the previous edge
            if edgeNum2==1 % if A is the first vertex, the previous edge is the last edge
                edgeNum2 = length(PG.VL);
            end
            normal2 = PG.normal(:,edgeNum2-1);
            normal2 = R*normal2;
            isConcave = ((-J*normal).'*normal2)<0; %checks if the vertex is concave
        end
        remove = [];
        for m = 1:length(possibleDirList) % go over the height decreasing directions:
            dir = possibleDirList(m);
            if isempty(normal2)
                % translation directions:
                if dir==1 && normal.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
                    remove = [remove,m];
                end
                if isempty(tangent2) % if the first contact is not on a vertex
                    if dir==3 && normal.'*tangent<-epsi % here the movement direction is tangent
                        remove = [remove,m];
                    end
                else
                    if dir==3 && normal.'*tangent2<-epsi % here the movement direction is tangent
                        remove = [remove,m];
                    end
                end

                if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<-epsi
                    remove = [remove,m];
                else
                    if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<epsi && fRel.'*normal<0
                        remove = [remove,m];
                    end
                end
                
            else
                if isConcave
                    if dir==1 && normal.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
                        remove = [remove,m];
                    else
                        if dir==1 % take into account the limit imposed by the second edge
                            if normal2.'*tangent>epsi
                                remove = [remove,m];
                            end
                        end
                    end
                    if isempty(tangent2) % if the first contact is not on a vertex
                        if dir==3 && (normal.'*tangent<-epsi) %the normal allows it
                            remove = [remove,m];
                        else
                            if dir==3 && normal2.'*tangent<-epsi % take into account the limit imposed by the second edge
                                remove = [remove,m];
                            end
                        end
                    else
                        if dir==3 && (normal.'*tangent2<-epsi) %the normal allows it
                            remove = [remove,m];
                        else
                            if dir==3 && normal2.'*tangent2<-epsi % take into account the limit imposed by the second edge
                                remove = [remove,m];
                            end
                        end
                    end
                    if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<-epsi
                        remove = [remove,m];
                    else
                        if any([2,4]==dir) % take into account the limit imposed by the second edge
                            if ((3-dir)*J*fRel).'*normal2<-epsi
                                remove = [remove,m];
                            end
                        end
                    end
                else % not concave
                    
                    if (dir==1 && normal.'*tangent>epsi) && normal2.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
                        remove = [remove,m];
                    end
                    
                    if isempty(tangent2) % if the first contact is not on a vertex
                        if (dir==3 && normal.'*tangent<-epsi) && normal2.'*tangent<-epsi % here the movement direction is tangent
                            remove = [remove,m];
                        end
                    else
                        if (dir==3 && normal.'*tangent2<-epsi) && normal2.'*tangent2<-epsi % here the movement direction is tangent
                            remove = [remove,m];
                        end
                    end
                    
                    if any([2,4]==dir) && (((3-dir)*J*fRel).'*normal<-epsi) && ((3-dir)*J*fRel).'*normal2<-epsi
                        remove = [remove,m];
                    end
                end
            end
        end
        possibleDirList(remove) = [];
end
end

function possibleDirList = genPosDirListMiddlePoint(PG,A,finger,cont)
% function finds cardinal directions that are allowed. (only for
% middle-points)
% at this point nodes have parameters: 4(SS not min),5(SS_min),7(DS)
epsi = 1E-5;
%% 
cur_s = A(1);
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
%% check whether this is an DS node
possibleDirList = 1:4;
A = A([1,3]);
while epsi<5E-3 && length(A)==2
    A = findOnCont(A,cont,finger,epsi);
    epsi = epsi*10;
end
if length(A)==2
    return
end
%% if this is a DS node, check for free directions
% DS_nodes - need to remove directions that go through the other finger:
fRel = R*(PG.get('1Pos',A(6))-PG.get('1Pos',A(1))); %vector from first contact to second
edgeNum2 = PG.get('edgeNum',A(6));
normal = PG.normal(:,edgeNum2);

J = [0 -1;1 0]; % rotation matrix for 90 degrees (CCW)
normal = R*normal;
normal2 = [];
if any(PG.S(PG.VL)==A(6)) % if A (second contact) is on a vertex - add the tangent of the previous edge
    if edgeNum2==1 % if A is the first vertex, the previous edge is the last edge
        edgeNum2 = length(PG.VL);
    end
    normal2 = PG.normal(:,edgeNum2-1);
    normal2 = R*normal2;
    isConcave = ((-J*normal).'*normal2)<0; %checks if the vertex is concave
end
remove = [];
for m = 1:length(possibleDirList) % go over the directions, remove directions that penetrate:
    dir = possibleDirList(m);
    if isempty(normal2)
        % translation directions:
        if dir==1 && normal.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
            remove = [remove,m];
        end
        if isempty(tangent2) % if the first contact is not on a vertex
            if dir==3 && normal.'*tangent<-epsi % here the movement direction is tangent
                remove = [remove,m];
            end
        else
            if dir==3 && normal.'*tangent2<-epsi % here the movement direction is tangent
                remove = [remove,m];
            end
        end
        
        if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<-epsi
            remove = [remove,m];
        else
            if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<epsi && fRel.'*normal<0
                remove = [remove,m];
            end
        end
        
    else
        if isConcave
            if dir==1 && normal.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
                remove = [remove,m];
            else
                if dir==1 % take into account the limit imposed by the second edge
                    if normal2.'*tangent>epsi
                        remove = [remove,m];
                    end
                end
            end
            if isempty(tangent2) % if the first contact is not on a vertex
                if dir==3 && (normal.'*tangent<-epsi) %the normal allows it
                    remove = [remove,m];
                else
                    if dir==3 && normal2.'*tangent<-epsi % take into account the limit imposed by the second edge
                        remove = [remove,m];
                    end
                end
            else
                if dir==3 && (normal.'*tangent2<-epsi) %the normal allows it
                    remove = [remove,m];
                else
                    if dir==3 && normal2.'*tangent2<-epsi % take into account the limit imposed by the second edge
                        remove = [remove,m];
                    end
                end
            end
            if any([2,4]==dir) && ((3-dir)*J*fRel).'*normal<-epsi
                remove = [remove,m];
            else
                if any([2,4]==dir) % take into account the limit imposed by the second edge
                    if ((3-dir)*J*fRel).'*normal2<-epsi
                        remove = [remove,m];
                    end
                end
            end
        else % not concave
            
            if (dir==1 && normal.'*tangent>epsi) && normal2.'*tangent>epsi % this is condition to remove, since the movement direction is -tangent
                remove = [remove,m];
            end
            
            if isempty(tangent2) % if the first contact is not on a vertex
                if (dir==3 && normal.'*tangent<-epsi) && normal2.'*tangent<-epsi % here the movement direction is tangent
                    remove = [remove,m];
                end
            else
                if (dir==3 && normal.'*tangent2<-epsi) && normal2.'*tangent2<-epsi % here the movement direction is tangent
                    remove = [remove,m];
                end
            end
            
            if any([2,4]==dir) && (((3-dir)*J*fRel).'*normal<-epsi) && ((3-dir)*J*fRel).'*normal2<-epsi
                remove = [remove,m];
            end
        end
    end
end
possibleDirList(remove) = [];

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

function intercept = findOnCont(intercept,cont,finger,epsi)
% epsi = 1E-4;
% if epsi>1E-2
%     disp('suspicious')
% end
for m =1:numel(cont)
    seg = cont{m};
    intercept = intercept(1,:);
    delta = seg([finger,4],:)-repmat(intercept.',1,size(seg,2));
    delta2 = delta.*delta;
    dist = sqrt(delta2(1,:)+delta2(2,:));
    if min(dist)<epsi
        index = find(dist==min(dist),1,'first');
        temp = [seg(1:4,index).',m,index];
        intercept = [temp(finger),temp(3:6),temp(3-finger),7];
        return
    end
end
end
%% checking functions for building edges

function connect = checkContNeighbor(nodes,intercept,B)
%this function finds whether there are other nodes between intercept and B
connect = false;
cur_seg=intercept(end-3);
next_seg = B(end-3);

if cur_seg==next_seg
    contourind = [intercept(end-2),B(end-2)];
    maxind = max(contourind);
    minind = min(contourind);
    for m = 1:3
        pos_between = nodes{m}(:,5:6);
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
connect1 = true;
intercept = [];
interceptv = [];
delta_theta = wrapToPi(B(3)-A(3));
if delta_theta==0
    connect = true;
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
for m=1:size(PG.vertex,2)
    % first find the point on the edge i that will intercept f2
    v1 = vertex(:,m);
    t1 = tangent(:,m);
    vertexTof1O = v1-f1O;
    % construct the polynom in alpha  - the relative distance on the
    % edge from the vertex to the point
    a = norm(t1)^2;
    b = 2*t1.'*vertexTof1O;
    c = norm(vertexTof1O)^2-sigma^2;
    alphav = 0;
    if b==0
        alphav = sqrt(-c/a);
    else
        alphav(1) = (-b+sqrt(b^2-4*a*c))/(2*a); 
        alphav(2) = (-b-sqrt(b^2-4*a*c))/(2*a);
    end
    for i = 1:length(alphav)
        alpha = alphav(i);
        if alpha>=0 && alpha<1 && isreal(alpha) % if the point is inside the edge
            int_point=v1+alpha*t1;
            rel_int_point = int_point-f1O;
            int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point relative to the inter-finger angle
            int_theta = wrapToPi(hand_theta-(int_theta+cur_theta)); % angle of the intercept relative to cur_theta
            if int_theta*delta_theta>0 && (abs(int_theta)-abs(delta_theta)<-1E-5) && abs(int_theta)>1E-5 % if the rotation would intercept the edge before completion
                    interceptv(end+1,:) = [cur_s,int_theta];
                    connect1 = false;
            end
        end
    end
end
if ~connect1
    mintercept = min(abs(interceptv(:,2))); % minimum delta-theta to intercept
    interceptInd = find(abs(interceptv(:,2))==mintercept,1,'first'); % index of first edge to intercept
    intercept = interceptv(interceptInd,:); %#ok<*FNDSB>
    if size(intercept,1)>1
        intercept = intercept(1,:);
    end
    intercept(2) = wrapToPi(intercept(2)+cur_theta);
end
connect2 = false;
if abs(delta_theta-pi())<1E-8||abs(delta_theta+pi())<1E-8 % if we're connecting nodes on the parallel equi-height line (the one representing horizontal orientation of the edge), than both half-strip/directions need to be checked
    [connect2,~] = check_rotation_Exactly_PI(PG,f1,f2,A,B);
end
if connect1||connect2
    connect = true;
    intercept = [];
else
    connect = false;
end
end

function [connect,intercept] = check_rotation_Exactly_PI(PG,f1,f2,A,B)
connect = true;
intercept = [];
interceptv = [];
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
for m=1:size(PG.vertex,2)
    % first find the point on the edge i that will intercept f2
    v1 = vertex(:,m);
    t1 = tangent(:,m);
    vertexTof1O = v1-f1O;
    % construct the polynom in alpha  - the relative distance on the
    % edge from the vertex to the point
    a = norm(t1)^2;
    b = 2*t1.'*vertexTof1O;
    c = norm(vertexTof1O)^2-sigma^2;
    if b==0
        alphav = sqrt(-c/a);
    else
        alphav(1) = (-b+sqrt(b^2-4*a*c))/(2*a);
        alphav(2) = (-b-sqrt(b^2-4*a*c))/(2*a);
    end
    for i = 1:length(alphav)
        alpha = alphav(i);
        if alpha>=0 && alpha<1 && isreal(alpha) % if the point is inside the edge
            int_point=v1+alpha*t1;
            rel_int_point = int_point-f1O;
            int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point relative to the inter-finger angle
            int_theta = wrapToPi(hand_theta-(int_theta+cur_theta)); % angle of the intercept relative to cur_theta
            if int_theta*delta_theta>0 && abs(int_theta)<abs(delta_theta) && abs(int_theta)>1E-5
                interceptv(end+1,:) = [cur_s,int_theta];
                connect = false;
            end
        end
    end
end
if ~connect
    mintercept = min(abs(interceptv(:,2))); % minimum delta-theta to intercept
    interceptInd = find(abs(interceptv(:,2))==mintercept); % index of first edge to intercept
    intercept = interceptv(interceptInd,:); %#ok<*FNDSB>
    if size(intercept,1)>1
        disp('intercept error 1');
    end
    intercept(2) = wrapToPi(intercept(2)+cur_theta);
end
end

function [connect,intercept] = check_translation(PG,f1,f2,A,B)
connect = true;
intercept = [];
epsi = 1E-8;
delta_s = B(1)-A(1);
if delta_s==0
    return
end
cur_theta = A(3);
R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
cur_edge = PG.get('edgeNum',A(1));
if delta_s<0 && abs(PG.S(PG.VL(cur_edge))-A(1))<1E-5 %if we're going backwards from a vertex
    cur_edge = cur_edge-1; % we're actually moving on the previous edge
    if cur_edge == 0 % if the previous edge is the last edge, fix the numbering
        cur_edge = length(PG.VL)-1;
    end
end
if delta_s<0 && abs(PG.S(PG.VL(end))-A(1))<1E-5
    cur_edge = length(PG.VL)-1;
end
trans = R*PG.tangent(:,cur_edge); 
trans = trans/norm(trans);
cur_position = PG.get('1Pos',A(1));
vertex = R*(PG.vertex-cur_position);
tangent = R*PG.tangent;
f2_rel = f2-f1;
interceptv = [];
for m=1:size(PG.vertex,2)
    v1 = vertex(:,m);
    t1 = tangent(:,m);
    mat = [t1 -trans];
    if abs(det(mat))>epsi % don't check the current edge or any edges parallel to it (there will be other inersections if any exist)
        p = mat\(f2_rel-v1); % p(1) is the relative location on the intercepting edge, p(2) is how much deltaS to intercept
        if (p(1)>=-epsi&&p(1)<=1+epsi)&&(abs(p(2))-abs(delta_s)<-1E-5&&p(2)*delta_s>0) % find if the edge intercepts the finger after start of movement and before end.
            interceptv(end+1,:) = [A(1)+p(2),cur_theta];
            connect = false;
        end
    end
end
if ~connect
    mintercept = min(abs(interceptv(:,1)-A(1))); % minimum deltaS to intercept
    interceptInd = find(abs(interceptv(:,1)-A(1))==mintercept); % index of first edge to intercept
    intercept = interceptv(interceptInd,:);
    if size(intercept,1)>1
        if any(intercept(1,:)~=intercept(2,:))
            disp('intercept error 2');
        end
    end
end
end

function [checkWater,thetaWater] = checkwatershed(PG,ref,A,B)
checkWater = true;
% this function checks if there is a high watershed curve between A and B
% this is only relevant after translation, and before checking rotation
s = A(1);
contact_pos = PG.get('1Pos',s);
contact_to_com = PG.com-contact_pos;
thetaWater = atan2(contact_to_com(1),contact_to_com(2));

Delta1 = wrapToPi(B(3)-A(3));
Delta2 = wrapToPi(thetaWater-A(3));
if Delta1*Delta2>0 && abs(Delta1)>abs(Delta2) % if true, the path crosses the watershed
    s_ref = ref(1);
    theta_ref = ref(3);
    contact_pos_ref = PG.get('1Pos',s_ref);
    contact_to_com_ref = PG.com-contact_pos_ref;
    R = [cos(theta_ref) -sin(theta_ref); sin(theta_ref) cos(theta_ref)];
    contact_to_com_ref = R*contact_to_com_ref;
    R = [cos(thetaWater) -sin(thetaWater); sin(thetaWater) cos(thetaWater)];
    contact_to_com = R*contact_to_com;
    if contact_to_com_ref(2)<contact_to_com(2)
        checkWater = false;
    end
end
end

function escape = escape_check(PG,cur_node_par,cur_type,f1,f2)
fingers = [f1,f2];
escape = true;
theta_index = [4,4,4,3,3,3];
cur_theta = cur_node_par(theta_index(cur_type));
% assign the supporting finger, 1 if DS, the relevant finger if SS
finger = 1;
if cur_type>3
    finger = cur_node_par(end);
end

% pose the object at the node's configuration
R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
contact_position = PG.get('1Pos',cur_node_par(1));
vertex = R*(PG.vertex-contact_position)+fingers(:,finger);

% find vertices of the object above the lower finger
low = find(fingers(2,:)==min(fingers(2,:)));
high = 3-low;
lowf = fingers(:,low);
highf = fingers(:,high);
left = find(fingers(1,:)==min(fingers(1,:)));
right = 3-left;
leftf = fingers(:,left);
rightf = fingers(:,right);

posvert = find(vertex(2,:)>lowf(2));
% if any of them are located between the two fingers, it is not yet an
% escape node, as it is possible the object is supported, or may fall to a
% stable grasp.
vertex = [vertex, vertex(:,1)];
    for i=1:length(posvert)
        vert = vertex(:,posvert(i));
        x1 = vert(1);
        y1 = vert(2);
        x2 = vertex(1,posvert(i)+1);
        y2 = vertex(2,posvert(i)+1);
        % is the vertice between (x-axis) the fingers
        c1 = (x1>=leftf(1,:))&(x1<=rightf(1,:));
        % are the vertice and the next one both above the min-finger, to
        % either side of the fingers
        c2 = false;
        if i<length(posvert)
        c2a1 = (x1<=leftf(1,:))&(x2>=rightf(1,:));
        c2a2 = (x2<=leftf(1,:))&(x1>=rightf(1,:));
        c2a = c2a1||c2a2;
        c2b = (posvert(i+1)==posvert(i)+1);
        c2 = c2a&&c2b;
        else % consider the periodicity of the vertices numbering
            if posvert(i)==size(vertex,2)-1 && posvert(1)==1
                c2a1 = (x1<=leftf(1,:))&(x2>=rightf(1,:));
                c2a2 = (x2<=leftf(1,:))&(x1>=rightf(1,:));
                c2 = c2a1||c2a2;
            end
        end
        % is the vertice to one side of the left-finger, the next vertice to the
        % other, and the edge between them crossing between the fingers
        c3a1 = (x1<=leftf(1,:))&(x2>leftf(1,:));
        c3a2 = (leftf(2)-y1)/(y2-y1)*(x2-x1)+x1>leftf(1);
        c3a = c3a1&&c3a2;
        c3b = y1>leftf(2);
        c3 = c3a&&c3b;
        % is the vertice to one side of the right-finger, the next vertice to the
        % other, and the edge between them crossing between the fingers
        c4a1 = (x1>=rightf(1,:))&(x2<rightf(1,:));
        c4a2 = (rightf(2)-y1)/(y2-y1)*(x2-x1)+x1<rightf(1);
        c4a = c4a1&&c4a2;
        c4b = y1>rightf(2);
        c4 = c4a&&c4b;
        if c1||c2||c3||c4
            escape = false;
        end
    end


%     escape = false;


end

function dupeInd = findSSNodeSTheta(nodes,s,theta,finger)
dupeInd = find(nodes(:,1)==s & nodes(:,3)==theta & nodes(:,end)==finger);
end
function dupeInd = findDSNodeSTheta(nodes,s,theta,finger)
dupeInd = find(nodes(:,finger)==s & nodes(:,3)==theta);
end

function check = checkrange(t1,t2,tm)
check = false;
% this function checks if angle tm is between t1 and t2, in the smaller
% range between the two angles.


Delta1 = wrapToPi(t2-t1);
Delta2 = wrapToPi(tm-t1);
if Delta1*Delta2>0 && abs(Delta1)>abs(Delta2)
    check = true;
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

% function relevantFingers = findRelevantFingers(nodes,cur_node,next_node,type_list,index_list)
% node = [cur_node,next_node];
% types = type_list(node);
% if all(types<4) % both nodes are DS, both planes need to be checked
%     relevantFingers = [1,2];
% else % one of the nodes (at least) is SS, need to find the relevant plane
%     ssind = find(types>3,1,'first');
%     ssnode = node(ssind);
%     sstype = types(ssind);
%     ssnode_par = nodes{sstype}(index_list(ssnode),:);
%     relevantFingers = ssnode_par(end);
% end
% end


% function check = inBetweenTheta(theta1,theta2,betweentheta)
% % this functions checks if betweentheta is an angle inside the smallest
% % range connecting theta1 and theta2, not including theta1 and theta2
% % themselves
% check = false;
% theta2 = wrapToPi(theta2-theta1);
% between = wrapToPi(betweentheta-theta1);
% if between*theta2>0 && abs(between)<abs(theta2)
%     check = true;
% end
% end

% function ind = nodeindex(Nodes,
