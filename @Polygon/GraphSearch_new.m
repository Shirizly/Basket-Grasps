function [Open_list,Closed_list,A,index_list,type_list,path_list] = GraphSearch_new(PG,cont,nodes,start_node,f1,f2)
%
% PG.S,PG.X are discretization of the boundary
% contv are the (s1,s2) of each point on the equisigma contour
% start_node = [node_type,Node_index]

% O and A are lists of neigbours, closed (checked) and open
% NB is a matrix of neighborhoods

% DS_nodes are [s1,s2,h,theta,index of contour segment, index inside
% contour segment].
% SS_nodes are [s,h,theta], with SS_min being [s,h,theta,type], type = 1 if
% node is hook, 0 if not.
%     labels = {'ds_max','ds_min','ds_virtual','ss_max','ss_min','ss_saddle'}
debug = 0;

%% Generate an array of all graph nodes:
index_list = []; %index of every node relative to its type
type_list = []; %type of every node (parallel to index)
type_index_list = []; %start of every group of nodes in index_list
h_list = []; % list of all heights of nodes.
h_index = [3,3,3,2,2,2];
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


%%
maxnb = size(index_list,2);
path_list = zeros(size(index_list,2),1);
A = false(maxnb);
cur_type = start_node(1);
cur_node = start_node(2)-1+type_index_list(cur_type);
Open_list=cur_node;
Closed_list=[];
labels = [];
escape = false;
check_functions={@not_connect,@connect_contour,@connect_contour,@connect_monotone_rotation,@connect_monotone_rotation,@not_connect;
    @connect_contour,@connect_monotone_DS_DS,@connect_contour,@connect_monotone_rotation,@connect_monotone_DSmin_SS,@connect_monotone_DSmin_SS;
    @connect_contour,@connect_contour,@connect_contour,@connect_monotone_rotation,@connect_monotone_rotation,@connect_monotone_DSmin_SS;
    @not_connect,@not_connect,@connect_monotone_rotation,@not_connect,@connect_monotone_rotation,@connect_monotone_rotation;
    @connect_monotone_rotation,@connect_monotone_DSmin_SS,@connect_monotone_rotation,@connect_monotone_rotation,@connect_monotone_SSmin_SSmin,@connect_monotone_SSmin_SSsaddle;
    @not_connect,@connect_monotone_DSmin_SS,@connect_monotone_DSmin_SS,@connect_monotone_rotation,@connect_monotone_SSmin_SSsaddle,@connect_monotone_rotation};
cur_type = type_list(cur_node);
node_par = nodes{cur_type}(index_list(cur_node),:);
escape = escape_check(PG,node_par,cur_type);
while ~escape
    possible_neighbors = [];
    switch cur_type

        case {1,2,3} %DS_max%DS_min%DS_virtual - separated to allow change later to more efficient form
            for finger = 1:2 % find neighbors in both s-theta planes
                node_par = nodes{cur_type}(index_list(cur_node),:);
                smin = PG.S(PG.VL(find(PG.S(PG.VL)<node_par(finger),1,'last')));
                smax = PG.S(PG.VL(find(PG.S(PG.VL)>node_par(finger),1,'first')));
                if isempty(smin)
                    smintemp = PG.S(PG.VL(end-1));
                    smaxtemp = PG.S(PG.VL(end));
                    possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger)];
                    possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                    smin = 0;
                end
                if isempty(smax)
                    smintemp = 0;
                    smaxtemp = PG.S(PG.VL(2));
                    possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger)];
                    possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                    smax = PG.S(end);
                end
                possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smin,smax,type_index_list,finger)];
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smin,smax,type_index_list)];
                
                possible_neighbors(possible_neighbors==cur_node)=[];
                if debug
                    labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels);
                end
            end
        case {4,5} %SS_max %SS_min
            %% find neigboring strips:
            node_par = nodes{cur_type}(index_list(cur_node),:);
            finger = node_par(end);
            smin = PG.S(PG.VL(find(PG.S(PG.VL)<node_par(1),1,'last')));
            smax = PG.S(PG.VL(find(PG.S(PG.VL)>node_par(1),1,'first')));
            if isempty(smin)
                smintemp = PG.S(PG.VL(end-1));
                smaxtemp = PG.S(PG.VL(end));
                possible_neighbors = search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger);
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                smin = 0;
            end
            if isempty(smax)
                smintemp = 0;
                smaxtemp = PG.S(PG.VL(2));
                possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger)];
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                smax = PG.S(end);
            end
                possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smin,smax,type_index_list,finger)];
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smin,smax,type_index_list)];
            
                possible_neighbors(possible_neighbors==cur_node)=[];
                if debug
                    labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels);
                end
        case 6 % SS_saddle
            node_par = nodes{cur_type}(index_list(cur_node),:);
            finger = node_par(end);
            smin = PG.S(PG.VL(find(PG.S(PG.VL)<node_par(1),1,'last')));
            smax = PG.S(PG.VL(find(PG.S(PG.VL)>node_par(1),1,'first')));
            if isempty(smin)
                smintemp = PG.S(PG.VL(end-1));
                smaxtemp = PG.S(PG.VL(end));
                possible_neighbors = search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger);
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                smin = 0;
            end
            if isempty(smax)
                smintemp = 0;
                smaxtemp = PG.S(PG.VL(2));
                possible_neighbors = search_SSnodes(nodes,smintemp,smaxtemp,type_index_list,finger);
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smintemp,smaxtemp,type_index_list)];
                smax = PG.S(end);
            end
                possible_neighbors = [possible_neighbors;search_SSnodes(nodes,smin,smax,type_index_list,finger)];
                possible_neighbors = [possible_neighbors;search_DSnodes(nodes,smin,smax,type_index_list)];
                possible_neighbors(possible_neighbors==cur_node)=[];
                if debug
                    labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels);
                end
    end
    neighbors = [];
    for i = 1:length(possible_neighbors)
%         disp(i);
        poss_type = type_list(possible_neighbors(i));
        check_function = check_functions{cur_type,poss_type};
        connect = check_function(PG,nodes,f1,f2,cur_node,possible_neighbors(i),type_list,index_list,cont);
        if connect
            neighbors = [neighbors,possible_neighbors(i)];
        end
    end
   
    A(cur_node,neighbors) = true;
    A(neighbors,cur_node) = true; % add connection on both directions
    for i = 1:length(neighbors)
        if path_list(neighbors(i),1)==0 
            % keep info of the first node leading to each new node, to help
            % reconstruct the path
            path_list(neighbors(i),1) = cur_node;
        end
    end
    % Add the neighbors of the current node to the open list:
    Open_list=union(Open_list,setdiff(find(A(cur_node,:)).',[Open_list;Closed_list]));
    
    % Remove current node from open list:
    Open_list=setdiff(Open_list,cur_node);
    
    % Add the current node to the end of the closed list:
    Closed_list=[Closed_list;cur_node];
    
    % Find the next node with smallest height and define it as the new current
    % node:
    [~,min_ind] = min(h_list(Open_list));
    cur_node=Open_list(min_ind);
    cur_type = type_list(cur_node);
    node_par = nodes{cur_type}(index_list(cur_node),:);
    escape = escape_check(PG,node_par,cur_type);
end
% Remove current node from open list:
Open_list=setdiff(Open_list,cur_node);

% Add the current node to the end of the closed list:
Closed_list=[Closed_list;cur_node];

end

%%
function possible_neighbors = search_SSnodes(nodes,smin,smax,type_index_list,finger)
% Need to consider the existence of first vertex nodes as last vertex nodes
% as well
possible_neighbors = [];
for i=4:6
    pos = find(nodes{i}(:,1)>=smin & nodes{i}(:,1)<=smax & nodes{i}(:,end)==finger);
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index]; 
end
max_s = max(nodes{4}(:,1));
if smax == max_s
    for i=4:5
        pos = find(nodes{i}(:,1)==0 & nodes{i}(:,end)==finger);
        pos_index = type_index_list(i)-1+pos;
        possible_neighbors = [possible_neighbors;pos_index];
    end
end

end
 
function possible_neighbors = search_DSnodes(nodes,smin,smax,type_index_list)
possible_neighbors = [];
for i=1:3 % no need to add ds_max, they only connect to DS nodes
    pos = find(nodes{i}(:,1)>=smin & nodes{i}(:,1)<=smax);
    pos_index = type_index_list(i)-1+pos;
    possible_neighbors = [possible_neighbors;pos_index];
end

end

%%
function connect = connect_monotone_rotation(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont)
% this function attempts to connect 2 nodes over s-theta plane, using only
% pure rotation and pure translation
connect = false;

%% check whether there is really only a need for pure rotation:
cur_type = type_list(cur_node);
cur_node_par = nodes{cur_type}(index_list(cur_node),:);
next_type = type_list(next_node);
next_node_par = nodes{next_type}(index_list(next_node),:);

%% if one of the nodes is SS and the other is DS, and the ss node relates to the second finger, switch the order of s1,s2 for the DS node
if cur_type <4 &&next_type>3 % cur is DS, next is SS
    if next_node_par(end)==2 % the finger index points to the correct s_index
        cur_node_par([1,2]) = cur_node_par([2,1]);
        temp = f2;
        f2 = f1;
        f1 = temp;
    end
end
if cur_type >3 &&next_type<4 % cur is SS, next is DS
    if cur_node_par(end)==2 % the finger index points to the correct s_index
        next_node_par([1,2]) = next_node_par([2,1]);
        temp = f2;
        f2 = f1;
        f1 = temp;
    end
end
    

%% check if same s, if so, check rotation
if cur_node_par(1)==next_node_par(1)
    connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
end

end

function connect = connect_monotone_SSmin_SSmin(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont)
% this function attempts to connect 2 SSmin nodes, on neighboring edges of
% a strip
connect = false;
theta_index = [4,4,4,3,3,3];
%% check whether the nodes are 2 SS_min with a saddle in between
ss_saddle = nodes{6};
cur_type = type_list(cur_node);
cur_node_par = nodes{cur_type}(index_list(cur_node),:);
next_type = type_list(next_node);
next_node_par = nodes{next_type}(index_list(next_node),:);
if cur_node_par(1)==0 && next_node_par(1)>=PG.S(PG.VL(end-1))
    cur_node_par(1)=PG.S(PG.VL(end));
end
first_s = min([cur_node_par(1),next_node_par(1)]);
last_s = max([cur_node_par(1),next_node_par(1)]);
len_sad = size(ss_saddle,1)/2;
possible_saddle1 = find(ss_saddle(1:len_sad,1)>first_s,1,'first');
possible_saddle2 = find(ss_saddle(1:len_sad,1)<last_s,1,'last');
if possible_saddle1==possible_saddle2 %there is a saddle between them, no connection
    return
end

%% if there isn't, find the lower node:
lower_node_h = min([cur_node_par(2),next_node_par(2)]);
if cur_node_par(2)==lower_node_h
    connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
    connect = connect&&check_translation(PG,f1,f2,next_node_par,cur_type,cur_node_par);
else
    connect = check_translation(PG,f1,f2,cur_node_par,cur_type,next_node_par);
    connect = connect&&check_rotation(PG,f1,f2,next_node_par,next_type,cur_node_par,cur_type);
end
end

function connect = connect_monotone_DSmin_SS(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont)
% this function attempts to connect 2 nodes over s-theta plane, using only
% pure rotation and pure translation
connect = false;
%% some preparations:
theta_index = [4,4,4,3,3,3];
cur_type = type_list(cur_node);
cur_node_par = nodes{cur_type}(index_list(cur_node),:);
next_type = type_list(next_node);
next_node_par = nodes{next_type}(index_list(next_node),:);
if cur_node_par(1)==0 && next_node_par(1)>=PG.S(PG.VL(end-1))
    cur_node_par(1)=PG.S(PG.VL(end));
end
%% if one of the nodes is SS and the other is DS, and the ss node relates to the second finger, switch the order of s1,s2 for the DS node, and switch locations of the fingers
if cur_type <4 &&next_type>3 % cur is DS, next is SS
    if next_node_par(end)==2 % the finger index points to the correct s_index
        cur_node_par([1,2]) = cur_node_par([2,1]);
        temp = f2;
        f2 = f1;
        f1 = temp;
    end
end
if cur_type >3 &&next_type<4 % cur is SS, next is DS
    if cur_node_par(end)==2 % the finger index points to the correct s_index
        next_node_par([1,2]) = next_node_par([2,1]);
        temp = f2;
        f2 = f1;
        f1 = temp;
    end
end

%% check whether there is only a need for pure rotation or translation:

% check if same theta
if cur_node_par(theta_index(cur_type))==next_node_par(theta_index(next_type))
    connect = check_translation(PG,f1,f2,cur_node_par,cur_type,next_node_par);
    return
end
% check if same s
if cur_node_par(1)==next_node_par(1)
    connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
    return
end
%% if there isn't, find the horizontal edge's theta:
cur_edge = find(PG.S(PG.VL)<mean([cur_node_par(1) next_node_par(1)]),1,'last');
horiz_theta = -atan2(PG.tangent(2,cur_edge),PG.tangent(1,cur_edge));
%% if horiz_theta is between the 2 nodes
theta1 = cur_node_par(theta_index(cur_type));
theta2 = next_node_par(theta_index(next_type));
if inBetweenTheta(theta1,theta2,horiz_theta)
    return
end

%% if horiz_theta isn't between the 2 nodes
lower_node_h = min([cur_node_par(2),next_node_par(2)]);
if cur_node_par(theta_index(cur_type)-1)==lower_node_h
    connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
    connect = connect&&check_translation(PG,f1,f2,next_node_par,cur_type,cur_node_par);
else
    connect = check_translation(PG,f1,f2,cur_node_par,cur_type,next_node_par);
    connect = connect&&check_rotation(PG,f1,f2,next_node_par,next_type,cur_node_par,cur_type);
end
end

function connect = connect_monotone_SSmin_SSsaddle(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont)
% this function attempts to connect 2 nodes over s-theta plane, using only
% pure rotation and pure translation
connect = false;
theta_index = [4,4,4,3,3,3];
%% check whether the saddle is the lower one
cur_type = type_list(cur_node);
cur_node_par = nodes{cur_type}(index_list(cur_node),:);
next_type = type_list(next_node);
next_node_par = nodes{next_type}(index_list(next_node),:);
if cur_node_par(1)==0 && next_node_par(1)>=PG.S(PG.VL(end-1))
    cur_node_par(1)=PG.S(PG.VL(end));
end
ss_saddle = nodes{6};
if cur_type ==6
    both_saddles = find(ss_saddle(:,1)==cur_node_par(1));
    if min(ss_saddle(both_saddles,2))<cur_node_par(2) % if our saddle is the higher one, don't connect to min
        return
    end
    % if it is the lower one, connect through translation then rotation:
    connect = check_translation(PG,f1,f2,cur_node_par,cur_type,next_node_par);
    connect = connect&&check_rotation(PG,f1,f2,next_node_par,cur_type,cur_node_par,next_type);
    return
end
%else 
    both_saddles = find(ss_saddle(:,1)==next_node_par(1));
    if min(ss_saddle(both_saddles,2))<next_node_par(2) % if our saddle is the higher one, don't connect
        return
    end
    % if it is the lower one, connect through rotation then translation:
    connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
    connect = connect&&check_translation(PG,f1,f2,next_node_par,cur_type,cur_node_par);
end

function connect = not_connect(varargin)
connect = false;
end

function connect = connect_contour(PG,nodes,f1,f2,cur_node,next_node,type_list,index_list,cont)
%% need to write
% check if the possible neighbor is directly next to the current one on a
% contour. If not, use a function similar to connect_monotone_DSmin_SS to
% try connecting them over each SS space, than check to make sure the
% connection went over the higher plane.
connect = false;
cur_type = type_list(cur_node);
cur_node_par = nodes{cur_type}(index_list(cur_node),:);
cur_seg=cur_node_par(5);
next_type = type_list(next_node);
next_node_par = nodes{next_type}(index_list(next_node),:);
next_seg = next_node_par(5);
seg = cont{cur_node_par(5)};
contourind = [cur_node_par(6),next_node_par(6)];
maxind = max(contourind);
minind = min(contourind);

if cur_seg==next_seg 
    for i = 1:3
        pos_between = nodes{i}(:,5:6);
        pos_between(pos_between(:,1)~=cur_seg,:)=[];
        in_between = find(pos_between(:,2)<maxind & pos_between(:,2)>minind);
        if ~isempty(in_between)
            return
        end
    end   
    connect = true;
end

%% add bit for if the cur_node is at the end of the segment
end

function connect = connect_monotone_DS_DS(varargin)
connect = false;
end

%%
function connect = check_rotation(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type)
    connect = true;
    theta_index = [4,4,4,3,3,3];
    delta_theta = wrapToPi(next_node_par(theta_index(next_type))-cur_node_par(theta_index(cur_type)));
    if delta_theta==0
        return
    end
    cur_s = cur_node_par(1);
    cur_theta = cur_node_par(theta_index(cur_type));
    vertex = PG.vertex;
    axle = PG.get('1Pos',cur_s);
    tangent = PG.tangent;
    f2_rel = f2-f1;
    hand_theta = atan2(f2_rel(2),f2_rel(1));
    sigma = norm(f2_rel);
    for i=1:size(PG.vertex,2)
        v1 = vertex(:,i);
        t1 = tangent(:,i);
        A = v1-axle;
        a = norm(t1)^2;
        b = 2*t1.'*A;
        c = norm(A)^2-sigma^2;
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
        if alpha>=0 && alpha<1 && isreal(alpha)
            int_point=v1+alpha*t1;
            rel_int_point = int_point-axle;
            int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point in the original position of the polygon
            int_theta = wrapToPi(hand_theta-(int_theta+cur_theta));
            if int_theta*delta_theta>0 && abs(int_theta)<abs(delta_theta)
                if abs(delta_theta-pi())<1E-8||abs(delta_theta+pi())<1E-8
                    connect = check_rotation_Exactly_PI(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type);
                    return
                else
                    connect = false;
                    return
                end
            end
        end
    end
end

function connect = check_rotation_Exactly_PI(PG,f1,f2,cur_node_par,cur_type,next_node_par,next_type)
    connect = true;
    theta_index = [4,4,4,3,3,3];
    delta_theta = -wrapToPi(next_node_par(theta_index(next_type))-cur_node_par(theta_index(cur_type)));
    if delta_theta==0
        return
    end
    cur_s = cur_node_par(1);
    cur_theta = cur_node_par(theta_index(cur_type));
    vertex = PG.vertex;
    axle = PG.get('1Pos',cur_s);
    tangent = PG.tangent;
    f2_rel = f2-f1;
    hand_theta = atan2(f2_rel(2),f2_rel(1));
    sigma = norm(f2_rel);
    for i=1:size(PG.vertex,2)
        v1 = vertex(:,i);
        t1 = tangent(:,i);
        A = v1-axle;
        a = norm(t1)^2;
        b = 2*t1.'*A;
        c = norm(A)^2-sigma^2;
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
        if alpha>=0 && alpha<1 && isreal(alpha)
            int_point=v1+alpha*t1;
            rel_int_point = int_point-axle;
            int_theta = atan2(rel_int_point(2),rel_int_point(1)); % angle of the intercept point in the original position of the polygon
            int_theta = wrapToPi(hand_theta-(int_theta+cur_theta));
            if int_theta*delta_theta>0 && abs(int_theta)<abs(delta_theta)
                connect = false;
                return
            end
        end
    end
end

function connect = check_translation(PG,f1,f2,cur_node_par,cur_type,next_node_par)
connect = true;
    delta_s = next_node_par(1)-cur_node_par(1);
    if delta_s==0
        return
    end
    theta_index = [4,4,4,3,3,3];
    cur_theta = cur_node_par(theta_index(cur_type));
    R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
    cur_edge = PG.get('edgeNum',cur_node_par(1));
    trans = R*PG.tangent(:,cur_edge);
    trans = trans/norm(trans);
    cur_position = PG.get('1Pos',cur_node_par(1));
    vertex = R*(PG.vertex-cur_position);
    tangent = R*PG.tangent;
    f2_rel = f2-f1;
    for i=1:size(PG.vertex,2)
        v1 = vertex(:,i);
        t1 = tangent(:,i);
        mat = [t1 -trans];
        if abs(det(mat))>1E-8 % don't check the current edge or any edges parallel to it (there will be other inersections if any exist)
            p = mat\(f2_rel-v1);
            if (p(1)>=0&&p(1)<=1)&&(p(2)<=abs(delta_s)&&p(2)*delta_s>0) % find if the edge intercepts the finger after start of movement and before end.
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
    cur_edge = R*PG.get('edgeNum',cur_node_par(1));
    cur_position = PG.get('1Pos',cur_node_par(1));
    vertex = R*(PG.vertex-cur_position);
    if any(vertex(2,:)>0)
        escape = false;
    end
end


%%
function labels = mark_possible_neighbors(possible_neighbors,nodes,index_list,type_list,labels)
theta_index = [4,4,4,3,3,3];
if ~isempty(labels)
    delete(labels)
end
hold on
for i=1:length(possible_neighbors)
    node = possible_neighbors(i);
    node_ind = index_list(node);
    node_type = type_list(node);
    node_par = nodes{node_type}(node_ind,:);
    s =node_par(1);
    theta = node_par(theta_index(node_type));
    h = node_par(theta_index(node_type)-1);
    labels(2*i-1) = plot(s,theta,'ko','markerSize',14,'lineWidth',3);
    labels(2*i) = text(s,theta+0.2,[num2str(i) ',' num2str(node) ',' num2str(h)] ,'fontSize',10);
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
%% old tries
% function possible_neighbors = search_DSnodes(PG,nodes,node_par,type,type_index_list,cont)
%     start_theta = node_par(3);
%     possible_neighbors = [];
%     %% search up
%     
%     ds_virtual = nodes{3};
%     epsi = 1E-3;
%     on_edge_nodes = find(ds_virtual(:,1)-node_par(1)<epsi);
%     thetas = ds_virtual(on_edge_nodes,4);
%     % upward pass
%     thetas_up = wrapTo2Pi(thetas-start_theta);
%     thetas_up_sorted = sort(thetas_up);
%     counter_up = 1;
%     flag = [1,1];
%     while flag(1)&&flag(2)
%         cur_node = find(thetas_up==thetas_up_sorted(counter_up)); % out of on_edge_nodes
%         cur_seg_num = on_edge_nodes(cur_node,5);
%         cur_seg_index = on_edge_nodes(cur_node,6);
%         cur_seg = cont{cur_seg_num};
%         % find which direction the contour goes forwards (2 right,1
%         % left)
%         direction_forwards = (cur_seg(1,cur_seg_index)<cur_seg(1,cur_seg_index+1))+1;
%         if flag(direction_forwards)
%             
%         
%     end
% 
% end