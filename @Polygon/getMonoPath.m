function [path_construct] = getMonoPath(PG,s1,dir,f1,f2)
% receives a single-support grasp (parametrized by the boundary
% point s1), and the position of the 2 fingers (f1,f2). dir is the
% direction over the boundary, positive means in general COM goes
% left, negative means going right.
syms spar; % s-parameter for the path segments
syms d_theta; % parameter for rotation segments
epsi = 1E-3;
%% define the starting parameters
start_edge = PG.get('edgeNum',s1);
edge = dir*PG.tangent(:,start_edge);
edge_angle = wrapTo2Pi(atan2(edge(2),edge(1)));
start_theta = (dir>0)*(2*pi()-edge_angle)+(dir<0)*(pi()-edge_angle);
cur_theta = start_theta;
cur_edge = start_edge;
cur_spar = s1;
state = 1; % parameter specifying which case for the loop, 1 = saddlepoint
%% perform the loop
flag = 1;
path = cell(0);
spar_end = 0; %these are necessary for end-of-loop checks
theta_end = 0;
intercept_edge = 0;
stop_cond = 0;
while flag
    COM_pos = PG.get('comLocationSS',f1,cur_spar,cur_theta);
    cur_x = COM_pos(1);
    cur_y = COM_pos(2);
    switch state %here we define the movement until a stopping
        % criteria is reached (other than a 2-finger grasp)
        case 1 % edge is horizontal
            % move according to direction until next vertex
            x = matlabFunction(cur_x+dir*spar);
            y = cur_y;
            theta = cur_theta;
            spar_start = cur_spar;
            end_vertex = cur_edge+(dir+1)/2;
            end_vertex = PG.wrap(end_vertex,dir);
            spar_end = PG.S(PG.VL(end_vertex));
            segment = {x,y,theta,spar_start,spar_end};
            [intercept,s2,segment] = doubleSupportReached_translation(PG,segment,f1,f2);
        case 2 % edge that allows downwards movement
            % move along edge until next vertex
            edge = PG.tangent(:,cur_edge);
            x = matlabFunction(cur_x+dir*edge(1)*spar);
            y = matlabFunction(cur_y+dir*edge(2)*spar);
            theta = cur_theta;
            spar_start = cur_spar;
            end_vertex = cur_edge+(dir+1)/2;
            end_vertex = PG.wrap(end_vertex,dir);
            spar_end = PG.S(PG.VL(end_vertex));
            segment = {x,y,theta,spar_start,spar_end};
            [intercept,s2,segment] = doubleSupportReached_translation(PG,segment,f1,f2);
        case 3 % next edge requires rotation (points down rather than up)
            % rotate body until stopping criteria: the rotation is
            % always the same sign as dir
            radius_to_COM = PG.get('comLocationSS',f1,cur_spar,cur_theta)-f1;
            rotation = [cos(d_theta) -sin(d_theta);sin(d_theta) cos(d_theta)]; %rotation matrix for the body's coordinate system
            pos = f1+rotation*radius_to_COM;
            x = matlabFunction(pos(1));
            y = matlabFunction(pos(2));
            theta = matlabFunction(d_theta);
            theta_start = cur_theta;
            %stopping criteria:
            %1. next edge is horizontal:
            edge = double(dir*subs(rotation,d_theta,cur_theta)*PG.tangent(:,cur_edge));
            edge_angle = wrapTo2Pi(atan2(edge(2),edge(1)));
            d_theta_required = (dir>0)*(2*pi()-edge_angle)+(dir<0)*(pi()-edge_angle);
            %2. rotation brings center of mass past a minimum -
            %single-support basket grasp (hook-hang):
            angle_to_COM = wrapTo2Pi(atan2(radius_to_COM(2),(radius_to_COM(1))));
            d_theta_to_equlilbrium = (dir>0)*(3*pi()/2-angle_to_COM)+(dir<0)*(-pi()/2-wrapToPi(angle_to_COM));
            % choose the lesser rotation as the limit
            theta_end = cur_theta+dir*min(abs(d_theta_required),abs(dir*d_theta_to_equlilbrium)); %take the shorter angle difference
            segment = {x,y,theta,theta_start,theta_end,dir};
            if dir*min(abs(d_theta_required),abs(dir*d_theta_to_equlilbrium))==d_theta_to_equlilbrium
                % if reached a hook grasp, stop and mark it;
                flag = 0;
                stop_cond = 'hook';
                com_pos = PG.get('comLocationSS',f1,s1,theta_end);
                info = com_pos(2);
            end
            [intercept,s2,segment] = doubleSupportReached_rotation(PG,segment,f1,f2);
        case 4 %Free fall - add later
            rotation = [cos(cur_theta) -sin(cur_theta);sin(cur_theta) cos(cur_theta)]; %rotation matrix for the body's coordinate system
            axle = PG.get('1Pos',cur_spar);
            vertices = rotation*(PG.vertex-repmat(axle.',1,size(PG.vertex,2)));
            fall_dist = max(vertices(2,:));
            [interceptss,spar_end] = freefall_intercept(PG,fall_dist,vertices);
            if interceptss == false
                flag = 0;
                stop_cond = 'escape';
                info = [];
                intercept_edge = 1; % meaningless, prevents the code from breaking later
                spar_end = 0; %same
            else
                [dir,intercept_edge] = checkdirdown(PG,spar_end,cur_theta);
            end
            segment = {cur_spar,spar_end,cur_theta};
            [intercept,s2,segment] = doubleSupportReached_translation(PG,segment,f1,f2);
            
    end
    if intercept %check if reached double support
        flag = 0;
        stop_cond = 'double';
        info = s2;
    end
    % prepare next run's other parameters:
    cur_edge = (state<3)*wrap(PG,cur_edge+dir,dir) + (state==3)*cur_edge + (state==4)*intercept_edge;
    cur_spar = (state~=3)*spar_end + (state==3)*cur_spar;
    cur_theta = (state==3)*theta_end + (state~=3)*cur_theta;
    % check next state - whether the next edge requires rotation
    state = 2; 
    rotation = [cos(cur_theta) -sin(cur_theta);sin(cur_theta) cos(cur_theta)]; %rotation matrix for the body's coordinate system
    next_edge = dir*rotation*PG.tangent(:,cur_edge);
    if next_edge(2)<-epsi
        state = 3;
    end
    if PG.fall(cur_spar,cur_theta,dir)
        state = 4;
    end
    path{numel(path)+1} = segment;
    if strcmp(stop_cond,'escape')
        path{end} = [];
    end
        
end
path_construct = {path,dir,stop_cond,info};
end

function [intercept,s1] = freefall_intercept(PG,fall_dist,vertices)
    n = size(vertices,2);
    intercept = false;
    s1 = [];
    fall = inf;
    for i = 1:n
        if i < n
            edge = vertices(:,i+1)-vertices(:,i);
        else
            edge = vertices(:,1)-vertices(:,n);
        end
            hit = [edge [0;-1]]\(-vertices(:,i)); %solving v+a*e+b*(g)=0
            %hit(1) is the distance along the edge
            %hit(2) is the fall distance to intercept
            if hit(2)>0 && hit(1)>=0 && hit(1)<=1
                fall = min(fall,hit(2));
                if fall == hit(2)
                    vspar = PG.S(PG.VL(i));
                    s1 = vspar + hit(1)*norm(edge);
                end
            end
    end
    if fall<fall_dist
        intercept = true;
    end
end
function [dir,edge_num] = checkdirdown(PG,s1,cur_theta)
    edge_num = PG.get('edgeNum',s1);
    [~,~,edge] = PG.get('edge',edge_num);
    rotation = [cos(cur_theta) -sin(cur_theta);sin(cur_theta) cos(cur_theta)]; %rotation matrix for the body's coordinate system
    edge = rotation*edge;
    dir = sign(edge(2));
end
