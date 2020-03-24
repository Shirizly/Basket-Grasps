classdef Polygon
% A polygon defined by a 2Xn matrix of vertex positions
properties
    vertex = []; % list of vertices
    nv = []; % number of vertices (or edges)
    com = [0;0]; % center of mass location
    J = [0 -1;1 0]; % rotation matrix by 90 degrees CCW
    res = 5; % resolution of the boundary parametrization
    gv = [0;-1]; % gravity direction
    epsilon = 5e-3;
    S = [];
    X = [];
    VL = [];
    tangent = [];
    normal = [];
    sv = [];
    elength = [];
end
methods
    %% Main functions
    function PG = Polygon(varargin)
            switch nargin
                case 0 % Empty object
                    PG; %#ok<VUNUS>
                case 2 % Input is: vertex, com
                    PG.vertex = varargin{1}; % vertex list
                    PG.com = varargin{2}; % center of mass
                    PG.nv = size(varargin{1},2);
            end
            if ~isempty(PG.vertex)
                extended_Polygon=[PG.vertex PG.vertex(:,1)];
                PG.sv(1) = 0;
                for k=1:PG.nv
                    PG.tangent(:,k)=(extended_Polygon(:,k+1)-extended_Polygon(:,k));
                    PG.normal(:,k)=PG.J*PG.tangent(:,k)/norm(PG.tangent(:,k));
                    PG.elength(k)=norm(extended_Polygon(:,k+1)-extended_Polygon(:,k));
                    PG.sv(k+1) = PG.sv(k)+PG.elength(k);
                end
            end
    end
    
    function [PG,s_tot,X_tot,vertex_list] = findBdyVariable(PG,res)
        % FINDBDYVARIABLE computes the contact space variable and (x,y) points along
        % the object's boundary, for graphical purposes.

        PG.res = res;
        vert = PG.vertex.';
        
        scale=1/res;
        
        s_tot=[];
        X_tot=[];
        vertex_list=1;
        for k=1:size(vert,1)
            if PG.elength(k)<scale
                num_points=2;
            else
                num_points=round(PG.elength(k)/scale)+1;
            end
            s=linspace(0,PG.elength(k),num_points);
            if k~=size(vert,1)
                s=s(1:end-1);
            end
            pos=vert(k,:).'*ones(size(s))+PG.tangent(:,k)*s/PG.elength(k);
            
            s_tot=[s_tot;(PG.sv(k)+s).'];
            X_tot=[X_tot;pos.'];
            vertex_list=[vertex_list;vertex_list(end)+num_points-1];
        end
        PG.VL = vertex_list;
        PG.S = s_tot;
        PG.X = X_tot;
    end
    
    function [ height_matrix ] = single_support_height(PG,s,theta)
        % single_support_height computes a matrix containing the height of
        % the center of mass for each point on the boundary held by a
        % finger at (0,0), and each possible angle.
        %%
        height_matrix=zeros(length(s),length(theta));
        for j=1:length(s)
            com_vector = PG.com-PG.get('1Pos',(s(j)));
            for i=1:length(theta)
                R = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
                com_pos = R*com_vector;
                height_matrix(j,i)=com_pos(2);
            end
        end
    end
    
    function [ dist_matrix ] = double_support_dist(PG,s,theta,f1,f2)
        % single_support_height computes a matrix containing the height of
        % the center of mass for each point on the boundary held by a
        % finger at (0,0), and each possible angle.
        %%
        n = size(PG.vertex,2);
        distance = zeros(size(PG.vertex,2),1);
        second = f2-f1;
        dist_matrix=zeros(length(s),length(theta));
        for j=1:length(s)
            axle = PG.get('1Pos',(s(j))).';
            for i=1:length(theta)
                R = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
                vertices = R*(PG.vertex-axle);
                count_intersect = 0;
                for k = 1:n
                    v1 = vertices(:,k);
                    if k<n
                        v2 = vertices(:,k+1);
                    else
                        v2 = vertices(:,1);
                    end
                    edge = v2-v1;
                    len = norm(edge);
                    edge = edge/norm(edge);
                    nm = [-edge(2);edge(1)];
                    nm = nm/norm(nm);
                    p = [edge,nm]\(second-v1);
                    if p(1)>=0 && p(1)<=len %the nearest point in inside the edge
                        distance(k) = abs(p(2));
                    else % the nearest point is a vertex
                        p(1) = p(1)*(p(1)<0)+(p(1)-len)*(p(1)>len); %take the distance from the nearest vertex along the edge
                        distance(k) = sqrt(p(1)^2+p(2)^2); % take the euclidean distance from the nearest vertex
                    end
                    
                    % check if the point is inside the polygon
                    epsi = 1E-8;
                    ray = [1;0];
                    if abs(ray.'*nm)>epsi 
                        p = [edge,-ray]\(second-v1);
                        if p(1)>=0 && p(1)<len && p(2)>0
                            count_intersect = count_intersect + 1;
                        end
                    else % in case the matrix is singular (the edge is horizontal)
                        if (v1(2)-second(2))*(v2(2)-second(2))<0 %if it's not clear if there's not an intersection (if the value is positive there isn't an intersection)
                            count_intersect = -inf; %don't use the ray casting method
                        end
                    end 
                end
                if count_intersect>=0
                dir = (mod(count_intersect,2)==0)*1-(mod(count_intersect,2)~=0)*1;
                else %just take the sign of the distance of a neighboring point
                    if i>1
                        dir = sign(dist_matrix(j,i-1));
                    else
                        dir = sign(dist_matrix(j-1,i));
                    end
                end
                dist_matrix(j,i)=dir*min(distance);
            end
        end
    end
    %% check functions
    function flag = checkConcave(PG,vertnum)
        flag = false;
        nor1 = PG.normal(vertnum);
        if vertnum>1
            edge2 = -PG.tangent(vertnum-1);
        else
            edge2 = -PG.tangent(end);
        end
        if nor1.'*edge2<0
            flag = true;
        end
    end

    %% equal height functions
    function [intercept,s2,segment] = doubleSupportReached_translation(PG,segment,x1,x2)
        intercept = false;
        s2 = [];
    end
    function [intercept,s2,segment] = doubleSupportReached_rotation(PG,segment,x1,x2)
        intercept = false;
        s2 = [];
    end
    function flag = escape(PG,cur_spar,cur_theta) %checks whether all vertex are outside of the forbidden quadrant (1 or 2)
        flag = true;
        cur_vert = wrap(PG,find(PG.S(PG.VL)==cur_spar),1);%current vertex
        range = 1:size(PG.vertex,2);
        range(cur_vert) = [];
        R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
        axle = PG.vertex(:,cur_vert);
        epsi = 1E-3;
        for i=range
            vert = R*(PG.vertex(:,i)-axle); %current location of vertex(i), relative to axle;
            if vert(2)>epsi %is there a vertex above the finger
                flag = false;
                return
            end
        end
    end
    
    function flag = fall(PG,cur_spar,cur_theta,dir) %checks whether the object would naturally freefall at this point
        flag = true;
        cur_vert = wrap(PG,find(PG.S(PG.VL)==cur_spar),1);%current vertex
        R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
        epsi = 1E-3;
        if dir>0
            norm = PG.normal(:,cur_vert);
            edge = PG.tangent(:,cur_vert);
        else
            if cur_vert>1
                norm = PG.normal(:,cur_vert-1);
                edge = -PG.tangent(:,cur_vert-1);
            else
                norm = PG.normal(:,end);
                edge = -PG.tangent(:,end);
            end
        end
        norm = R*norm;
        edge = R*edge;
        if norm(2)>-epsi || edge(2)<epsi %is the next normal y-axis positive or is the next edge pointing downwards
            flag = false;
        end
    end
    
    function vert_num = wrap(PG,vert_num,dir) %connect the first and n+1 vertex
        if vert_num>=length(PG.VL) && dir>0
            vert_num = 1;
        else
            if vert_num <= 1 && dir<0
                vert_num = length(PG.VL);
            end
        end
    end

    %% Contour functions
    
    function [cont,graspind,sig] = GetSigmaContour(PG,Sigma,sf1,sf2)
        % receives a single grasp identified by the pair sf1,sf2
        % the inter-finger distance matrix Sigma
        % Outputs: cont - cell array of the contour line of
        % the grasp, separated into the discontinuous parts
        % graspind - the index in cont{1} of the grasp itself
        maxS = max(PG.S);
        [f1pos,f2pos] = PG.get('2Pos',sf1,sf2);
        sig = inter_finger_distance(f1pos.',f2pos.');
        figure
        [c,~] = contour(PG.S,PG.S,Sigma,[0;sig],'showtext','on','linewidth',2,'linecolor','k');
        close
        
        ind = find(c(2,:)>maxS); % remove the out of bounds values
        c(:,ind) = [];
        
        difc(1,:) = diff(c(1,:));
        difc(2,:) = diff(c(2,:));
        normdifc = sqrt(difc(1,:).*difc(1,:)+difc(2,:).*difc(2,:));
        sep = find(normdifc>0.1);
        nparts = length(sep)+1;
        sep = [0,sep,size(c,2)];
        part = cell(nparts,1);
        graspind = [];
        epsi = 0.5;
        for i =1:nparts
            part{i} = [c(1,sep(i)+1:sep(i+1));c(2,sep(i)+1:sep(i+1))];
            %     figure
            %     hold on
            %     plot(part{i}(1,:),part{i}(2,:))
            %     axis([0,maxS,0,maxS])
            if isempty(graspind)
                possiblex = find(abs(part{i}(1,:)-sf1)<epsi);
                possibley = find(abs(part{i}(2,:)-sf2)<epsi);
                possiblegrasp = intersect(possiblex,possibley);
                mindist = 2*epsi^2;
                for j=1:length(possiblegrasp)
                    dist = (part{i}(1,possiblegrasp(j))-sf1)^2+(part{i}(2,possiblegrasp(j))-sf2)^2;
                    if mindist>dist
                        mindist=dist;
                        graspind = possiblegrasp(j);
                        partord = i;
                    end
                end
            end
        end
        % Identifying the particular contour for the grasp
        partsv = 1:nparts;
        partsv(partord) = [];
        epsi = 0.01;
        first = part{partord}(:,1); %last point in the last part of the contour;
        flag = 1;
        while flag
            last = part{partord(end)}(:,end); %last point in the last part of the contour;
            if last(1) == 0
                last(1) = maxS;
                
            else
                if last(1) >= maxS-epsi
                    last(1) = 0;
                end
            end
            if last(2) == 0
                last(2) = maxS;
            else
                if last(2) >= maxS-epsi
                    last(2) = 0;
                end
            end
            mindist = epsi;
            if norm(last-first)<epsi
                break
            end
            nextpart = [];
            for i=partsv
                dist = norm(last-part{i}(:,1));
                if dist<mindist
                    mindist = dist;
                    nextpart = i;
                end
                dist = norm(last-part{i}(:,end));
                if dist<mindist
                    mindist = dist;
                    nextpart = -i;
                end
            end
            if nextpart<0
                part{-nextpart} = flip(part{-nextpart},2);
                nextpart = -nextpart;
            end
            partord(end+1) = nextpart;
            nextind = find(partsv==nextpart);
            partsv(nextind) = [];
            if isempty(nextpart)
                error('Error!');
            end
        end
        cont = part(partord.');
    end
    
    function [cont] = GetSigmaContours(PG,Sigma,sig)
        % receives a specific inter finger distance (sig)
        % the inter-finger distance matrix Sigma
        % Outputs: cont - cell array of the contour lines of
        % the grasp, separated into the discontinuous parts
        maxS = max(PG.S);
        figure
        [c,~] = contour(PG.S,PG.S,Sigma,[0;sig],'showtext','on','linewidth',2,'linecolor','k');
        close
        
        ind = find(c(2,:)>maxS); % remove the out of bounds values
        c(:,ind) = [];
        
        difc(1,:) = diff(c(1,:));
        difc(2,:) = diff(c(2,:));
        normdifc = sqrt(difc(1,:).*difc(1,:)+difc(2,:).*difc(2,:));
        sep = find(normdifc>0.1);
        nparts = length(sep)+1;
        sep = [0,sep,size(c,2)];
        part = cell(nparts,1);
        epsi = 0.5;
        for i =1:nparts
            part{i} = [c(1,sep(i)+1:sep(i+1));c(2,sep(i)+1:sep(i+1))];
        end
        % Identifying the particular contour for the grasp
%         partsv = 1:nparts;
%         partsv(partord) = [];
%         epsi = 0.01;
%         first = part{partord}(:,1); %last point in the last part of the contour;
%         flag = 1;
%         while flag
%             last = part{partord(end)}(:,end); %last point in the last part of the contour;
%             if last(1) == 0
%                 last(1) = maxS;
%                 
%             else
%                 if last(1) >= maxS-epsi
%                     last(1) = 0;
%                 end
%             end
%             if last(2) == 0
%                 last(2) = maxS;
%             else
%                 if last(2) >= maxS-epsi
%                     last(2) = 0;
%                 end
%             end
%             mindist = epsi;
%             if norm(last-first)<epsi
%                 break
%             end
%             nextpart = [];
%             for i=partsv
%                 dist = norm(last-part{i}(:,1));
%                 if dist<mindist
%                     mindist = dist;
%                     nextpart = i;
%                 end
%                 dist = norm(last-part{i}(:,end));
%                 if dist<mindist
%                     mindist = dist;
%                     nextpart = -i;
%                 end
%             end
%             if nextpart<0
%                 part{-nextpart} = flip(part{-nextpart},2);
%                 nextpart = -nextpart;
%             end
%             partord(end+1) = nextpart;
%             nextind = find(partsv==nextpart);
%             partsv(nextind) = [];
%             if isempty(nextpart)
%                 error('Error!');
%             end
%         end
%         cont = part(partord.');
        cont = part;
    end
    
    function [s1,s2] = GetSigPoint(PG,Sigma,sig,rect)
        figure
        [c,~] = contour(PG.S,PG.S,Sigma,[0;sig],'showtext','on','linewidth',2,'linecolor','k');
        close
        ind1 = intersect(find(c(1,:)>rect(1)),find(c(1,:)<rect(2)));
        ind2 = intersect(find(c(2,:)>rect(3)),find(c(2,:)<rect(4)));
        ind = intersect(ind1,ind2);
        s1 = c(1,ind(1));
        s2 = c(2,ind(1));
    end
    
    %% Equilibrium identification functions
    function [s,h,theta] = SSEqcheck(PG) % find saddle points on edges
        s = -ones(1,size(PG.VL,1)-1);
        h =  -ones(1,size(PG.VL,1)-1);
        theta = -ones(1,size(PG.VL,1)-1);
        for i =1:size(PG.VL,1)-1 % run over the edges
            [v1,v2,nor] = PG.get('edge',i);
            edge = (v2-v1);
            p = [nor,edge]\(PG.com-v1);
            if p(2)>=0 && p(2)<=1
                s1 = PG.S(PG.VL(i));
                s2 = PG.S(PG.VL(i+1));
                if p(2)>1-norm(edge)/PG.res %in case the point is very close to the vertex
                    p(2) = 1-norm(edge)/PG.res; % take one discrete point before the vertex instead;
                end
                s(i) = s1+p(2)*(s2-s1);

                s(i) = PG.S(find(PG.S>=s(i),1,'first'));
                h(i) = p(1);
                theta(i) = wrapToPi(pi()/2-atan2(nor(2),nor(1)));
            end
        end
        remove = find(s==-1);
        s(remove) = [];
        h(remove) = [];
        theta(remove) = [];
        s = [s,s];
        h = [h,-h];
        theta = [theta,wrapToPi(theta+pi())];
    end
    function [s,h,t] = HookEqcheck(PG) % find single-support hook grasps - concave vertices with a stable grasp
        s = -ones(1,size(PG.VL,1)-1);
        h = -ones(1,size(PG.VL,1)-1);
        t = -ones(1,size(PG.VL,1)-1);
        for i =1:size(PG.VL,1)-1 % run over the vertices - check whether both normals are y-axis positive and x-axis different, after rotating to put the com directly beneath the vertex
            % first check if the vertex is concave
            if checkConcave(PG,i)
                [~,~,nor1] = PG.get('edge',i);
                if i>1
                    [~,~,nor2] = PG.get('edge',i-1);
                else
                    [~,~,nor2] = PG.get('edge',size(PG.vertex,2));
                end
                axle = PG.vertex(:,i);
                com_vector = (PG.com(:)-axle);
                cur_theta = -pi()/2-atan2(com_vector(2),com_vector(1));
                R = [cos(cur_theta) -sin(cur_theta); sin(cur_theta) cos(cur_theta)];
                com_vector = R*(PG.com(:)-axle);
                p = [R*nor1,R*nor2]\(com_vector);
                if (p(1)*p(2)>0 && p(1)<0)
                    s(i) = PG.S(PG.VL(i));
                    h(i) = com_vector(2);
                    t(i) = cur_theta;
                end
            end
        end
        remove = find(s==-1);
        s(remove) = [];
        h(remove) = [];
        t(remove) = [];
    end
    
    function [s,h,t] = SSMincheck(PG) % find minimal height single support grasps at vertices
        s = PG.S(PG.VL(1:end-1)).';
        h = zeros(1,size(PG.VL,1)-1);
        t = zeros(1,size(PG.VL,1)-1);
        for i =1:size(PG.VL,1)-1 % run over the vertices - find the necessary rotation to bring com h to minimum
            axle = PG.vertex(:,i);
            com_vector = (PG.com(:)-axle);
            rotate_theta = -pi()/2-atan2(com_vector(2),com_vector(1));
            R = [cos(rotate_theta) -sin(rotate_theta); sin(rotate_theta) cos(rotate_theta)];
            com_vector = R*(PG.com(:)-axle);
            h(i) = com_vector(2);
            t(i) = wrapToPi(rotate_theta);
        end
    end
    
    function saddlepoint = EV2check(PG,vertnum,edgenum)
        % get all the geometric parameters
        [v1,v2,edgenormal] = PG.get('edge',edgenum);
        if vertnum==length(PG.VL)
            [v,~,vertnormal2] = PG.get('edge',1);
        else
            [v,~,vertnormal2] = PG.get('edge',vertnum);
        end
        if vertnum==1
            [~,~,vertnormal1] = PG.get('edge',length(PG.VL)-1);
        else
            [~,~,vertnormal1] = PG.get('edge',vertnum-1);
        end
        % check if the edge normal can be countered by the vertex
        p = [vertnormal1 vertnormal2]\(-edgenormal);
        cond1 = all(p>0);
        
        % check if the edge normal connects the edge and the vertex
        p = [edgenormal v2-v1]\(v-v1);
        cond2 = and(p(1)>0,and(p(2)>=0,p(2)<=1));
        if and(cond1,cond2)
            saddlepoint = [PG.S(PG.VL(vertnum));PG.S(PG.VL(edgenum))+p(2)*(PG.S(PG.VL(edgenum+1))-PG.S(PG.VL(edgenum)))];
        else
            saddlepoint = [];
        end
        
    end
    
    function EQlist = Rectcheck(PG,rect)
        edges = [rect(3) rect(3) rect(1) rect(1)];
        count = 0;
        EQlist = [];
        for i=1:4
            vertnum = find(PG.S(PG.VL)==rect(i));
            edgenum = find(PG.S(PG.VL)==edges(i));
            saddlepoint = PG.EV2check(PG.S,PG.VL,vertnum,edgenum);
            if ~isempty(saddlepoint)
                EQlist(:,count+1) = saddlepoint;
                count = count+1;
                if i>2 % flipping the saddle point for when the vertex is s2
                    EQlist(:,end)=fliplr(EQlist(:,end));
                end
            end
            
        end
    end
    
    function [flag,type] = EEGraspCheck(PG,p1,p2,n1,n2,varargin)
        flag = 0;
        gravity = PG.gv;
        switch nargin
            case 5 %default
            case 6 %rotated object
                theta = varargin{1};
                gravity = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        end
        N = [n1 n2];
        if abs(det(N))<PG.epsilon % parallel edges
            cond2 = and(all(abs(N.'*PG.J*gravity)<(PG.epsilon/100)),any(N.'*gravity<0)); % normal directions all vertical, at least one pointing up
            cond11 = ((PG.com-p1).'*PG.J*n1)*((PG.com-p2).'*PG.J*n2)<0; % forces generate opposite torques
            cond12 = 1;
            if any(N(2,:)<0) % checks that a hammer handle case is properly handeled
                cond13 = and(N(2,1)<0,(abs((PG.com-p1).'*PG.J*n1))>abs(((PG.com-p2).'*PG.J*n2)));
                cond14 = and(N(2,2)<0,(abs((PG.com-p1).'*PG.J*n1))<abs(((PG.com-p2).'*PG.J*n2)));
                cond12 = or(cond14,cond13);
            end
            cond1 = and(cond12,cond11);
            type = 2; %parallel edges
        else % non parallel edges
            lamda = [n1 n2]\(-gravity);
            cond2 = all(lamda>0); % contact forces (+gravity) positively span the origin
            grmat = [n1 n2 gravity;
                -(PG.com-p1).'*PG.J*n1 -(PG.com-p2).'*PG.J*n2 0];
            cond1 = abs(det(grmat))<PG.epsilon; % grasp matrix is singular
            type = 1; %non parallel
        end
        if and(cond1,cond2)
            flag = 1;
        end
    end
    
    function flag = EVGraspCheck(PG,p1,p2,n11,n12,n2)
        flag = 0;
        gravity = PG.gv;
        G = [n11 n12 n2;
            -(PG.com-p1).'*PG.J*n11 -(PG.com-p1).'*PG.J*n12 -(PG.com-p2).'*PG.J*n2];
        if abs(det(G))<PG.epsilon
            cond1 = 0;
        else
            lamda = G\[-gravity;0];
            cond1 = and(all(lamda>=0),lamda(1)~=lamda(2));
        end
        if cond1
        flag = 1;
        end
    end
    
    function flag = VVGraspCheck(PG,p1,p2,n11,n12,n21,n22)
        Aeq = [1 0 1 0 0;
            0 1 0 1 0;
            -r1(2) r1(1) -r2(2) r2(1) -1];
        beq = [0; 1; 0];
        A = [t1(1)-mu(1)*n1(1), t1(2)-mu(1)*n1(2), 0, 0, 0;
            -t1(1)-mu(1)*n1(1), -t1(2)-mu(1)*n1(2), 0, 0, 0;
            0, 0, t2(1)-mu(2)*n2(1), t2(2)-mu(1)*n2(2), 0;
            0, 0, -t2(1)-mu(2)*n2(1), -t2(2)-mu(1)*n2(2), 0;];
        A = [n11,n12,0,0,0];
        b = [0; 0; 0; 0];
        C = [0;0;0;0;1];
        %solving problem
        [Xmin,fval1,ef1] = linprog(C,A,b,Aeq,beq);
        [Xmax,fval2,ef2] = linprog(-C,A,b,Aeq,beq);
        
        
    end

    function [s1,s2,linesegs] = Eqcheck(PG)
        s1 = 0;
        s2 = 0;
        se = PG.S;
        se(PG.VL) = [];
        Xe = PG.X;
        Xe(PG.VL,:) = [];
        startline = [];
        linesegs = zeros(1,4);
        for i=2:length(se) %edge-edge check
            for j = 1:i
                p1 = Xe(i,:).';
                p2 = Xe(j,:).';
                n1 = PG.get('normal',se(i));
                n2 = PG.get('normal',se(j));
%                 if and(se(i)>18.123,se(j)>4.123)
%                     1
%                 end
                [flag,type] = PG.EEGraspCheck(p1,p2,n1,n2);
                if and(flag,type==2)
                    if isempty(startline)
                        startline = [se(i),se(j)];
                        endline = [se(i),se(j)];
                    else
                        endline = [se(i),se(j)];
                    end
                end
                if flag&&type==1
                    s1(end+1) = se(i);
                    s2(end+1) = se(j);
                    s1(end+1) = se(j);
                    s2(end+1) = se(i);
                end
                if and(~flag,~isempty(startline)) %end of a line segment
                    linesegs(end+1,:) = [startline,endline];
                    startline = [];
                end
            end
        end
        
        for i=1:length(PG.VL)-1 %vertex-edge check
            for j=1:length(se)
                p1 = PG.X(PG.VL(i),:).';
                p2 = Xe(j,:).';
                n11 = PG.get('normal',PG.S(PG.VL(i)));
                n12 = PG.get('normal',PG.S(PG.VL(i)+1));
                n2 = PG.get('normal',se(j));
                if PG.EVGraspCheck(p1,p2,n11,n12,n2)
                    s1(end+1) = PG.S(PG.VL(i));
                    s2(end+1) = se(j);
                    s1(end+1) = se(j);
                    s2(end+1) = PG.S(PG.VL(i));
                end 
            end
        end
        
%         for i=2:length(PG.VL)-1 %vertex-vertex check
%             for j=1:i-1
%                 p1 = X(PG.VL(i),:).';
%                 p2 = X(PG.VL(j),:).';
%                 n11 = PG.get('normal',s(PG.VL(i)));
%                 n12 = PG.get('normal',s(PG.VL(i)+1));
%                 n21 = PG.get('normal',s(PG.VL(j)));
%                 n22 = PG.get('normal',s(PG.VL(j)+1));
%                 if PG.VVGraspCheck(p1,p2,n11,n12,n21,n22)
%                     s1(end+1) = s(PG.VL(i));
%                     s2(end+1) = s(PG.VL(j));
%                     s1(end+1) = s(PG.VL(j));
%                     s2(end+1) = s(PG.VL(i));
%                 end 
%             end
%         end
        s1 = s1(2:end);
        s2 = s2(2:end);
        linesegs = linesegs(2:end,:);
    end

    function status = Statuscheck(PG,s1,s2,S,X,VL)
        status = zeros(size(s1));
        for i=1:length(s1)
            sind1 = find(PG.S==s1(i));
            sind2 = find(PG.S==s2(i));
            vert1 = find(PG.VL==sind1);
            vert2 = find(PG.VL==sind2);
            if or(~isempty(vert1),~isempty(vert2)) %check if V-V or V-E case
                status(i) = 0;
                if isempty(vert2)
                    n1 = PG.get('normal',s1(i));
                    n2 = PG.get('normal',s1(i)+PG.epsilon);
                    p1 = X(sind1-1,:).';
                    p2 = X(sind1+1,:).';
                    lamda = [n1 -n2]\(p2-p1);
                    status(i) = sign(lamda(1));
                else
                    if isempty(vert1)
                    n1 = PG.get('normal',s2(i));
                    n2 = PG.get('normal',s2(i)+PG.epsilon);
                    p1 = X(sind2-1,:).';
                    p2 = X(sind2+1,:).';
                    lamda = [n1 -n2]\(p2-p1);
                    status(i) = sign(lamda(1));
                    end
                end    
            else
                p1 = X(sind1,:).';
                p2 = X(sind2,:).';
                n1 = PG.get('normal',s1(i));
                n2 = PG.get('normal',s2(i));
                N = [n1 n2];
                if det(N) == 0
                    status(i) = 0;
                else
                    flamda = N\(-PG.gv);
                    ro = [n1 -n2]\(p2-p1);
                    p = p1+ro(1)*n1;
                    status(i) = sign(-flamda.'*ro+(PG.com(2)-p(2)));
                end
            end
        end
    end

    
end
end



%%

%     function flag = VVGraspCheck(PG,p1,p2,n11,n12,n21,n22) %theory based
%     try
%         flag = 0;
%         cond1=0;
%         cond2=0;
%         N = [n11 n12 n21 n22];
%         v = [p1 p2];
%         p = zeros(2,4);
%         check = 0;
%         for i=1:2
%             p(:,i*2-1) = [N(:,i*2-1),-PG.gv]\(PG.com-v(:,i)); % finds parametrization for projection of the vertex on gravity
%             p(:,i*2) = [N(:,i*2),-PG.gv]\(PG.com-v(:,i));
%             range(:,i) = sort(p(2,i*2-1:i*2)).';
%             if p(1,i*2-1)*p(1,i*2)<0
%                 range(:,i) = range([2,1],i);
%                 check = check+2^i;
%             end
%         end
%         switch check
%             case 0 % works only for convex bodies ATM
%                 interrange = [max(range(1,:)) min(range(2,:))];
%                 if interrange(1)<=interrange(2)
%                     highestjointpoint = min(range(2,:));
%                     cond1 = or(p1(2)<highestjointpoint,p2(2)<highestjointpoint);
%                     cond2 = (PG.com(1)-p1(1))*(PG.com(1)-p2(1))<0;
%                 end                    
%             case 2
% %                 interrange = [max(range(1,:)) min(range(2,:))];
%             case 4
%                 
%             case 6
%                 
%         end
%         if and(cond1,cond2)
%             flag = 1;
%         end
%     end


%%
%     function fig = drawPolygonSseparate(PG,s1,s2,S,X)
%         % draw the poygon and add the each possible grasp
% %         fig = PG.drawPolygon();
% %         hold on
%         dist = 0;
%         preX1 = [0,0];
%         preX2 = [0,0];
%         colorv = [0,0,1]; %rand(1000,3);
%         curcol = 1;
%         c1 = 0;
%         for i=1:2:length(s1)
%             sind1 = find(S==s1(i));
%             sind2 = find(S==s2(i));
%             dist = sqrt(norm(X(sind1,:)-preX1)^2+norm(X(sind2,:)-preX2)^2);
%             if dist>PG.epsilon*500
%                 fig = PG.drawPolygon();
%                 hold on
% %                 curcol = curcol+1;
%             end
%             size1 = 4;
%             size2 = 4;
%             if all(preX1 == X(sind1,:))&&c1>4
%                 size1 = 20;
%                 c1 = c1+1;
%             else
%                 c1 = 0;
%             end
% %             plot(X(sind1,1),X(sind1,2),'.b');
% %             plot(X(sind2,1),X(sind2,2),'.b'); 
%             plot(X(sind1,1),X(sind1,2),'.','Color',colorv(curcol,:),'MarkerSize',size1);
%             plot(X(sind2,1),X(sind2,2),'.','Color',colorv(curcol,:),'MarkerSize',size2); 
%             preX1 = X(sind1,:);
%             preX2 = X(sind2,:);
%         end
%         hold off 
%     end
%     
%% equal height functions
%     function new_segment = starting_seg(PG,cur_edge,cur_theta,h,s_start,s_end,spar,dir)
%     % generate a segment of equal-height path starting from the saddlepoint
%     % path begins at one of 4 directions, according to dir, and lasts from
%     % astart to aend.
%     switch dir
%         case {1,-1} % theta changing path
%             s = spar;
%             t = (pi() - 2*atan2(h,-s));
%             %         y = h*cos(t)-a.*sin(t);
%             x = -s.*cos(t)-h*sin(t);
%         case {2,-2} % theta = 0 path
%             s = spar;
%             t = 0*s;
%             x = -s;
%     end
%     t = subs(t,spar,spar-s_start);
%     t = t + cur_theta;
%     x = subs(x,spar,spar-s_start);
%     new_segment = {t,x,s_start,s_end};
%     end
%     
%     function new_segment = middle_seg(PG,cur_edge,cur_theta,h,X0,s_start,s_end,spar,dir)
%     % generate a segment of equal-height path starting from the saddlepoint
%     % path begins at one of 4 directions, according to dir, and lasts from
%     % astart to aend.
%     [~,~,nm] = PG.get('edge',cur_edge);
%     edge = -PG.J*nm;
%     t0 = atan2(dir*edge(2),dir*edge(1))+cur_theta;
%     phi0 = atan2(h,X0);
%     df = phi0-t0;
%     r = (h^2+X0^2)^0.5;
%     b = r*cos(df);
%     disc = b^2-X0^2;
%     s = spar;
%     
%     stop_cond1 = b-sqrt(disc); %stop point on the edge
%     
%     stop_cond2 = dir*(s_end-s_start);% stop point is the next vertex
%     if disc<0
%         stop_cond1 = stop_cond2;
%     end
%     stop_cond = min(stop_cond1,stop_cond2);
%     c = (r^2+s^2-2*r*s*cos(df))^0.5;
%     beta = atan2(r*sin(df),r*cos(df)-s);
%     tgal = asin(h/c)-beta;
%     y1 = r*sin(tgal+df)-s*sin(tgal);
%     x1 = -s*cos(tgal)+r*cos(tgal+df);
%     t1 = tgal-(t0-cur_theta);
%     if stop_cond<stop_cond2
%         tstart = double(subs(t1,spar,stop_cond)); %definitions for the second part
%         X0 = real(double(subs(x1,spar,stop_cond)));
%         t0 = real(atan2(edge(2),edge(1))+tstart);
%         phi0 = atan2(h,X0);
%         df = phi0-t0;
%         r = (h^2+X0^2)^0.5;
%         s = -spar;
%         
%         c = (r^2+s^2-2*r*s*cos(df))^0.5;
%         beta = atan2(r*sin(df),r*cos(df)-s);
%         tgal = pi()-asin(h/c)-beta;
%         y2 = r*sin(tgal+df)-s*sin(tgal);
%         x2 = -s*cos(tgal)+r*cos(tgal+df);
%         t2 = tgal-(t0-tstart);
%         
%         t1 = subs(t1,spar,spar-s_start);
%         x1 = subs(x1,spar,spar-s_start);
%         t2 = subs(t2,spar,spar-s_start);
%         x2 = subs(x2,spar,spar-s_start);
%         
%         new_segment = {t1,x1,s_start,s_start+stop_cond;
%             t2,x2,s_start+stop_cond,s_start};
%     else
%     
%     t1 = subs(t1,spar,spar-s_start);
%     x1 = subs(x1,spar,spar-s_start);
%     new_segment = {t1,x1,s_start,s_start+dir*stop_cond};
%     end
%     end
%   
%  function path_list = getEqualPaths(PG,S,X,VL)
%         syms spar
%         [sss,ssh]= PG.SSEqcheck(S,X,PG.VL);
%         for i=1:length(sss)
%             
%             start_edge = find(PG.S(PG.VL)>sss(i),1,'first') - 1;
%             h = ssh(i);
%             [~,~,nm] = PG.get('edge',start_edge);
%             start_theta = pi()/2-atan2(nm(2),nm(1));
%             dirv = [1,-1,2,-2];
%             paths = cell(4,1);
%             for j = 1:4
%                 cur_edge = start_edge;
%                 cur_theta = start_theta;
%                 flag = true;
%                 path = cell(0);
%                 dir = dirv(j);
%                 while flag
%                     if numel(path) == 0 %% new path
%                         s_start = sss(i);
%                         s_end = (dir>0)*PG.S(PG.VL(cur_edge+1))+(dir<0)*PG.S(PG.VL(cur_edge));
%                         new_segment = PG.starting_seg(cur_edge,cur_theta,h,s_start,s_end,spar,dir);
%                         path{1} = new_segment;
%                         cur_edge = cur_edge+sign(dir);
%                         cur_edge = cur_edge - (cur_edge>PG.nv)*PG.nv + (cur_edge<1)*PG.nv;
%                         cur_theta = double(subs(new_segment{1},spar,new_segment{4}));
%                         dir = sign(dir);
%                     
%                     else %% middle of the path
%                         s_start = path{end}{4};
%                         if and(s_start==PG.S(PG.VL(end)),dir>0)
%                             s_start=0;
%                         else
%                             if and(s_start==0,dir<0)
%                                 s_start=PG.S(PG.VL(end));
%                             end
%                         end
%                            
%                         s_end = (dir>0)*PG.S(PG.VL(cur_edge+1))+(dir<0)*PG.S(PG.VL(cur_edge));
%                         X0 = double(subs(path{end}{2},spar,path{end}{4}));
%                         [new_segment] = PG.middle_seg(cur_edge,cur_theta,h,X0,s_start,s_end,spar,dir);
%                         if numel(new_segment)==4
%                             path{numel(path)+1} = new_segment;
%                         else
%                             path{numel(path)+1} = new_segment(1,:);
%                             path{numel(path)+1} = new_segment(2,:);
%                             dir = -dir;
%                         end
%                         flag = false;
%                     end
%                     relevantss = find(ssh==ssh(i)); %saddle points of the same height
%                     if s_start>s_end %check if the path crossed a saddlepoint
%                         flag1 = isempty(find(and(sss(relevantss)>s_start,sss(relevantss)<s_end)));
%                     else
%                         flag1 = isempty(find(and(sss(relevantss)<s_start,sss(relevantss)>s_end)));
%                     end
%                     flag = and(flag,flag1);
%                 end
%                 paths{j} = path;
%             end
%             
%             path_list{i} = {paths,ssh(i)};
%         end
%     end