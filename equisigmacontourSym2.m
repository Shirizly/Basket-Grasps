function [Contours,DataSet] = equisigmacontourSym2(PG,BGSSegEnds,steps,direction,DataSet)
% function reveives a line segment of basket grasp curve, and the number of
% steps (rectangles) to explore the sigma-contour of the segment in the
% desired direction. The dataset holds the results of past computations to
% reduce computation time.
% the output is a data structure that holds information regarding the
% sigma-contours that cross the basket grasp curve. If the basket grasp
% segment received has qualitatively different contours (i.e they cross
% different rectangles along their path), the basket grasp segment is
% subsivided. The contour of the endpoints of each subdivision are given.

%note - the sigma-contours themselves are piecewise-smooth, so are
%described in cell arrays by a collection of segments. The information
%regarding each segment is only its endpoints at the boundary of the
%rectangles.

FullSet = 0;
s1st = BGSSegEnds(1,1);
s2st = BGSSegEnds(2,1);
s1end = BGSSegEnds(1,2);
s2end = BGSSegEnds(2,2);
% add function to flip st and end if needed, so that we're working from
% high sigma to low. Will be useful later for fixing special corner cases
% more easily - identify when the corner is a local maximum with regard to
% the already identified edge.
FullRangeS1 = [s1st,s1end];
FullRangeS2 = [s2st,s2end];
CurrentRangeS1 = [s1st,s1end];
CurrentRangeS2 = [s2st,s2end];
Contours = cell(0);
nextFeature = [0,0];
dictionary = DataSet{1};
Data = DataSet{2};
if isempty(Data)
    Data = cell(0);
end
while ~FullSet
    % find the contour for one end of the basket grasp set
    complete = 0;
    contour = cell(0);
    prevedge = [0;0];
    pathpoints = [];
    s10 = CurrentRangeS1(1);
    s20 = CurrentRangeS2(1);
    theta0 = 0;% Generate the contour for one end of the current set
    while ~complete %construct the contour from segments in each rectangle
        % the contour passes through
        
        if s10>PG.S(end)||s20>PG.S(end)
            disp('s error');
        end
        segment1 = [];
        if ~isempty(dictionary)
            [tf,index] = ismember(dictionary,[s10,s20,direction],'rows');
            if any(tf)
                segment1 = Data{find(index)};
            end
        end
        if isempty(segment1)
            segment1 = equisigmaSegment2(PG,s10,s20,theta0,prevedge,direction);
            dictionary(end+1,:) = [s10,s20,direction];
            Data{end+1} = segment1;
        end
        send = segment1{end};
        if size(pathpoints,2)+1==steps
            complete = 1;
            s = send;
            vertex = findvertex(PG,s);
            nextFeature(vertex(1)) = -vertex(2);
            nextFeature(3-vertex(1)) = PG.get('edgeNum',s(3-vertex(1)));
            segment1{end+1} = nextFeature;
        end
        pathpoints = [pathpoints,send.'];
        contour{end+1} = segment1;
        s10 = send(1);
        s20 = send(2);
        theta0 = segment1{1}(end);
        prevedge = segment1{2};
    end
    Found = 0;
    theta0 = 0;
    %Generate the contour for the other end of the current set
    while ~Found
        
        contour2 = cell(0);
        s10 = CurrentRangeS1(end);
        s20 = CurrentRangeS2(end);
        prevedge = [0;0];
        complete = 0;
        split = 0;
        pathpoints2 = [];
        count = 1;
        while ~complete&&~split %construct the contour from segments in each rectangle
            % the contour passes through
            %Generate the segments for both ends of the current set
            if s10>PG.S(end)||s20>PG.S(end)
                disp('s error');
            end
            segment1 = [];
        if ~isempty(dictionary)
            [tf,index] = ismember(dictionary,[s10,s20,direction],'rows');
            if any(tf)
                segment1 = Data{find(index)};
            end
        end
        if isempty(segment1)
            segment1 = equisigmaSegment2(PG,s10,s20,theta0,prevedge,direction);
            dictionary(end+1,:) = [s10,s20,direction];
            Data{end+1} = segment1;
        end
            send = segment1{end};
            if size(pathpoints2,2)+1==steps
                complete = 1;
                Found = 1;
                s = segment1{end};
                vertex = findvertex(PG,s);
                nextFeature(vertex(1)) = -vertex(2);
                nextFeature(3-vertex(1)) = PG.get('edgeNum',s(3-vertex(1)));
                segment1{end+1} = nextFeature;
            end
            pathpoints2 = [pathpoints2,send.'];
            contour2{end+1} = segment1;
            s10 = send(1);
            s20 = send(2);
            theta0 = segment1{1}(end);
            % Check if both segments reach the same rectangle edge:
            if ~sameEdge(PG,send.',pathpoints(:,count))
                Found = 0;
                split = true;
                sigmaS = sigmaSfind(PG,pathpoints(:,count),send.');
                splitpoint = splitpointFind(PG,[CurrentRangeS1;CurrentRangeS2],sigmaS);
            end
            count = count+1;
            prevedge = segment1{2};
        end
        if split
            CurrentRangeS1 = [CurrentRangeS1(1),splitpoint(1)];
            CurrentRangeS2 = [CurrentRangeS2(1),splitpoint(2)];
        end
    end
    Contours{end+1} = {contour,contour2,CurrentRangeS1,CurrentRangeS2};
    if (CurrentRangeS1(end) == FullRangeS1(end))&&(CurrentRangeS2(end) == FullRangeS2(end))
        FullSet = 1;
    else
        CurrentRangeS1 = [CurrentRangeS1(end)+2E-4*sign(FullRangeS1(end)-CurrentRangeS1(end)),FullRangeS1(end)];
        CurrentRangeS2 = [CurrentRangeS2(end)+2E-4*sign(FullRangeS2(end)-CurrentRangeS2(end)),FullRangeS2(end)];
    end
%     catch me
%         disp(me)
%     end
end
DataSet = {dictionary,Data};    
end

function check = sameEdge(PG,send1,send2)
%function checks if the rectangle edge reached by the segments is the same
%for both ends of the set checked here. If yes, return 1
% If not, return 0
check = 0;
vertex1 = findvertex(PG,send1);
vertex2 = findvertex(PG,send2);
if vertex1==vertex2
    check = 1;
end
end
function splitpoint = splitpointFind(PG,BGSSegEnds,sigmaS)
%function receives a linear set of points in s1-s2 denoted by its
%endpoints, and a sigma value lying within the set
% function returns the s1-s2 point on the set with sigma==sigmaS
s10 = BGSSegEnds(1,1);
s20 = BGSSegEnds(2,1);
dir1 = sign(BGSSegEnds(1,2)-BGSSegEnds(1,1));
dir2 = sign(BGSSegEnds(2,2)-BGSSegEnds(2,1));
beta = abs((BGSSegEnds(2,2)-BGSSegEnds(2,1))/(BGSSegEnds(1,2)-BGSSegEnds(1,1)));
[p10,p20] = PG.get('2pos',s10,s20);
if beta~=inf && beta~=-inf && beta~=0
    edge1 = PG.get('edgeNum',s10);
    edge2 = PG.get('edgeNum',s20);
    t1 = PG.tangent(:,edge1);
    t2 = PG.tangent(:,edge2);
    t1 = t1/norm(t1);
    t2 = t2/norm(t2);
    syms S
    sigma = norm(p10+dir1*S*t1-(p20+dir2*beta*S*t2));
    sol = double(vpasolve(sigma==sigmaS));
    sol = sol(imag(sol)==0);
    sol = sol(sol>=0);
    sol = min(sol);
    slightbefore = 1-1E-3;%separate the set just before the corner's sigma
    splitpoint = [s10+slightbefore*dir1*sol,s20+slightbefore*dir2*beta*sol];
else
    sboth = [s10,s20];
    vertex = findvertex(PG,sboth);
    vertnum = vertex(2);
    edgenum = PG.get('edgeNum',sboth(3-vertex(1)));
    t1 = PG.tangent(:,edgenum);
    t1 = t1/norm(t1);
    dir = [dir1,dir2];
    dire = dir(3-vertex(1));
    p = [p10,p20];
    pe = p(:,3-vertex(1));
    pv = p(:,vertex(1));
    syms S
    sigma = norm(pe+dire*S*t1-pv);
    sol = double(vpasolve(sigma==sigmaS));
    sol = sol(imag(sol)==0);
    sol = sol(sol>=0);
    sol = min(sol);
    slightbefore = 1-1E-3;%separate the set just before the corner's sigma
    if beta == 0
    splitpoint = [s10+slightbefore*dir1*sol,s20+slightbefore*dir2*beta*sol];
    end
    if beta == inf || beta == -inf
    splitpoint = [s10,s20+slightbefore*dir2*sol];
    end
end
end

function sigmaS = sigmaSfind(PG,intersect1,intersect2)
% find the sigma value of the point at which the contours jump. 
% first, find the relevant rectangle corner:
sigma1 = calculateSigma(PG,intersect1);
sigma2 = calculateSigma(PG,intersect2);
%check the slope of sigma along the first edge's intersection
vertex = findvertex(PG,intersect1);
vertex1 = vertex(1);
epsi = zeros(2,1);
epsi(3-vertex1) = 1E-5;
sigma1plus = calculateSigma(PG,intersect1+epsi);
if sigma1>sigma2
    if sigma1plus<sigma1
        corner = zeros(2,1);
        corner(vertex1) = intersect1(vertex1);
        corner(3-vertex1) = PG.S(PG.VL(find(PG.S(PG.VL)>intersect1(3-vertex1),1,'first')));
    else
        corner = zeros(2,1);
        corner(vertex1) = intersect1(vertex1);
        corner(3-vertex1) = PG.S(PG.VL(find(PG.S(PG.VL)<intersect1(3-vertex1),1,'last')));
    end
else
    if sigma1plus>sigma1
        corner = zeros(2,1);
        corner(vertex1) = intersect1(vertex1);
        corner(3-vertex1) = PG.S(PG.VL(find(PG.S(PG.VL)>intersect1(3-vertex1),1,'first')));
    else
        corner = zeros(2,1);
        corner(vertex1) = intersect1(vertex1);
        corner(3-vertex1) = PG.S(PG.VL(find(PG.S(PG.VL)<intersect1(3-vertex1),1,'last')));
    end
end
cornersigma = calculateSigma(PG,corner);
% second, check whether the rectangle edge of the first intersect contains a
% extrema point along it towards the relevant corner, if it does, that is
% the split value:
vertnum = vertex(2);
edgenum = PG.get('edgeNum',intersect1(3-vertex1));
ne = PG.normal(:,edgenum);
vp = PG.vertex(:,vertnum);
ve1 = PG.vertex(:,edgenum);
J = [0 1;-1 0];
te = J*ne;
sol = [te ne]\(vp-ve1);
if sol(1)>0 && sol(1)<PG.S(PG.VL(edgenum+1))-PG.S(PG.VL(edgenum))
    srange = [intersect1(3-vertex1),corner(3-vertex1)];
    smin = PG.S(PG.VL(edgenum))+sol(1);
    if smin < max(srange) && smin > min(srange)
        sigmaS = calculateSigma(PG,[intersect1(vertex1),smin]);
        return
    end
end
% Third, check whether the rectangle edge of the second intersect contains a
% extrema point along it towards the relevant corner, if it does, that is
% the split value:
vertex = findvertex(PG,intersect2);
vertex2 = vertex(1);
vertnum = vertex(2);
edgenum = PG.get('edgeNum',intersect2(3-vertex2));
ne = PG.normal(:,edgenum);
vp = PG.vertex(:,vertnum);
ve1 = PG.vertex(:,edgenum);
J = [0 1;-1 0];
te = J*ne;
sol = [te ne]\(vp-ve1);
if sol(1)>0 && sol(1)<PG.S(PG.VL(edgenum+1))-PG.S(PG.VL(edgenum))
    srange = [intersect2(3-vertex2),corner(3-vertex2)];
    smin = PG.S(PG.VL(edgenum))+sol(1);
    if smin < max(srange) && smin > min(srange)
        sigmaS = calculateSigma(PG,[intersect2(vertex2),smin]);
        return
    end
end

% Fourth, if no extrema was found, return the sigma value of the relevant
% corner:
sigmaS = cornersigma;
end

function sigma = calculateSigma(PG,s)
[O1,O2] = PG.get('2pos',s(1),s(2));
sigma = norm(O2-O1);
end

function vertex = findvertex(PG,s)
%recieves 2 s values, returns which of them repr. a vertex, and which
%vertex it is. assumes it is not both.
vertex1 = find(PG.S(PG.VL)==s(1));
vertex2 = find(PG.S(PG.VL)==s(2));
if ~isempty(vertex1)
    vertex = [1,cleanvertex(vertex1,PG.nv)];
    return
end
if ~isempty(vertex2)
    vertex = [2,cleanvertex(vertex2,PG.nv)];
    return
end
vertex1 = find(abs(PG.S(PG.VL)-s(1))==min(abs((PG.S(PG.VL)-s(1)))));
vertex2 = find(abs(PG.S(PG.VL)-s(2))==min(abs((PG.S(PG.VL)-s(2)))));
if min(abs((PG.S(PG.VL)-s(2))))<min(abs((PG.S(PG.VL)-s(1))))
    vertex = [2,cleanvertex(vertex2,PG.nv)];
else
    vertex = [1,cleanvertex(vertex1,PG.nv)];
end
if min(min(abs((PG.S(PG.VL)-s(2)))),min(abs((PG.S(PG.VL)-s(1)))))>1E-6
    disp('not on rectangle edge');
end
end

function cleanedvertnum = cleanvertex(vertnum,nv)
cleanedvertnum = vertnum;
if vertnum == nv+1
    cleanedvertnum = 1;
end
end