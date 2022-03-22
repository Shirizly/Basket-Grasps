function EqCurves = EqCurveFinder(PG)
EqCurves = cell(0);
numvert = PG.nv;
%% find edge-edge equilibrium curves:
for i=1:numvert
    normal1 = PG.normal(:,i);
    tangent1 = PG.tangent(:,i);
    p1 = PG.vertex(:,i);
    s1 = PG.S(PG.VL(i));
    for j=i+1:numvert
        normal2 = PG.normal(:,j);
        tangent2 = PG.tangent(:,j);
        p2 = PG.vertex(:,j);
        s2 = PG.S(PG.VL(j));
        SegEnds = EqEEFinder(p1,s1,normal1,tangent1,p2,s2,normal2,tangent2);
        if ~isempty(SegEnds)
            EqCurves{end+1} = SegEnds;
        end
    end
end
%% find edge-vertex equilibrium curves:
%% tbd
end

function SegEnds = EqEEFinder(p1,s1,normal1,tangent1,p2,s2,normal2,tangent2)
% p_i, s_i, normal_i, tangent_i are the start of the ith edge (position and s-value), its normalized
% normal and non-normalized tangent (with length equal to the edge's)
% SegEnds = [s1start,s1end;s2start,s2end] describes the line of equilibrium
SegEnds = [];
% first, check if the edges can support an equilibrium at all:
if normal1(1)*normal2(1)>0 || all([normal1(2),normal2(2)]<0) % each of these (reversed) conditions is necessary
    return
end

%% find the intersections of the normals with the vertical line [0;t]
hln = [0;1];
if abs(det([normal1 hln]))<1E-10 || abs(det([normal2 hln]))<1E-10
    return
end
linsol11 = [hln -normal1]\p1;
t1(1:2,1) = linsol11;
linsol12 = [hln -normal1]\(p1+tangent1);
t1(1:2,2) = linsol12;
linsol21 = [hln -normal2]\p2;
t2(1:2,1) = linsol21;
linsol22 = [hln -normal2]\(p2+tangent2);
t2(1:2,2) = linsol22;

%% sort the intersection according to height, for each edge
s1start = s1;
s1end = s1+norm(tangent1);
signs1 = 1;
if t1(1,1)>t1(1,2)
    t1 = t1(:,[2,1]);
    s1start = s1+norm(tangent1);
    s1end = s1;
    signs1 = -1;
end
s2start = s2;
s2end = s2+norm(tangent2);
signs2 = 1;
if t2(1,1)>t2(1,2)
    t2 = t2(:,[2,1]);
    s2start = s2+norm(tangent2);
    s2end = s2;
    signs2 = -1;
end
s1st = s1start;
s2st = s2start;
s1en = s1end;
s2en = s2end;
% now each edge's intersctions are in ascending order
%% find the intersection of the edges' intersection ranges

% second, divide into cases:
% %if all the intersections are on the positive side of each edge:
% if all(t1(2,:)>0) && all(t2(2,:)>0) || all(t1(2,:)<0) && all(t2(2,:)<0)
    % find the intersection of ranges, and if necessary, decrease the s_i:
    if t1(1,1)>t2(1,1)
        s2start = s2start + signs2*(t1(1,1)-t2(1,1))/(t2(1,2)-t2(1,1))*norm(tangent2);
    else
        s1start = s1start + signs1*(t2(1,1)-t1(1,1))/(t1(1,2)-t1(1,1))*norm(tangent1);
    end
    if t1(1,2)<t2(1,2)
        s2end = s2end - signs2*(t2(1,2)-t1(1,2))/(t2(1,2)-t2(1,1))*norm(tangent2);
    else
        s1end = s1end - signs1*(t1(1,2)-t2(1,2))/(t1(1,2)-t1(1,1))*norm(tangent1);
    end
    if any([s1start,s1end]>max(s1st,s1en)) || any([s1start,s1end]<min(s1st,s1en)) || any([s2start,s2end]>max(s2st,s2en)) || any([s2start,s2end]<min(s2st,s2en))
        return
    else
        SegEnds = [s1start,s1end;s2start,s2end];
        return
    end  
% end
end
