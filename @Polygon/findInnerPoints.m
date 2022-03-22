function Distance = findInnerPoints(PG,otherFingerRelative,thetav)
Distance = zeros(length(PG.S),length(thetav));
for i = 1:length(PG.S)
    s = PG.S(i);
    for j = 1:length(thetav)
        theta = thetav(j);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        contact = PG.get('1Pos',s);
        vertex = R*(PG.vertex-contact);
        mindist = p_poly_dist1(otherFingerRelative(1), otherFingerRelative(2), vertex(1,:), vertex(2,:));
        Distance(i,j) = mindist;
%         if mindist>0.1
%             disp('weird')
%         end
%         check = inpolygon(otherFingerRelative(1),otherFingerRelative(2),vertex(1,:),vertex(2,:));
%         mindist = 100;
%         for k = 1:size(vertex,2)
%             OFRTV = otherFingerRelative-vertex(:,k);
%             tangent = R*PG.tangent(:,k); %other finger relative to vertex
%             normal = R*PG.normal(:,k);
%             p = [tangent normal]\OFRTV;
%             dist = abs(p(2)) + 1000*~(p(1)>=0&&p(1)<=1); %distance is what got calculated if the intersection of the normal and the edge is inside the edge, otherwise arbitrarily high
%             mindist = min(mindist,dist);
%         end
        
%         Distance(i,j) = (-2*check+1)*mindist;

    end
end
end


% function:	p_poly_dist
% Description:	distance from pont to polygon whose vertices are specified by the
%              vectors xv and yv
% Input:  
%    x - point's x coordinate
%    y - point's y coordinate
%    xv - vector of polygon vertices x coordinates
%    yv - vector of polygon vertices x coordinates
% Output: 
%    d - distance from point to polygon (defined as a minimal distance from 
%        point to any of polygon's ribs, positive if the point is outside the
%        polygon and negative otherwise)
% Routines: p_poly_dist.m
% Revision history:
%    7/9/2006  - case when all projections are outside of polygon ribs
%    23/5/2004 - created by Michael Yoshpe 
% Remarks:
function d = p_poly_dist1(x, y, xv, yv) 
% If (xv,yv) is not closed, close it.
xv = xv(:);
yv = yv(:);
Nv = length(xv);
if ((xv(1) ~= xv(Nv)) || (yv(1) ~= yv(Nv)))
    xv = [xv ; xv(1)];
    yv = [yv ; yv(1)];
    Nv = Nv + 1;
end
% linear parameters of segments that connect the vertices
A = -diff(yv);
B =  diff(xv);
C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);
% find the projection of point (x,y) on each rib
AB = 1./(A.^2 + B.^2);
vv = (A*x+B*y+C);
xp = x - (A.*AB).*vv;
yp = y - (B.*AB).*vv;
% find all cases where projected point is inside the segment
idx_x = (((xp>=xv(1:end-1)) & (xp<=xv(2:end))) | ((xp>=xv(2:end)) & (xp<=xv(1:end-1))));
idx_y = (((yp>=yv(1:end-1)) & (yp<=yv(2:end))) | ((yp>=yv(2:end)) & (yp<=yv(1:end-1))));
idx = idx_x & idx_y;
% distance from point (x,y) to the vertices
dv = sqrt((xv(1:end-1)-x).^2 + (yv(1:end-1)-y).^2);
if(~any(idx)) % all projections are outside of polygon ribs
   d = min(dv);
else
   % distance from point (x,y) to the projection on ribs
   dp = sqrt((xp(idx)-x).^2 + (yp(idx)-y).^2);
   d = min(min(dv), min(dp));
end
if(inpolygon(x, y, xv, yv)) 
   d = -d;
end
end