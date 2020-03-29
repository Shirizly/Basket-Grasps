function Distance = findInnerPoints(PG,otherFingerRelative,thetav)
Distance = zeros(length(PG.S),length(thetav));
for i = 1:length(PG.S)
    s = PG.S(i);
    for j = 1:length(thetav)
        theta = thetav(j);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        contact = PG.get('1Pos',s);
        vertex = R*(PG.vertex-contact);
        check = inpolygon(otherFingerRelative(1),otherFingerRelative(2),vertex(1,:),vertex(2,:));
        mindist = 100;
        for k = 1:size(vertex,2)
            OFRTV = otherFingerRelative-vertex(:,k);
            tangent = R*PG.tangent(:,k); %other finger relative to vertex
            normal = R*PG.normal(:,k);
            p = [tangent normal]\OFRTV;
            dist = abs(p(2)) + 1000*~(p(1)>=0&&p(1)<=1); %distance is what got calculated if the intersection of the normal and the edge is inside the edge, otherwise arbitrarily high
            mindist = min(mindist,dist);
        end
        Distance(i,j) = (-2*check+1)*mindist;
    end
end
end