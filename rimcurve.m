function Rim = rimcurve(PG,basepos,pp)
% input: Polygon, positions of fingers, parametrization of
% puncture point on the rim
% output: data structure detailing the segments of the rim
Rim = cell(0);
vertex = PG.vertex;
nvert = size(vertex,2);
%% First, create list of c-space positions (x1,theta1,x2,theta2) at which each vertex crosses the rim height
vertpara = 100*ones(2,4,nvert); % 100 value marks empty column
Uf = zeros(1,2);
for finger = 1:2
    Uf(finger) = pp(3)-basepos(2,finger); % to work with height relative to contacting finger
    for v = 1:nvert
        vert = vertex(:,v);
        rad = norm(vert);
        if rad > Uf(finger) % vertex crosses the height at two angles
            thetavert = 3*pi()/2-atan2(vert(2),vert(1))); % angle at which com is directly above vertex
            alpha = pi()/2-asin(Uf(finger)/rad); % angle needed to rotate to bring com to rim height, relative to directly above
            theta1 = wrapToPi(thetavert+alpha);
            theta2 = wrapToPi(thetavert-alpha);
            x1 = basepos(1,finger)-sin(alpha)*rad;
            x2 = basepos(1,finger)+sin(alpha)*rad;
            vertpara(finger,:,v) = [x1,theta1,x2,theta2].';
        end
    end
end
%% Second, create segments of pp-height curves in object c-space
segment = cell(0);
vertpara(:,:,end+1) = vertpara(:,:,1);
for finger = 1:2
for i=1:nvert
    

end
end
end

