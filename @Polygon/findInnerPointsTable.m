function Distance = findInnerPointsTable(PG,ContactHeight,thetav)
Distance = zeros(length(PG.S),length(thetav));
vert = PG.vertex;
for i = 1:length(PG.S)
    s = PG.S(i);
    for j = 1:length(thetav)
        theta = thetav(j);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        contact = PG.get('1Pos',s);
        vertex = R*(vert-contact); %relative to contact
        
        convertex = vertex(:,PG.CVL(1:end-1)); %relative to contact
        mindisttable = min(convertex(2,:))+ContactHeight;
        
        Distance(i,j) = mindisttable;
    end
end
end