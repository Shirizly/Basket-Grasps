angle1 = 1;
angle2 = 2;
vert1 = PG.vertex(:,1);
vert2 = PG.vertex(:,2);
rad1 = norm(vert1);
rad2 = norm(vert2);
theta = angle1:0.05:angle2;
alpha1 = atan2(vert1(2),vert1(1));
x1 = cos(alpha1-theta)*rad1;
y1 = sin(alpha1-theta)*rad1;
alpha2 = atan2(vert2(2),vert2(1));
x2 = cos(alpha2-theta)*rad2;
y2 = sin(alpha2-theta)*rad2;
relative = [x2-x1;y2-y1];
dist = relative(1,:).*relative(1,:)+relative(2,:).*relative(2,:);
plot(x1,y1,x2,y2,(x1+x2)/2,(y1+y2)/2)
hold on
for i=1:length(x1)
    plot([x1(i),x2(i)],[y1(i),y2(i)],'b')
end
figure
plot(theta,atan2(relative(2,:),relative(1,:)))
