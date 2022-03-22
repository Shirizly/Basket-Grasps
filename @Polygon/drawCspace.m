function drawCspace(PG,s1,s2,theta)
f1 = [0;0];
f1rel = PG.get('1Pos',s1);
f2rel = PG.get('1Pos',s2);
axle = f1;
d = f1-f1rel;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
f2 = f1 + R*(f2rel-f1rel);
fingers = [f1,f2];

sig = norm(f2-f1);
[PG,S,X,VL] = PG.findBdyVariable(PG.res);
Sigma=inter_finger_distance(X,X);
[cont_original] = PG.GetSigmaContours(Sigma,sig);
cont = PG.CleanContour(cont_original,fingers);

figure
% hold on

s = [s1,s2];
col = [1,0,0;0,0,1];
xc = [];
yc = [];
thetac = [];
for finger = 1:2
    for thetav = theta-pi():0.1:theta+pi()
        comrel = PG.get('comLocationSS',fingers(:,finger),0,thetav); 
        R = [cos(thetav) -sin(thetav); sin(thetav) cos(thetav)];
        pos = R*PG.vertex;
        xs = pos(1,:) + comrel(1);
        ys = pos(2,:) + comrel(2);
%         plot3(xs,thetav*ones(size(xs)),ys,'color',col(finger,:));
%         plot3(xs([1,end]),thetav*ones(1,2),ys([1,end]),'color',col(finger,:));
        xc(end+1,:) = xs;
        yc(end+1,:) = ys;
        thetac(end+1,:) = thetav*ones(size(xs));
    end
    hold on
    surf(xc,yc,thetac);
    surf(xc(:,[1,end]),yc(:,[1,end]),thetac(:,[1,end]));
    alpha 0.5
    xc = [];
    yc = [];
    thetac = [];
end
    
    con = cont{3};
    xds = zeros(size(con,2),1);
    yds = zeros(size(con,2),1);
    tds = zeros(size(con,2),1);
    for j=1:size(con,2)
        c = con(:,j);
        comrel = PG.get('comLocationSS',fingers(:,1),c(1),c(4));
        xds(j) = comrel(1);
        yds(j) = comrel(2);
        tds(j) = c(4);
    end
    plot3(xds,yds,tds,'.k','linewidth',3)
    zlim([-pi/4 pi/4])    
end

% for finger = 1:2
%     for s = PG.S.'
%         thetav = -pi()/3:0.1:pi()/3;
%         comrel = PG.get('comLocationSS',fingers(:,finger),s,0);
%         x = fingers(1,finger) + cos(thetav)*comrel(1) -sin(thetav)*comrel(2);
%         y = fingers(2,finger) + sin(thetav)*comrel(1) + cos(thetav)*comrel(2);
%         plot3(x,thetav,y,'color',col(finger,:));
%     end
% end