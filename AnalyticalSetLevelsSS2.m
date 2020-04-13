

com = [1;3];
t = [-1;1];
t = t/norm(t);
send = 3;
figure
hold on
cstart = -5;
cend = 5;
cvector = [cstart:0.2:cstart+2,cstart+2.1:0.5:cend-2.1,cend-2:0.2:cend];
for c=cvector
    theta = [-pi():0.0001:pi()].';
    vtheta = [sin(theta),-cos(theta)];
    s = (vtheta*com-c)./(vtheta*t);
    p = (c-cstart)/(cend-cstart);
    col = [max(0,1-abs(3*(p-1))),max(0,1-abs(3*(p-2/3))),max(0,1-abs(3*(p-1/3)))];
    plot(s,theta,'.','Color',col,'lineWidth',1,'markerSize',1)
end
J = [0,1;-1,0];
sh = [t,J*t]\com;
h = sh(2);
theta0 = atan2(t(2),t(1));
theta = [-pi():0.0001:pi()].';
xtheta = [cos(theta),sin(theta)];
sx = (xtheta*com)./(xtheta*t);
plot(sx,theta,'.k','markerSize',15)
for c=[h,-h]
    theta = [-pi():0.0001:pi()].';
    vtheta = [sin(theta),-cos(theta)];
    s = (vtheta*com-c)./(vtheta*t);
    p = (c-cstart)/(cend-cstart);
    col = [max(0,1-abs(3*(p-1))),max(0,1-abs(3*(p-2/3))),max(0,1-abs(3*(p-1/3)))];
    plot(s,theta,'.','Color',col)
    plot(s,theta,'.','Color',col,'markerSize',6)
    thetatemp = (theta0-(1-sign(c))*pi()/2);
    plot([0,send],ones(2,1)*thetatemp,'Color',col,'lineWidth',3)
    plot(sh(1),thetatemp,'+b','markerSize',15,'lineWidth',8)
    plot(0,thetatemp,'+k','markerSize',15,'lineWidth',8)
    plot(send,thetatemp,'+k','markerSize',15,'lineWidth',8)
end

plot([0,1],[-0.67,-0.67])
smin = find(sx == min(abs(sx)));
plot(sx(smin),theta(smin),'+r','markerSize',15,'lineWidth',8)
plot(sx(smin),theta(smin)-pi(),'+g','markerSize',15,'lineWidth',8)
smax = find(sx == send - min(abs(sx-send)));
plot(sx(smax),theta(smax),'+r','markerSize',15,'lineWidth',8)
plot(sx(smax),theta(smax)-pi(),'+g','markerSize',15,'lineWidth',8)
hold off
xlim([0,send]);
ylim([-pi(),pi()]);
title('Generic strip','fontSize',24)
xlabel(['s'],'fontSize',22)
ylabel('\theta','fontSize',22)