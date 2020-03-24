

A = 1;
B = 1;
send = 3;
theta0 = 0;
figure
hold on
for c=-3:0.1:2
   theta = -pi():0.001:pi();
   thetag = theta-theta0;
   s = A*tan(thetag)+B*ones(size(thetag))-c./cos(thetag);
   p = (c+3)/5;
   col = [max(0,1-abs(3*(p-1))),max(0,1-abs(3*(p-2/3))),max(0,1-abs(3*(p-1/3)))];
   plot(s,theta,'Color',col)

end

for c=[0.999999,-0.999999]
   theta = -pi():0.0005:pi();
   thetag = theta-theta0;
   s = A*tan(thetag)+B*ones(size(thetag))-c./cos(thetag);
   p = (c+3)/5;
   col = [max(0,1-abs(3*(p-1))),max(0,1-abs(3*(p-2/3))),max(0,1-abs(3*(p-1/3)))];
   plot(s,theta,'.','Color',col,'markerSize',12)
   plot([0,send],ones(2,1)*sign(c)*pi()/2,'Color',col,'lineWidth',4)
end



hold off
xlim([0,send]);
ylim([-pi(),pi()]);