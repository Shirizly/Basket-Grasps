function [y1,y2,y3,y4] = EEheight1(dp,nj,nk,vj,vk)
dp1 = dp(1);
dp2 = dp(2);
nj1 = nj(1);
nj2 = nj(2);
nk1 = nk(1);
nk2 = nk(2);
vj1 = vj(1);
vj2 = vj(2);
vk1 = vk(1);
vk2 = vk(2);
Avector = EEheightAvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
a1 = Avector(1);
a2 = Avector(2);
a3 = Avector(3);
a4 = Avector(4);
a5 = Avector(5);
p1 = a1 - a3;
p2 = 2*a5 - 2*a2;
p3 = 4*a4 - 2*a1;
p4 = 2*a2 + 2*a5;
p5 = a1 + a3;
p = (8*p1*p3-3*p2^2)/(8*p1^2);

q = (p2^3-4*p1*p2*p3+8*p1^2*p4)/(8*p1^3);

del0 = p3^2-3*p2*p4+12*p1*p5;

del1 = 2*p3^3-9*p2*p3*p4+27*p2^2*p5+27*p1*p4^2-72*p1*p3*p5;

Q = ((del1+sqrt(del1^2-4*del0^3))/2)^(1/3);

S = 0.5*sqrt(-2*p/3+1/(3*p1)*(Q+del0/Q));


x1 = -p2/(4*p1)-S+0.5*sqrt(-4*S^2-2*p+q/S);
x2 = -p2/(4*p1)-S-0.5*sqrt(-4*S^2-2*p+q/S);
x3 = -p2/(4*p1)+S+0.5*sqrt(-4*S^2-2*p-q/S);
x4 = -p2/(4*p1)+S-0.5*sqrt(-4*S^2-2*p-q/S);
theta1 = NaN;
if isreal(x1)
theta1 = 2*atan2(x1,1);
end
theta2 = NaN;
if isreal(x2)
theta2 = 2*atan2(x2,1);
end
theta3 = NaN;
if isreal(x3)
theta3 = 2*atan2(x3,1);
end
theta4 = NaN;
if isreal(x4)
theta4 = 2*atan2(x4,1);
end


Dvector = EEheightDvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
Terms1 = [ cos(theta1)^2, cos(theta1)*sin(theta1), cos(theta1), sin(theta1)^2, sin(theta1), 1];
Terms2 = [ cos(theta2)^2, cos(theta2)*sin(theta2), cos(theta2), sin(theta2)^2, sin(theta2), 1];
Terms3 = [ cos(theta3)^2, cos(theta3)*sin(theta3), cos(theta3), sin(theta3)^2, sin(theta3), 1];
Terms4 = [ cos(theta4)^2, cos(theta4)*sin(theta4), cos(theta4), sin(theta4)^2, sin(theta4), 1];

y1 = Terms1*Dvector.';
y2 = Terms2*Dvector.';
y3 = Terms3*Dvector.';
y4 = Terms4*Dvector.';
end
