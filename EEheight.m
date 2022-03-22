function [sol1,sol2,sol3,sol4] = EEheight(dp,nj,nk,vj,vk,dels1,dels2)
% Function returns the obect configuration at Edge-Edge maximum height equilibrium
% stances. 
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
% transform the quadratic equation in cosine and sine to 4th order poly.
Avector = EEheightAvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
a1 = Avector(1);
a2 = Avector(2);
a3 = Avector(3);
a4 = Avector(4);
a5 = Avector(5);
% find the coefficients for the canonical form of 4th order poly.
p1 = a1 - a3;
p2 = 2*a5 - 2*a2;
p3 = 4*a4 - 2*a1;
p4 = 2*a2 + 2*a5;
p5 = a1 + a3;
% solve the 4th order poly.
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
% check each solution for realness, and if real, calculate theta
eps = 1E-7;
theta1 = NaN;
if abs(imag(x1))<eps
theta1 = 2*atan2(real(x1),1);
end
theta2 = NaN;
if abs(imag(x2))<eps
theta2 = 2*atan2(real(x2),1);
end
theta3 = NaN;
if abs(imag(x3))<eps
theta3 = 2*atan2(real(x3),1);
end
theta4 = NaN;
if abs(imag(x4))<eps
theta4 = 2*atan2(real(x4),1);
end



% calculate the second derivative of the height according to theta
Cvector = EEheightCvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
Terms1 = [ cos(theta1)^2, cos(theta1)*sin(theta1), cos(theta1), sin(theta1)^2, sin(theta1)];
Terms2 = [ cos(theta2)^2, cos(theta2)*sin(theta2), cos(theta2), sin(theta2)^2, sin(theta2)];
Terms3 = [ cos(theta3)^2, cos(theta3)*sin(theta3), cos(theta3), sin(theta3)^2, sin(theta3)];
Terms4 = [ cos(theta4)^2, cos(theta4)*sin(theta4), cos(theta4), sin(theta4)^2, sin(theta4)];

diff2d1 = Terms1*Cvector.';
diff2d2 = Terms2*Cvector.';
diff2d3 = Terms3*Cvector.';
diff2d4 = Terms4*Cvector.';

diffd1 = Terms1*Avector.';
diffd2 = Terms2*Avector.';
diffd3 = Terms3*Avector.';
diffd4 = Terms4*Avector.';
%% find if the contact is within the actual edge and check feasibility
% also find the COM position relative to finger 1
%old way that is problematic
% [S1vector,S2vector] = EEheightSvectors(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
% Terms1 = [cos(theta1), sin(theta1), 1];
% Terms2 = [cos(theta2), sin(theta2), 1];
% Terms3 = [cos(theta3), sin(theta3), 1];
% Terms4 = [cos(theta4), sin(theta4), 1];
% 
% S11 = Terms1*S1vector.';
% S12 = Terms2*S1vector.';
% S13 = Terms3*S1vector.';
% S14 = Terms4*S1vector.';
% S21 = Terms1*S2vector.';
% S22 = Terms2*S2vector.';
% S23 = Terms3*S2vector.';
% S24 = Terms4*S2vector.';

e = [0,1];
J = [0 1; -1 0];
A = [-J*nj J*nk];
B = [J*nj J*nk];
dv = vj-vk;
sv = vj+vk;
R = [cos(theta1) -sin(theta1);sin(theta1) cos(theta1)];
feasicond1 = ((e*R*nj>0)||(e*R*nk>0))&&((e*J*R*nj)*(e*J*R*nk)<0);
if ~any(isnan(R))
lamda = [R*nj R*nk]\[0;1];
feasicond1 = feasicond1&& all(lamda>0);
end
s1 = A\(R.'*dp+dv);
d1 = 0.5*(dp-R*(sv+B*s1));
R = [cos(theta2) -sin(theta2);sin(theta2) cos(theta2)];
feasicond2 = ((e*R*nj>0)||(e*R*nk>0))&&((e*J*R*nj)*(e*J*R*nk)<0);
if ~any(isnan(R))
lamda = [R*nj R*nk]\[0;1];
feasicond2 = feasicond2&& all(lamda>0);
end
s2 = A\(R.'*dp+dv);
d2 = 0.5*(dp-R*(sv+B*s2));
R = [cos(theta3) -sin(theta3);sin(theta3) cos(theta3)];
feasicond3 = ((e*R*nj>0)||(e*R*nk>0))&&((e*J*R*nj)*(e*J*R*nk)<0);
if ~any(isnan(R))
lamda = [R*nj R*nk]\[0;1];
feasicond3 = feasicond3&& all(lamda>0);
end
s3 = A\(R.'*dp+dv);
d3 = 0.5*(dp-R*(sv+B*s3));
R = [cos(theta4) -sin(theta4);sin(theta4) cos(theta4)];
feasicond4 = ((e*R*nj>0)||(e*R*nk>0))&&((e*J*R*nj)*(e*J*R*nk)<0);
if ~any(isnan(R))
lamda = [R*nj R*nk]\[0;1];
feasicond4 = feasicond4 && all(lamda>0);
end
s4 = A\(R.'*dp+dv);
d4 = 0.5*(dp-R*(sv+B*s4));
% remove solutions that don't fit the two criterions:
if  any(s1<0) || s1(1)>dels1 || s1(2)>dels2 || diff2d1>0 || ~feasicond1
    theta1 = NaN;
end
if  any(s2<0) || s2(1)>dels1 || s2(2)>dels2 || diff2d2>0 || ~feasicond2
    theta2 = NaN;
end
if  any(s3<0) || s3(1)>dels1 || s3(2)>dels2 || diff2d3>0 || ~feasicond3
    theta3 = NaN;
end
if  any(s4<0) || s4(1)>dels1 || s4(2)>dels2 || diff2d4>0 || ~feasicond4
    theta4 = NaN;
end
% find the translation vector for each solution - switched to more accurate
% version
% Dmatrix = EEheightDmatrix(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
% Terms1 = [ cos(theta1)^2, cos(theta1)*sin(theta1), cos(theta1), sin(theta1)^2, sin(theta1), 1];
% Terms2 = [ cos(theta2)^2, cos(theta2)*sin(theta2), cos(theta2), sin(theta2)^2, sin(theta2), 1];
% Terms3 = [ cos(theta3)^2, cos(theta3)*sin(theta3), cos(theta3), sin(theta3)^2, sin(theta3), 1];
% Terms4 = [ cos(theta4)^2, cos(theta4)*sin(theta4), cos(theta4), sin(theta4)^2, sin(theta4), 1];

% d1 = diag(repmat(Terms1,2,1)*Dmatrix.');
% d2 = diag(repmat(Terms2,2,1)*Dmatrix.');
% d3 = diag(repmat(Terms3,2,1)*Dmatrix.');
% d4 = diag(repmat(Terms4,2,1)*Dmatrix.');
% if any(~isnan([theta1,theta2,theta3,theta4]))
%     disp('check')
% end

% prepare the solution vectors for output
sol1 = [d1;theta1];
sol2 = [d2;theta2];
sol3 = [d3;theta3];
sol4 = [d4;theta4];
end
