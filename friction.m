%% Calculating friction according to rotation experiment
theta_crit = 40/180*pi();
%% calculate finger positions based on finger opening, measure com
% sigma = 0.112; %units are meters 
% syms sc
% object = sc*[-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].';
% com = [0;0];
% object = object-com;
% %% measure the contact points
% c1 = 0.0145;
% c2 = 0.0495;
% dir1 = (object(:,1)-object(:,end));
% dir1 = dir1/sc;
% dir1 = dir1/norm(dir1);
% dir2 = (object(:,3)-object(:,2));
% dir2 = dir2/sc;
% dir2 = dir2/norm(dir2);
% c1m = object(:,1)+c1*dir1;
% c2m = object(:,2)+c2*dir2;
% eq = norm(c2m-c1m)-sigma;
% scale = double(solve(eq));
%%
scale = 0.0205;
object = double(scale*[-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].');
com = [0;0.005];
object = object-com;
com = [0;0];
PG = Polygon(object,com); 
PG.drawPolygon()
%% define finger normals and tangents
N1 = PG.normal(:,7);
N2 = PG.normal(:,2);
J = [0 1;-1 0];
t1 = J*N1;
t2 = J*N2;

c1m = PG.vertex(:,1)-c1*t1;
c2m = PG.vertex(:,2)+c2*t2;
norm(c1m-c2m)
%% find normals intersection point
sol = [N1 -N2]\(c2m-c1m);
inter = c1m + sol(1)*N1;

%% construct linear equations in the reaction forces
syms f1 f2 n1 n2 
% theta_crit = 0;
mg = [-sin(theta_crit);-cos(theta_crit)];
eq12 = f1*t1 + n1*N1 + f2*t2 + n2*N2 + mg;
mge = [mg;0];
mgtorque = cross(mge,[inter;0]);
eq3 = f1*sol(1) + f2*sol(2) + mgtorque(3);
F1 = n1*N1 + f1*t1;
F2 = n2*N2 + f2*t2;
F1e = [F1;0];
F2e = [F2;0];
eq4 = cross(F1e,[-c1m;0])+cross(F2e,[-c2m;0]);
fric1min = 10;
fric2min = 10;
minSumReact = 10;
x = [f1;f2;n2];
for i = 0.01:0.001:1.5
    EQ = subs([eq12;eq3],n1,i);
    % turn system of linear equations to matrix form
    
    A = double(jacobian(EQ,x));
    b = double(A*x-EQ);%[-mg;-mgtorque(3)]
    
    % solve matrix equation for the reaction forces
    xsol = A\b;
%     SumReact = sqrt(i^2+xsol(1).*xsol);
    
    
    if xsol(3)>0 
        fric1 = abs(xsol(1)/i);
        fric2 = abs(xsol(2)/xsol(3));
        if max(fric1,fric2)<max(fric1min,fric2min)
            fric1min = fric1;
            fric2min = fric2;
            forces = [xsol(1:2);i;xsol(3)];
        end
    end
end
max(fric1min,fric2min)
f1 = forces(1);
f2 = forces(2);
n1 = forces(3);
n2 = forces(4);
%% draw object with reaction forces
PG.drawPolygon();
hold on
scale = 0.3;
temp = [PG.com,PG.com+scale*mg];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--r')
temp = [PG.com,PG.com-scale*mg];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--r')
temp = [c1m,c2m];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'ok')
temp = [c1m,c1m+scale*(n1*N1+f1*t1)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--b')
temp = [c2m,c2m+scale*(n2*N2+fric1min/fric2min*f2*t2)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--b')
temp = [c1m,c1m+scale*(n1*N1)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'b')
temp = [c2m,c2m+scale*(n2*N2)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'b')
temp = [c1m,c1m+scale*(n1*N1-f1*t1)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--b')
temp = [c2m,c2m+scale*(n2*N2-fric1min/fric2min*f2*t2)];
xtemp = temp(1,:);
ytemp = temp(2,:);
plot(xtemp,ytemp,'--b')
camroll(theta_crit*180/pi())
ax = gca;
set(ax,'Visible','off')
