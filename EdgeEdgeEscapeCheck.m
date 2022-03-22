% clearvars

%% Define the polygon
% %% Define the polygon
clc
object = 0;
switch object
    case 0 % pentagon
        theta = 0;%-pi()/2;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
        PG = Polygon(R*[3,-4;3,4;-1.5,5;-3,0;-1.5,-5].'+[0;3],[10;0]);
        subFolderName = 'Examples for research proposal\Pentagon'; % where to save results
        res = 50;
        f2 = [8;1];
        f1 = [0;0];
        dp = f2-f1;
    case 1 % cup
       PG = Polygon([-3,2;-2,-2;2,-2;2,1;4,1;4,2].'+[0;0],[0;0]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results
       f2 = [1;1];
       f1 = [0;0];
       dp = f2-f1;
%        basepos = [f1,f2];
       res = 50;
    case 2 % bottle
        load object5.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [0.06;0.01];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 1500;
        starting_node = [2,21];
    case 3 % Elon
        load object.mat
        %     object = [object(:,7:8) object(:,1:6)];
        com = [0.045,0.07].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Elon'; % where to save results
    case 4 % Crown
        object = [-2,-1;2,-1;3,3;1,0;0,1;-1,0;-4,1].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 50;
         starting_node = [2,4];    
    case 5 % new bottle
        load NewBottle.mat
        %     object = [object(:,7:8) object(:,1:6)];
%         com = [0.045,0.07].';
%         object = object + 0.001*randn(size(object));
%         theta =  1.515;
%         R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
%         object = R*object;
        com(2) = max(object(2,:))-com(2);
        object(2,:) = max(object(2,:))-object(2,:);
        object(2,1) = 0;
        object = round(object./20,1);
        com = round(com./20,0);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Bottle'; % where to save results
        f2 = [4;2];
        f1 = [0;0];
        basepos = [f1,f2];
        res = 50;
        starting_node = [2,24];
        case 6 % Gun
%         load gun2.mat
%         object = Round_object;
%         object = 0.8*[1,10;26,10;21,1;29,9;31,2;28,0;32,0;32,12;35,16;32,18;33,16;32,14;31,14;31,16;1,16;0.5,13.5;0,13;0.5,12.5].';
        
        object = 0.8*[1,10;26,10;21,1;29,9;32,0;32,12;36,12;40,8;41,9;35,14;31,16;1,16].';%32,14;31,14
        com = 0.8*[27;10];
        theta =  0.423;
        R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
         object = R*object;
         com = R*com;
%          object = round(object./20,1);
%          com = round(com./20,0);
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Gun'; % where to save results
        f2 = [18.2;16.09];
        f1 = [11.4;13.9];
        basepos = [f1,f2];
        res = 30;
        starting_node = [2,6];
        
        case 7 % Crown with more edges
        object = [-2,-1;2,-1;2.25,1;3.25,3;1,0;0,1;-1,0;-4,1;-2.75,0].';
        com = [0,0].';
        PG = Polygon(object,com);
        subFolderName = 'Examples for research proposal\Crown'; % where to save results
        f2rel = [2.7519;2.0074];
        f1rel = [-2.8470;-0.1530];
        f1 = [0;0];
        f2 = f2rel-f1rel;
        res = 50;
         starting_node = [2,4];  
         
             case 8 % modified cup
       
       PG = Polygon([-3,2;-2,-2;2,-2;2,-0.5;4,1.5;4,2].'+[0;3],[0;3]);
       subFolderName = 'Examples for research proposal\Cup'; % where to save results

%        basepos = [f1,f2];
       res = 50;
       
    case 9
        load DrawnObject.mat
        object = Round_object;
        PG = Polygon(object,[0;0]);
        res = 50;
end
[PG,S,X,VL] = PG.findBdyVariable(res);
%%
edge1 = 5;
edge2 = 2;
dels1 = PG.S(PG.VL(edge1+1))-PG.S(PG.VL(edge1));
dels2 = PG.S(PG.VL(edge2+1))-PG.S(PG.VL(edge2));

nj = PG.normal(:,edge1);
nj = nj/norm(nj);
nk = PG.normal(:,edge2);
nk = nk/norm(nk);
epsi = 1E-7;
J = [0 -1;1 0];
% if abs(nj.'*J*nk)<epsi
%     continue
% end
vj = PG.vertex(:,edge1);
vk = PG.vertex(:,edge2);
%%
syms theta
% theta = sol1(3,1)
% d = sol1(1:2,1)
dp = dp(:,1);%[5;1];
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
J = [0 1; -1 0];
A = [-J*nj J*nk];
B = [J*nj J*nk];
% dP = dp(:,1);
dv = vj-vk;
sv = vj+vk;


dels = A\(R.'*dp+dv);
d = 0.5*(dp-R*(sv+B*dels));
diffd = diff(d(2),theta);
diff2d = diff(diffd,theta);

%%
thetav = -pi():0.1:pi;
figure
hold on
plot(0,0,'or')
plot(dp(1),dp(2),'or')
for i=1:length(thetav)
    dv(:,i) = subs(d,theta,thetav(i));
    sf = subs(dels,theta,thetav(i));
    PG.drawPolygonMoved([0;0],thetav(i),dv(:,i),'k')
    pause(0.1)
end
%%
% figure
% plot(thetav,dv(2,:))
[~,thetamaxind] = findpeaks(dv(2,:));
thetamax = thetav(thetamaxind)
%%
figure
hold on
plot(0,0,'or')
plot(dp(1),dp(2),'or')
for i=1:length(thetamax)
thetam = double(thetamax(i))
dmax = double(subs(d,theta,thetamax(i)))
% sf = subs(dels,theta,thetamax(i))
PG.drawPolygonMoved([0;0],thetamax(i),dmax,'k')
end
%%
[sol1,sol2,sol3,sol4] = EEheight(dp,nj,nk,vj,vk,dels1,dels2)
%%
fun = matlabFunction(diffd);
theta0 = 0;
% feval(fun,0)
fzero(fun,theta0)
%%
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
[Avector2,Aterms] = coeffs(diffd,[cos(theta),sin(theta)])
Avector = EEheightAvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2)
a1 = Avector(1);
a2 = Avector(2);
a3 = Avector(3);
a4 = Avector(4);
a5 = Avector(5);
%% find the coefficients for the canonical form of 4th order poly.
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