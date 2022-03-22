syms dp1 dp2 vj1 vj2 vk1 vk2 tk1 tk2 theta Xcm Ycm a

dp = [dp1;dp2];
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
vj = [vj1;vj2];
vk = [vk1;vk2];
tk = [tk1;tk2];
p2 = dp;
d=[Xcm;Ycm];
eq1 = d+R*vj;
eq2 = d+R*(vk+a*tk)-p2;

sol = solve([eq1;eq2],[Xcm,Ycm,theta,a]);
%%
clc
% sol.Xcm
sol.Ycm
collect(sol.Ycm)
simplify(sol.Ycm)
% sol.theta
% collect(sol.theta)
simplify(sol.theta)
% sol.a
%% condition for equilibrium
syms nj11 nj12 nj21 nj22
nj1 = [nj11;nj12];
nj2 = [nj21;nj22];
J = [0 -1; 1 0];
% Xcmv = sol.Xcm(1);
av1 = sol.a(1);
thetav1 = sol.theta(1);
nk = J*tk;
Rv = subs(R,theta,thetav1);
supp2 = Rv*(vk+av1*tk); %relative to c.o.m 
alphabeta1 = [Rv*nk [0;-1]]\(-supp2);
% alphav1 = simplify(alphabeta1(1));
betav1 = simplify(alphabeta1(2))
supp1 = Rv*vj; %relative to c.o.m
supp1ToT = [0;betav1]-supp1;
angleFromSupp1 = simplify(atan2(supp1ToT(2),supp1ToT(1)))


av2 = sol.a(2);
thetav2 = sol.theta(2);
Rv = subs(R,theta,thetav2);
supp2 = Rv*(vk+av2*tk); %relative to c.o.m 
alphabeta2 = [Rv*nk [0;-1]]\(-supp2);
% alphav = simplify(alphabeta2(1));
betav2 = simplify(alphabeta2(2))
supp1 = Rv*vj; %relative to c.o.m
supp1ToT = [0;betav2]-supp1;
angleFromSupp2 = simplify(atan2(supp1ToT(2),supp1ToT(1)))
%%
syms anglefromsupp
njeff = [cos(anglefromsupp);sin(anglefromsupp)];
J = [0 -1;1 0];
nk = J*tk;
lamda = [njeff R*nk]\[0;1]