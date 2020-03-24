function new_point = correct_position(PG,s1,s2,sig)
% this function corrects errors in s2 given a point s1,s2 near the
% sigma-level-set of a polygon, calculated by the contour function
%(the contour function generates a noisy curve, close to the real one, but
%with enough noise to generate significant errors
[fs1,fs2] = PG.get('2Pos',s1,s2);
sig_cur = norm(fs2-fs1);
t2 = PG.tangent(:,PG.get('edgeNum',s2));
t2 = t2/norm(t2);
B = 2*(fs2-fs1).'*t2;
C = sig_cur^2-sig^2;
delta_s = (-B+sign(B)*sqrt(B^2-4*C))/2;
if isreal(delta_s)
    new_point = s2+delta_s;
else
    new_point = s2;
end
    