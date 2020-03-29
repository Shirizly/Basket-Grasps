function [ss_max,ss_min,ss_saddle] = SSNodes(PG)
[sss,ssh,sstheta]= PG.SSEqcheck();
[hooks,~,~] = PG.HookEqcheck();
[mins,minh,mintheta] = PG.SSMincheck();
exclude = [];
for i =1:length(hooks)
exclude = find(mins==hooks(i));
end
type = zeros(size(mins));
type(exclude) = ones(size(exclude));
ss_min = [mins.',minh.',mintheta.',type.'];
ss_saddle = [sss.',ssh.',sstheta.'];
ss_max = [mins.',-minh.',wrapToPi(mintheta.'+pi())];
ss_virtual = zeros((length(PG.VL)-1)*4,3);
for i = 1:length(PG.VL)-1
    normal = PG.normal(:,i);
    thetaSaddle = wrapToPi(pi()/2-atan2(normal(2),normal(1))); %angle of the higher saddle line
    ss_virtual(4*i-3:4*i-2,[1,3]) = [PG.S(PG.VL([i,i+1])),thetaSaddle*ones(2,1)];
    ss_virtual(4*i-1:4*i,[1,3]) = [PG.S(PG.VL([i,i+1])),wrapToPi(thetaSaddle+pi())*ones(2,1)];
%     if abs(thetaSaddle - pi())<1E-8 || abs(thetaSaddle + pi())<1E-8
%     % dealing with virtual nodes exactly on the edges of the plane
%         ss_virtual(4*i-3:4*i-2,[1,3]) = [PG.S(PG.VL([i,i+1])),-thetaSaddle*ones(2,1)];
%     end
%     if abs(thetaSaddle)<1E-8
%         ss_virtual(4*i-3:4*i-2,[1,3]) = [PG.S(PG.VL([i,i+1])),-thetaSaddle*ones(2,1)];
%     end
        
end
for i = 1:size(ss_virtual,1)
    dist = PG.get('comLocationSS',[0;0],ss_virtual(i,1),ss_virtual(i,3));
    ss_virtual(i,2) = dist(2);
end
ss_saddle = [ss_saddle;ss_virtual];
end