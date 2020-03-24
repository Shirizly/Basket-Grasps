function [ss_max,ss_min,ss_saddle] = SSNodes(PG)
[sss,ssh,sstheta]= PG.SSEqcheck();
[hooks,hookh,hooktheta] = PG.HookEqcheck();
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
end