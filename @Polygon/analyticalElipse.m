function curve = analyticalElipse(PG,rect,sigma)

curve = sigma;
edge1 = PG.get('edgeNum',mean(rect(1,1:2)));
edge2 = PG.get('edgeNum',mean(rect(1,3:4)));
t1 = PG.tangent(:,edge1);
t2 = PG.tangent(:,edge2);
v11 = PG.X(PG.VL(edge1)).';
v21 = PG.X(PG.VL(edge2)).';
A = [t1,-t2];
b = v11-v21;
center_rel = -(A.'*A)\A.'*b;
center = center_rel + [rect(1,1);rect(1,3)];
s1 = zeros(5,1);
s2 = zeros(5,1);
sig = zeros(5,1);
ord1 = [1,2,2,1,1];
ord2 = [1,1,2,2,1];
for i = 1:5
    s1(i) = rect(ord1(i));
    s2(i) = rect(ord2(i));
    [f1,f2] = PG.get('2pos',s1(i),s2(i));
    sig(i) = norm(f1-f2);
end
pair = zeros(2,2);
counter = 1;
for i=1:4
    alpha = (sig(i)-sigma)/(sig(i)-sig(i+1));
    if alpha>=0 && alpha<1
        fixedp = (mod(i,2)==0)*s1(i)+(mod(i,2)==1)*s2(i);
        changep = (mod(i,2)==1)*s1(i)+(mod(i,2)==0)*s2(i);
        newp = correct_position(PG,fixedp,changep,sigma);
        pair(counter) = [fixedp,newp];
        counter = counter+1;
    end
end
curve = {pair,center};
end