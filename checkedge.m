function ssind = checkedge(contv,VL,ssind,sss)
lim1 = find(VL<=sss,1,'last');
lim2 = find(VL>sss,1,'first');
remove = [];
for i=1:length(ssind)
    check = and(contv(ssind(i))>=VL(lim1),contv(ssind(i))<VL(lim2));
    if ~check
        remove = [remove i];
    end
end
ssind(remove) = [];
end