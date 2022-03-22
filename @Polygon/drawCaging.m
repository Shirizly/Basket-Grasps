function [finger1,finger2] = drawCaging(PG,BGS,dep,col)
for finger = 1:2
    sOri = BGS(finger);
    pOri(:,finger) = PG.get('1Pos',sOri);
    sCurve{finger} = findCurveS(PG,sOri,pOri(2,finger)-dep);
end
regionXc = cell(2,1); regionYc = cell(2,1);
for finger = 1:2
    regionX = zeros(1,length(sCurve{finger})); regionY = zeros(1,length(sCurve{finger}));
    sc = sCurve{finger};
    for j = 1:length(sc)
        pos = PG.get('1Pos',sc(j));
        regionX(j) =  pos(1);
        regionY(j) = pos(2);
    end
    regionXc{finger} = regionX;
    regionYc{finger} = regionY;
end

finger1X = [regionXc{1},regionXc{2}(end:-1:1)+pOri(1,1)-pOri(1,2)];
finger1Y = [regionYc{1},regionYc{2}(end:-1:1)+pOri(2,1)-pOri(2,2)];
finger1 = [finger1X;finger1Y];
plot(finger1X,finger1Y,'Color',col,'lineWidth',1)
finger2X = [regionXc{2},regionXc{1}(end:-1:1)+pOri(1,2)-pOri(1,1)];
finger2Y = [regionYc{2},regionYc{1}(end:-1:1)+pOri(2,2)-pOri(2,1)];
plot(finger2X,finger2Y,'Color',col,'lineWidth',1)
finger2 = [finger2X;finger2Y];
end

function sCurve = findCurveS(PG,sOri,hRim)
% this function finds the first point on the boundary, starting from sOri,
% that is at rim height.
BG = PG.get('1Pos',sOri);
nextvert = find(PG.S(PG.VL)>sOri,1,'first');
nextvertpos = PG.get('1Pos',PG.S(PG.VL(nextvert)));
dir = 1-2*(nextvertpos(2)>=BG(2));
if dir==-1
    range = [nextvert:-1:1,length(PG.VL)-1:-1:nextvert-1];
else
    range = [nextvert-1:1:length(PG.VL)-1,1:1:nextvert];
end


for i = 2:length(range(2:end))
    m = range(i);prevm = range(i-1);
    
    if PG.vertex(2,m)<=hRim
        svert1 = PG.S(PG.VL(prevm));
        if m==1 && dir == 1
            m = length(PG.VL);
        end
        svert2 = PG.S(PG.VL(m));
        pvert1 = PG.get('1Pos',svert1);
        pvert2 = PG.get('1Pos',svert2);
        alpha = (hRim-pvert1(2))/(pvert2(2)-pvert1(2));
        sEnd = svert1 + (svert2-svert1)*alpha;
        sCurve = [sOri,PG.S(PG.VL(range(2:i-1))).',sEnd];
        if sEnd == inf || sEnd == -inf
            disp('inf Error');
        end
        return
    end
end   
end
